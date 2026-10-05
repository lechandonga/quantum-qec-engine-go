// Package noise implements stochastic Pauli noise channels, a small
// density-matrix reference implementation for exact checks, and a
// deterministic parallel batch sampler.
//
// All randomness flows through *math/rand.Rand values derived from a single
// caller-supplied seed, so results are reproducible bit-for-bit regardless of
// the worker count used for batch sampling.
package noise

import (
	"errors"
	"fmt"
	"math"
	"math/rand"

	"github.com/lechandonga/quantum-qec-engine-go/internal/qcircuit"
)

var (
	// ErrInvalidStrength is returned when a noise strength is outside [0,1].
	ErrInvalidStrength = errors.New("noise: invalid noise strength")
	// ErrResourceLimit is returned before any work starts when a request
	// exceeds the configured budget.
	ErrResourceLimit = errors.New("noise: resource limit exceeded")
	// ErrInvalidQubit is returned when the model names an unknown qubit.
	ErrInvalidQubit = errors.New("noise: qubit index out of range")
)

// Model describes independent per-qubit Pauli channels. The three channels
// compose sequentially per qubit (bit flip, then phase flip, then
// depolarising); at zero strengths the composition is the identity channel.
//
// Strengths are Pauli-error probabilities:
//   - BitFlip p:       rho -> (1-p) rho + p X rho X.
//   - PhaseFlip p:     rho -> (1-p) rho + p Z rho Z.
//   - Depolarising p:  rho -> (1-p) rho + p/3 (X rho X + Y rho Y + Z rho Z).
//
// Qubits selects the qubits the channel acts on; nil means every qubit.
type Model struct {
	BitFlip      float64
	PhaseFlip    float64
	Depolarising float64
	Qubits       []int
}

// Validate checks probability ranges and qubit indices for an n-qubit system.
func (m Model) Validate(n int) error {
	for name, p := range map[string]float64{
		"bit-flip": m.BitFlip, "phase-flip": m.PhaseFlip, "depolarising": m.Depolarising,
	} {
		if math.IsNaN(p) || p < 0 || p > 1 {
			return fmt.Errorf("%w: %s strength %v", ErrInvalidStrength, name, p)
		}
	}
	for _, q := range m.Qubits {
		if q < 0 || q >= n {
			return fmt.Errorf("%w: qubit %d for %d-qubit system", ErrInvalidQubit, q, n)
		}
	}
	return nil
}

// targetQubits resolves the qubit list for an n-qubit system.
func (m Model) targetQubits(n int) []int {
	if m.Qubits != nil {
		return m.Qubits
	}
	out := make([]int, n)
	for i := range out {
		out[i] = i
	}
	return out
}

// IsNoNoise reports whether the model is the identity channel.
func (m Model) IsNoNoise() bool {
	return m.BitFlip == 0 && m.PhaseFlip == 0 && m.Depolarising == 0
}

// SampleError draws one Pauli fault per affected qubit from the model using
// rng in fixed order (bit flip draw, phase flip draw, depolarising draw for
// qubit 0, then qubit 1, ...). The returned circuit contains the sampled X/Y/Z
// gates and may be applied to a state.
func (m Model) SampleError(n int, rng *rand.Rand) (*qcircuit.Circuit, error) {
	if err := m.Validate(n); err != nil {
		return nil, err
	}
	c, err := qcircuit.NewCircuit(n)
	if err != nil {
		return nil, err
	}
	for _, q := range m.targetQubits(n) {
		if rng.Float64() < m.BitFlip {
			if err := c.X(q); err != nil {
				return nil, err
			}
		}
		if rng.Float64() < m.PhaseFlip {
			if err := c.Z(q); err != nil {
				return nil, err
			}
		}
		if rng.Float64() < m.Depolarising {
			switch rng.Intn(3) {
			case 0:
				err = c.X(q)
			case 1:
				err = c.Y(q)
			case 2:
				err = c.Z(q)
			}
			if err != nil {
				return nil, err
			}
		}
	}
	return c, nil
}

// Apply mutates st by appending one realisation of the noise model.
func (m Model) Apply(st *qcircuit.State, rng *rand.Rand) error {
	errC, err := m.SampleError(st.NumQubits(), rng)
	if err != nil {
		return err
	}
	return errC.RunInto(st)
}
