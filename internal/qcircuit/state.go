package qcircuit

import (
	"fmt"
	"math"
	"math/rand"
)

// State is a normalised pure state vector.
type State struct {
	n   int
	amp []complex128
}

// NewState returns the |0...0> state of n qubits.
func NewState(n int) *State {
	s := &State{n: n, amp: make([]complex128, 1<<n)}
	s.amp[0] = 1
	return s
}

// FromState validates and adopts an amplitude slice (copied). The slice must
// have length 2^n and be normalised within 1e-9.
func FromState(n int, amp []complex128) (*State, error) {
	if n < 1 || len(amp) != 1<<n {
		return nil, fmt.Errorf("%w: want %d amplitudes, got %d", ErrInvalidState, 1<<n, len(amp))
	}
	var norm float64
	for _, a := range amp {
		norm += real(a)*real(a) + imag(a)*imag(a)
	}
	if math.Abs(norm-1) > 1e-9 {
		return nil, fmt.Errorf("%w: squared norm %.12g", ErrInvalidState, norm)
	}
	out := make([]complex128, len(amp))
	copy(out, amp)
	return &State{n: n, amp: out}, nil
}

// Clone returns an independent copy of the state.
func (s *State) Clone() *State {
	out := &State{n: s.n, amp: make([]complex128, len(s.amp))}
	copy(out.amp, s.amp)
	return out
}

// NumQubits reports the state width.
func (s *State) NumQubits() int { return s.n }

// Amplitudes returns a defensive copy of the amplitude vector.
func (s *State) Amplitudes() []complex128 {
	out := make([]complex128, len(s.amp))
	copy(out, s.amp)
	return out
}

// Amplitude returns the amplitude of computational basis state |index>.
func (s *State) Amplitude(index uint64) (complex128, error) {
	if int(index) >= len(s.amp) {
		return 0, fmt.Errorf("%w: basis index %d for %d qubits", ErrQubitOutOfRange, index, s.n)
	}
	return s.amp[index], nil
}

// Probabilities returns a copy of the computational-basis measurement
// distribution; sum of entries is 1 up to floating-point precision.
func (s *State) Probabilities() []float64 {
	out := make([]float64, len(s.amp))
	var total float64
	for i, a := range s.amp {
		p := real(a)*real(a) + imag(a)*imag(a)
		out[i] = p
		total += p
	}
	// Renormalise tiny drift away so callers always see a proper distribution.
	if total > 0 && math.Abs(total-1) > 0 && math.Abs(total-1) < 1e-12 {
		for i := range out {
			out[i] /= total
		}
	}
	return out
}

// Norm returns the squared L2 norm (should be 1 for physical states).
func (s *State) Norm() float64 {
	var nrm float64
	for _, a := range s.amp {
		nrm += real(a)*real(a) + imag(a)*imag(a)
	}
	return nrm
}

// Run evolves a fresh copy of |0...0> through the circuit and returns the
// final state. The circuit is validated eagerly when gates are added, so Run
// only fails on a width mismatch between state and circuit.
func (c *Circuit) Run() (*State, error) {
	s := NewState(c.n)
	if err := c.RunInto(s); err != nil {
		return nil, err
	}
	return s, nil
}

// RunInto applies the circuit to an existing state, mutating it in place.
func (c *Circuit) RunInto(s *State) error {
	if s.n != c.n {
		return fmt.Errorf("%w: circuit has %d qubits, state has %d", ErrQubitOutOfRange, c.n, s.n)
	}
	for i := range c.gates {
		if err := s.applyGate(c.gates[i]); err != nil {
			return fmt.Errorf("gate %d (%s): %w", i, c.gates[i].Name, err)
		}
	}
	return nil
}

func (s *State) applyGate(g Gate) error {
	switch g.Name {
	case "SWAP":
		s.applySwap(g.Qubits[0], g.Qubits[1], g.Controls)
		return nil
	}
	if len(g.Qubits) != 1 {
		return fmt.Errorf("%w: gate %s expects 1 target", ErrControlConflict, g.Name)
	}
	base := g.Name
	switch base {
	case "CNOT", "CX":
		base = "X"
	case "CY":
		base = "Y"
	case "CZ":
		base = "Z"
	case "CH":
		base = "H"
	case "CRX":
		base = "RX"
	case "CRY":
		base = "RY"
	case "CRZ":
		base = "RZ"
	}
	m := singleMatrix(base, g.Param)
	if m == [4]complex128{} {
		return fmt.Errorf("%w: %s", ErrUnknownGate, g.Name)
	}
	s.applySingle(g.Qubits[0], g.Controls, m)
	return nil
}

func allControlsSet(index uint64, controls []int) bool {
	for _, q := range controls {
		if index&(1<<uint(q)) == 0 {
			return false
		}
	}
	return true
}

// applySingle applies a 2x2 unitary on qubit q when all controls are 1.
func (s *State) applySingle(q int, controls []int, m [4]complex128) {
	bit := uint64(1) << uint(q)
	size := len(s.amp)
	for i := 0; i < size; i++ {
		idx := uint64(i)
		if idx&bit != 0 {
			continue // each pair handled from the zero partner
		}
		j := idx | bit
		if !allControlsSet(idx, controls) {
			continue
		}
		a0 := s.amp[i]
		a1 := s.amp[j]
		s.amp[i] = m[0]*a0 + m[1]*a1
		s.amp[j] = m[2]*a0 + m[3]*a1
	}
}

// applySwap swaps amplitudes for the two target qubits when controls are 1.
func (s *State) applySwap(a, b int, controls []int) {
	ba := uint64(1) << uint(a)
	bb := uint64(1) << uint(b)
	size := len(s.amp)
	for i := 0; i < size; i++ {
		idx := uint64(i)
		if idx&ba != 0 || idx&bb == 0 {
			continue // handle only the a=0,b=1 corner of each quartet
		}
		other := (idx | ba) &^ bb
		if !allControlsSet(idx, controls) || !allControlsSet(other, controls) {
			continue
		}
		s.amp[idx], s.amp[other] = s.amp[other], s.amp[idx]
	}
}

// ProbabilityOf returns the probability that measuring qubit q yields 1.
func (s *State) ProbabilityOf(q int) (float64, error) {
	if q < 0 || q >= s.n {
		return 0, fmt.Errorf("%w: qubit %d for %d-qubit state", ErrQubitOutOfRange, q, s.n)
	}
	bit := uint64(1) << uint(q)
	var p float64
	for i := range s.amp {
		if uint64(i)&bit != 0 {
			a := s.amp[i]
			p += real(a)*real(a) + imag(a)*imag(a)
		}
	}
	return p, nil
}

// Measure samples qubit q in the computational basis using rng (use a fixed
// seed for reproducibility), collapses the state in place, and returns the
// observed bit 0/1.
func (s *State) Measure(q int, rng *rand.Rand) (int, error) {
	p1, err := s.ProbabilityOf(q)
	if err != nil {
		return 0, err
	}
	r := rng.Float64()
	bit := 0
	if r < p1 {
		bit = 1
	}
	if err := s.collapse(q, bit); err != nil {
		return 0, err
	}
	return bit, nil
}

// collapse projects qubit q onto the given outcome and renormalises.
func (s *State) collapse(q, outcome int) error {
	bit := uint64(1) << uint(q)
	norm := math.Sqrt(s.outcomeProb(q, outcome))
	if norm == 0 {
		return fmt.Errorf("qcircuit: cannot collapse qubit %d onto zero-probability outcome %d", q, outcome)
	}
	inv := 1 / norm
	for i := range s.amp {
		set := int((uint64(i)&bit)>>uint(q)) == 1
		if set != (outcome == 1) {
			s.amp[i] = 0
		} else {
			s.amp[i] *= complex(inv, 0)
		}
	}
	return nil
}

func (s *State) outcomeProb(q, outcome int) float64 {
	bit := uint64(1) << uint(q)
	var p float64
	for i := range s.amp {
		set := (uint64(i) & bit) != 0
		if set == (outcome == 1) {
			a := s.amp[i]
			p += real(a)*real(a) + imag(a)*imag(a)
		}
	}
	return p
}

// MeasureAll samples every qubit in the computational basis, collapses the
// state to the observed basis string, and returns the bits indexed by qubit.
func (s *State) MeasureAll(rng *rand.Rand) ([]int, error) {
	out := make([]int, s.n)
	for q := 0; q < s.n; q++ {
		b, err := s.Measure(q, rng)
		if err != nil {
			return nil, err
		}
		out[q] = b
	}
	return out, nil
}

// PauliExpectation returns <psi| P |psi> for a tensor-product Pauli string.
// pauli is indexed by qubit and may contain only 'I','X','Y','Z'.
func (s *State) PauliExpectation(pauli []byte) (float64, error) {
	if len(pauli) != s.n {
		return 0, fmt.Errorf("%w: want %d Pauli letters, got %d", ErrInvalidPauli, s.n, len(pauli))
	}
	var xBits, zBits uint64
	var phase complex128 = 1
	for q, ch := range pauli {
		switch ch {
		case 'I', 'i':
		case 'X', 'x':
			xBits |= 1 << uint(q)
		case 'Y', 'y':
			xBits |= 1 << uint(q)
			zBits |= 1 << uint(q)
			phase *= 1i // Y = i X Z convention: Y|j> contributes i relative to XZ
		case 'Z', 'z':
			zBits |= 1 << uint(q)
		default:
			return 0, fmt.Errorf("%w: %q at position %d", ErrInvalidPauli, ch, q)
		}
	}
	// P|j> = phase * sign(z) * |j xor x>, sign(z) = (-1)^{popcount(j & z)}.
	var sum complex128
	for j := range s.amp {
		idx := uint64(j)
		sign := 1.0
		if bitsOnesCount(idx&zBits)%2 == 1 {
			sign = -1
		}
		k := idx ^ xBits
		sum += complex(sign, 0) * phase * conj(s.amp[j]) * s.amp[k]
	}
	return real(sum), nil
}

func conj(z complex128) complex128 { return complex(real(z), -imag(z)) }

func bitsOnesCount(x uint64) int {
	// math/bits.OnesCount64 is standard library; avoid extra import churn.
	n := 0
	for x != 0 {
		x &= x - 1
		n++
	}
	return n
}
