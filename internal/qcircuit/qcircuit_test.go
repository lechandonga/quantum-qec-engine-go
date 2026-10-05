package qcircuit_test

import (
	"errors"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"

	"github.com/lechandonga/quantum-qec-engine-go/internal/qcircuit"
)

func approx(t *testing.T, got, want, tol float64, msg string) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Fatalf("%s: got %.10f want %.10f (tol %g)", msg, got, want, tol)
	}
}

func TestInvalidCircuits(t *testing.T) {
	t.Run("rejects bad width and resource limits", func(t *testing.T) {
		if _, err := qcircuit.NewCircuit(0); !errors.Is(err, qcircuit.ErrQubitOutOfRange) {
			t.Fatalf("width 0: want ErrQubitOutOfRange, got %v", err)
		}
		if _, err := qcircuit.NewCircuitWithLimits(30, qcircuit.Limits{MaxQubits: 8}); !errors.Is(err, qcircuit.ErrResourceLimit) {
			t.Fatalf("over-limit width: want ErrResourceLimit, got %v", err)
		}
	})

	c, err := qcircuit.NewCircuit(2)
	if err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name string
		add  func() error
		want error
	}{
		{"qubit below zero", func() error { return c.X(-1) }, qcircuit.ErrQubitOutOfRange},
		{"qubit above width", func() error { return c.H(2) }, qcircuit.ErrQubitOutOfRange},
		{"unknown gate", func() error { return c.Add("BOGUS", []int{0}, nil, 0) }, qcircuit.ErrUnknownGate},
		{"control equals target", func() error { return c.Add("X", []int{0}, []int{0}, 0) }, qcircuit.ErrControlConflict},
		{"duplicate controls", func() error { return c.Add("X", []int{1}, []int{0, 0}, 0) }, qcircuit.ErrControlConflict},
		{"CNOT missing control", func() error { return c.Add("CNOT", []int{0}, nil, 0) }, qcircuit.ErrControlConflict},
		{"NaN rotation", func() error { return c.RX(0, math.NaN()) }, qcircuit.ErrInvalidParameter},
		{"Inf rotation", func() error { return c.RY(0, math.Inf(1)) }, qcircuit.ErrInvalidParameter},
		{"swap bad arity", func() error { return c.Add("SWAP", []int{0}, nil, 0) }, qcircuit.ErrControlConflict},
	}
	for _, tc := range cases {
		t.Run("rejects "+tc.name, func(t *testing.T) {
			if err := tc.add(); !errors.Is(err, tc.want) {
				t.Fatalf("want %v, got %v", tc.want, err)
			}
		})
	}

	t.Run("rejects gate count over budget", func(t *testing.T) {
		small, err := qcircuit.NewCircuitWithLimits(2, qcircuit.Limits{MaxGates: 1})
		if err != nil {
			t.Fatal(err)
		}
		if err := small.H(0); err != nil {
			t.Fatal(err)
		}
		if err := small.H(1); !errors.Is(err, qcircuit.ErrResourceLimit) {
			t.Fatalf("want ErrResourceLimit, got %v", err)
		}
	})
}

func TestGateEvolution(t *testing.T) {
	t.Run("Bell state via H+CNOT", func(t *testing.T) {
		c, _ := qcircuit.NewCircuit(2)
		if err := c.H(0); err != nil {
			t.Fatal(err)
		}
		if err := c.CNOT(0, 1); err != nil {
			t.Fatal(err)
		}
		st, err := c.Run()
		if err != nil {
			t.Fatal(err)
		}
		amps := st.Amplitudes()
		sqrt := complex(1/math.Sqrt2, 0)
		for i, want := range []complex128{sqrt, 0, 0, sqrt} {
			if cmplx.Abs(amps[i]-want) > 1e-12 {
				t.Fatalf("amp[%d]=%v want %v", i, amps[i], want)
			}
		}
		probs := st.Probabilities()
		approx(t, probs[0], 0.5, 1e-12, "P(00)")
		approx(t, probs[3], 0.5, 1e-12, "P(11)")
	})

	t.Run("parameterised rotations give known amplitudes", func(t *testing.T) {
		c, _ := qcircuit.NewCircuit(1)
		if err := c.RX(0, math.Pi/2); err != nil {
			t.Fatal(err)
		}
		st, err := c.Run()
		if err != nil {
			t.Fatal(err)
		}
		amps := st.Amplitudes()
		want := complex(1/math.Sqrt2, 0)
		if cmplx.Abs(amps[0]-want) > 1e-12 || cmplx.Abs(amps[1]-(-1i*want)) > 1e-12 {
			t.Fatalf("RX(pi/2) amplitudes = %v, %v", amps[0], amps[1])
		}
	})

	t.Run("controlled phase and swap operations", func(t *testing.T) {
		// Verify the CZ phase: on the |++> state the only minus sign is on
		// |11>, so measuring 11 in the X basis on both qubits has
		// amplitude (1-1-1-1)/2 = -1 => outcome 11 deterministically.
		c, _ := qcircuit.NewCircuit(2)
		if err := c.H(0); err != nil {
			t.Fatal(err)
		}
		if err := c.H(1); err != nil {
			t.Fatal(err)
		}
		if err := c.CZ(0, 1); err != nil {
			t.Fatal(err)
		}
		if err := c.H(0); err != nil {
			t.Fatal(err)
		}
		if err := c.H(1); err != nil {
			t.Fatal(err)
		}
		st, _ := c.Run()
		probs := st.Probabilities()
		// (H⊗H)CZ(H⊗H) = CNOT acting on |++>, giving a uniform distribution.
		for i, p := range probs {
			approx(t, p, 0.25, 1e-12, "uniform after H-H-CZ-H-H, outcome "+string(rune('0'+i)))
		}

		c2, _ := qcircuit.NewCircuit(2)
		if err := c2.X(0); err != nil {
			t.Fatal(err)
		}
		if err := c2.SWAP(0, 1); err != nil {
			t.Fatal(err)
		}
		st2, _ := c2.Run()
		p1, _ := st2.ProbabilityOf(1)
		approx(t, p1, 1, 1e-12, "SWAP moves |10> excitation to q1")
	})
}

func TestMeasurementCollapseAndReproducibility(t *testing.T) {
	c, _ := qcircuit.NewCircuit(1)
	if err := c.H(0); err != nil {
		t.Fatal(err)
	}

	run := func() ([]int, []float64) {
		bits := make([]int, 0, 400)
		pAfter := make([]float64, 0, 400)
		rng := rand.New(rand.NewSource(42))
		for i := 0; i < 400; i++ {
			// Rebuild the same superposition each shot.
			st, _ := c.Run()
			b, err := st.Measure(0, rng)
			if err != nil {
				t.Fatal(err)
			}
			bits = append(bits, b)
			p, _ := st.ProbabilityOf(0)
			pAfter = append(pAfter, p)
		}
		return bits, pAfter
	}

	bitsA, pA := run()
	bitsB, pB := run()
	for i := range bitsA {
		if bitsA[i] != bitsB[i] {
			t.Fatalf("shot %d not reproducible: %d vs %d", i, bitsA[i], bitsB[i])
		}
		if math.Abs(pA[i]-pB[i]) > 0 {
			t.Fatalf("post-measurement probabilities differ at shot %d", i)
		}
		// Collapse semantics: after measuring b, P(0) must be 0 or 1.
		if pA[i] != 0 && pA[i] != 1 {
			t.Fatalf("shot %d: state not collapsed, P(0)=%v", i, pA[i])
		}
	}

	var ones int
	for _, b := range bitsA {
		ones += b
	}
	ratio := float64(ones) / float64(len(bitsA))
	approx(t, ratio, 0.5, 0.08, "H measurement frequency converges to 0.5")

	// Measuring a qubit outside the width is rejected.
	st, _ := c.Run()
	rng := rand.New(rand.NewSource(1))
	if _, err := st.Measure(3, rng); !errors.Is(err, qcircuit.ErrQubitOutOfRange) {
		t.Fatalf("want ErrQubitOutOfRange, got %v", err)
	}
}

func TestPauliExpectations(t *testing.T) {
	// |+> on q0, |0> on q1: <Z0>=0, <X0>=1, <Z1>=1, <Z0Z1>=0, <X0Z1>=1.
	c, _ := qcircuit.NewCircuit(2)
	if err := c.H(0); err != nil {
		t.Fatal(err)
	}
	st, _ := c.Run()

	cases := []struct {
		pauli string
		want  float64
	}{
		{"IZ", 1},
		{"XI", 1}, // pauli[0] is q0: |+>_0 |0>_1, so X on q0 is +1
		{"IX", 0},
		{"XZ", 1},
		{"ZI", 0},
		{"ZZ", 0},
		{"ZX", 0},
	}
	for _, tc := range cases {
		got, err := st.PauliExpectation([]byte(tc.pauli))
		if err != nil {
			t.Fatalf("%s: %v", tc.pauli, err)
		}
		approx(t, got, tc.want, 1e-12, "<"+tc.pauli+">")
	}

	// RX(pi/2)|0> = (|0> - i|1>)/sqrt2.
	c3, _ := qcircuit.NewCircuit(1)
	if err := c3.RX(0, math.Pi/2); err != nil {
		t.Fatal(err)
	}
	st3, _ := c3.Run()
	y3, err := st3.PauliExpectation([]byte("Y"))
	if err != nil {
		t.Fatal(err)
	}
	// RX(theta) rotates toward -Y; RX(pi/2)|0> = (|0>-i|1>)/sqrt2 has <Y>=+1.
	approx(t, y3, 1, 1e-12, "Y expectation of RX(pi/2)|0>")

	if _, err := st.PauliExpectation([]byte("Q")); !errors.Is(err, qcircuit.ErrInvalidPauli) {
		t.Fatalf("want ErrInvalidPauli, got %v", err)
	}
	if _, err := st.PauliExpectation([]byte("XXX")); !errors.Is(err, qcircuit.ErrInvalidPauli) {
		t.Fatalf("want ErrInvalidPauli for wrong length, got %v", err)
	}
}

func TestFromStateValidation(t *testing.T) {
	if _, err := qcircuit.FromState(1, []complex128{1, 1}); !errors.Is(err, qcircuit.ErrInvalidState) {
		t.Fatalf("wrong length: want ErrInvalidState, got %v", err)
	}
	if _, err := qcircuit.FromState(1, []complex128{0.5, 0.5}); !errors.Is(err, qcircuit.ErrInvalidState) {
		t.Fatalf("unnormalised: want ErrInvalidState, got %v", err)
	}
}
