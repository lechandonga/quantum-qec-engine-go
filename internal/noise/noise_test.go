package noise_test

import (
	"errors"
	"math"
	"math/rand"
	"testing"

	"github.com/lechandonga/quantum-qec-engine-go/internal/noise"
	"github.com/lechandonga/quantum-qec-engine-go/internal/qcircuit"
)

func bellCircuit(t *testing.T) *qcircuit.Circuit {
	t.Helper()
	c, err := qcircuit.NewCircuit(2)
	if err != nil {
		t.Fatal(err)
	}
	if err := c.H(0); err != nil {
		t.Fatal(err)
	}
	if err := c.CNOT(0, 1); err != nil {
		t.Fatal(err)
	}
	return c
}

func TestInvalidNoiseInputs(t *testing.T) {
	bad := []noise.Model{
		{BitFlip: -0.1},
		{BitFlip: 1.1},
		{PhaseFlip: math.NaN()},
		{Depolarising: 2},
		{Qubits: []int{7}},
	}
	for i, m := range bad {
		err := m.Validate(2)
		if err == nil {
			t.Fatalf("model %d should be rejected", i)
		}
		if i < 4 && !errors.Is(err, noise.ErrInvalidStrength) {
			t.Fatalf("model %d: want ErrInvalidStrength, got %v", i, err)
		}
		if i == 4 && !errors.Is(err, noise.ErrInvalidQubit) {
			t.Fatalf("model %d: want ErrInvalidQubit, got %v", i, err)
		}
	}

	c := bellCircuit(t)
	if _, err := noise.RunBatch(noise.BatchRequest{
		Circuit: c, Shots: 10, Seed: 1,
		Limits: noise.BatchLimits{MaxShots: 5},
	}); !errors.Is(err, noise.ErrResourceLimit) {
		t.Fatalf("over-limit shots: want ErrResourceLimit, got %v", err)
	}
	if _, err := noise.RunBatch(noise.BatchRequest{
		Circuit: c, Shots: 0, Seed: 1,
	}); !errors.Is(err, noise.ErrResourceLimit) {
		t.Fatalf("zero shots: want ErrResourceLimit, got %v", err)
	}
}

func TestZeroNoiseIsIdentity(t *testing.T) {
	c := bellCircuit(t)
	ideal, err := c.Run()
	if err != nil {
		t.Fatal(err)
	}

	st, err := c.Run()
	if err != nil {
		t.Fatal(err)
	}
	rng := rand.New(rand.NewSource(99))
	for i := 0; i < 1000; i++ {
		if err := (noise.Model{}).Apply(st, rng); err != nil {
			t.Fatal(err)
		}
	}
	got := st.Amplitudes()
	want := ideal.Amplitudes()
	for i := range got {
		if got[i] != want[i] {
			t.Fatalf("zero noise drifted amplitude[%d]: %v vs %v", i, got[i], want[i])
		}
	}
}

func TestExactAgreementSmallSystem(t *testing.T) {
	// Bit-flip channel p=0.2 on |0> gives P(1)=0.2 exactly.
	d, err := noise.FromPure([]complex128{1, 0})
	if err != nil {
		t.Fatal(err)
	}
	m := noise.Model{BitFlip: 0.2}
	if err := d.ApplyExact(m); err != nil {
		t.Fatal(err)
	}
	if tr := real(d.Trace()); math.Abs(tr-1) > 1e-12 {
		t.Fatalf("trace drift: Tr(rho)=%.12f", tr)
	}
	probs := d.Probabilities()
	if math.Abs(probs[1]-0.2) > 1e-12 {
		t.Fatalf("exact bit flip P(1)=%.4f want 0.2", probs[1])
	}

	// Depolarising p on |0>: <Z> = 1 - 4p/3.
	d2, _ := noise.FromPure([]complex128{1, 0})
	if err := d2.ApplyExact(noise.Model{Depolarising: 0.3}); err != nil {
		t.Fatal(err)
	}
	z, err := d2.PauliExpectation([]byte("Z"))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(z-(1-4*0.3/3)) > 1e-12 {
		t.Fatalf("depolarising <Z>=%.4f want %.4f", z, 1-4*0.3/3)
	}

	// Phase flip keeps diagonal probabilities but kills coherence:
	// H|0> then phase p gives <X>=1-2p.
	d3, _ := noise.FromPure([]complex128{1 / math.Sqrt2, 1 / math.Sqrt2})
	if err := d3.ApplyExact(noise.Model{PhaseFlip: 0.25}); err != nil {
		t.Fatal(err)
	}
	x, err := d3.PauliExpectation([]byte("X"))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(x-0.5) > 1e-12 {
		t.Fatalf("phase flip <X>=%.4f want 0.5", x)
	}

	// Exact density reference refuses large systems up front.
	if _, err := noise.FromPure(make([]complex128, 1<<8)); !errors.Is(err, noise.ErrResourceLimit) {
		t.Fatalf("large exact reference: want ErrResourceLimit, got %v", err)
	}
}

func TestSamplingConvergence(t *testing.T) {
	// Prepare |0> then sample bit-flip p=0.2 directly (identity circuit).
	c, err := qcircuit.NewCircuit(1)
	if err != nil {
		t.Fatal(err)
	}
	res, err := noise.RunBatch(noise.BatchRequest{
		Circuit: c,
		Model:   noise.Model{BitFlip: 0.2},
		Shots:   20000,
		Seed:    7,
		Workers: 1,
	})
	if err != nil {
		t.Fatal(err)
	}
	p := res.EmpiricalProbability(1)
	if math.Abs(p-0.2) > 0.02 {
		t.Fatalf("sampled bit-flip frequency %.4f not within 0.02 of 0.2", p)
	}

	// Bell circuit without noise: only outcomes 00 and 11, each 1/2.
	bell := bellCircuit(t)
	res2, err := noise.RunBatch(noise.BatchRequest{
		Circuit: bell, Model: noise.Model{}, Shots: 8000, Seed: 11, Workers: 4,
	})
	if err != nil {
		t.Fatal(err)
	}
	if got := res2.EmpiricalProbability(1) + res2.EmpiricalProbability(2); got != 0 {
		t.Fatalf("noiseless Bell produced anti-correlated outcomes, frequency %.4f", got)
	}
	if got := res2.EmpiricalProbability(0); math.Abs(got-0.5) > 0.03 {
		t.Fatalf("Bell P(00)=%.4f want ~0.5", got)
	}
}

func TestParallelDeterminism(t *testing.T) {
	c := bellCircuit(t)
	req := func(workers int) noise.BatchRequest {
		return noise.BatchRequest{
			Circuit: c,
			Model:   noise.Model{BitFlip: 0.03, PhaseFlip: 0.02, Depolarising: 0.01},
			Shots:   3000,
			Seed:    12345,
			Workers: workers,
		}
	}
	r1, err := noise.RunBatch(req(1))
	if err != nil {
		t.Fatal(err)
	}
	for _, w := range []int{2, 3, 8} {
		r, err := noise.RunBatch(req(w))
		if err != nil {
			t.Fatal(err)
		}
		for i := range r1.Shots {
			if r1.Shots[i].Index != r.Shots[i].Index {
				t.Fatalf("worker=%d shot %d: %d vs serial %d", w, i,
					r.Shots[i].Index, r1.Shots[i].Index)
			}
		}
	}
}
