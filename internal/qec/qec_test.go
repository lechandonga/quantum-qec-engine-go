package qec_test

import (
	"errors"
	"math"
	"testing"

	"github.com/lechandonga/quantum-qec-engine-go/internal/qec"
)

func TestCodeValidation(t *testing.T) {
	for _, d := range []int{1, 2, 4, 0, -3} {
		if _, err := qec.NewCode(d); !errors.Is(err, qec.ErrInvalidDistance) {
			t.Fatalf("distance %d: want ErrInvalidDistance, got %v", d, err)
		}
	}
	c, err := qec.NewCode(3)
	if err != nil {
		t.Fatal(err)
	}
	if c.NumQubits() != 9 || c.CorrectableWeight() != 1 {
		t.Fatalf("d=3: n=%d t=%d, want 9/1", c.NumQubits(), c.CorrectableWeight())
	}

	if _, err := qec.ExactLogicalFailure(3, qec.Rates{PX: -0.1, PZ: 0.01}); !errors.Is(err, qec.ErrInvalidProbability) {
		t.Fatalf("negative p: want ErrInvalidProbability, got %v", err)
	}
	if _, err := qec.Sample(c, qec.Rates{PX: 0.5, PZ: 0.5}, qec.SampleConfig{Shots: 0, Seed: 1}); !errors.Is(err, qec.ErrResourceLimit) {
		t.Fatalf("zero shots: want ErrResourceLimit, got %v", err)
	}
}

func TestSyndromesAndCorrection(t *testing.T) {
	c, _ := qec.NewCode(3)

	t.Run("no errors corrected", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome != qec.Corrected {
			t.Fatalf("no-error outcome = %s, want corrected", r.Outcome)
		}
	})

	t.Run("single bit flip corrected", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		e.X[c.Qubit(0, 1)] = true
		syn, err := c.Extract(e)
		if err != nil {
			t.Fatal(err)
		}
		if !syn.BitSyndrome[0][0] || !syn.BitSyndrome[0][1] {
			t.Fatalf("single interior X must flip both adjacent edges, got %v", syn.BitSyndrome[0])
		}
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome != qec.Corrected || r.BitLogicalFlip {
			t.Fatalf("single bit flip: outcome=%s bitFlip=%v", r.Outcome, r.BitLogicalFlip)
		}
	})

	t.Run("single phase flip corrected", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		e.Z[c.Qubit(1, 2)] = true
		syn, _ := c.Extract(e)
		if !syn.PhaseOuter[0] || !syn.PhaseOuter[1] {
			t.Fatalf("phase fault in block 1 must flip both outer edges, got %v", syn.PhaseOuter)
		}
		if !syn.PhaseInner[1][1] || syn.PhaseInner[1][0] {
			t.Fatalf("phase fault at block1 pos2 must flip inner edge 1, got %v", syn.PhaseInner[1])
		}
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome != qec.Corrected || r.PhaseLogicalFlip {
			t.Fatalf("single phase flip: outcome=%s phaseFlip=%v", r.Outcome, r.PhaseLogicalFlip)
		}
	})

	t.Run("Y fault is both bit and phase", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		q := c.Qubit(2, 0)
		e.X[q] = true
		e.Z[q] = true
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome != qec.Corrected {
			t.Fatalf("single Y fault: outcome=%s", r.Outcome)
		}
	})
}

func TestUncorrectableErrorsFlagged(t *testing.T) {
	c, _ := qec.NewCode(3)

	t.Run("bit flips beyond the correction radius never silently succeed", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		e.X[c.Qubit(0, 0)] = true
		e.X[c.Qubit(0, 1)] = true
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		// Weight 2 == (d+1)/2 is the minimum-weight undetectable logical
		// flip: its syndrome is zero and the minimum-weight explanation of
		// that syndrome is identity, so the decoder does nothing. The
		// important guarantee is that it must never read "corrected".
		if r.Outcome == qec.Corrected {
			t.Fatalf("weight-%d logical-weight error reported corrected", e.Weight())
		}
	})

	t.Run("two bad blocks at phase level never silently succeed", func(t *testing.T) {
		e := qec.NewErrorSet(c)
		e.Z[c.Qubit(0, 0)] = true
		e.Z[c.Qubit(1, 0)] = true
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome == qec.Corrected {
			t.Fatal("weight-2 phase-level logical error reported corrected")
		}
	})

	t.Run("heavy mixed error never silently succeeds", func(t *testing.T) {
		// Put 2 X errors in each of two blocks and a Z fault in two blocks.
		e := qec.NewErrorSet(c)
		e.X[c.Qubit(0, 0)] = true
		e.X[c.Qubit(0, 1)] = true
		e.X[c.Qubit(1, 0)] = true
		e.X[c.Qubit(1, 1)] = true
		e.Z[c.Qubit(0, 2)] = true
		e.Z[c.Qubit(1, 2)] = true
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome == qec.Corrected {
			t.Fatalf("weight-%d error reported as corrected", e.Weight())
		}
	})
}

func TestDecodeDeterminism(t *testing.T) {
	c, _ := qec.NewCode(5)
	e := qec.NewErrorSet(c)
	e.X[c.Qubit(0, 2)] = true
	e.Z[c.Qubit(3, 1)] = true
	first, err := c.Decode(e)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 10; i++ {
		r, err := c.Decode(e)
		if err != nil {
			t.Fatal(err)
		}
		if r.Outcome != first.Outcome {
			t.Fatalf("decode nondeterministic: %s vs %s", r.Outcome, first.Outcome)
		}
	}
}

func TestExactMonotoneLogicalFailure(t *testing.T) {
	rates := qec.Rates{PX: 0.001, PZ: 0.001}
	distances := []int{3, 5, 7, 9}
	prev := math.Inf(1)
	for _, d := range distances {
		f, err := qec.ExactLogicalFailure(d, rates)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("exact logical failure d=%d: %.6e", d, f)
		if !(f < prev) {
			t.Fatalf("logical failure not decreasing: d=%d f=%.3e >= previous %.3e", d, f, prev)
		}
		prev = f
	}
}

func TestSampledTrendAndReproducibility(t *testing.T) {
	rates := qec.Rates{PX: 0.001, PZ: 0.001}
	distances := []int{3, 5}
	const shots = 200000

	run := func(workers int) []float64 {
		ratesOut := make([]float64, len(distances))
		for i, d := range distances {
			c, err := qec.NewCode(d)
			if err != nil {
				t.Fatal(err)
			}
			rep, err := qec.Sample(c, rates, qec.SampleConfig{
				Shots: shots, Seed: 20261006, Workers: workers,
			})
			if err != nil {
				t.Fatal(err)
			}
			ratesOut[i] = rep.LogicalFailureRate
			t.Logf("workers=%d d=%d corrected=%d possible=%d flip=%d logicalFailure=%.5f",
				workers, d, rep.Corrected, rep.PossibleLogicalFlip, rep.LogicalFlip, rep.LogicalFailureRate)
		}
		return ratesOut
	}

	serial := run(1)
	parallel := run(4)
	for i := range serial {
		if serial[i] != parallel[i] {
			t.Fatalf("d=%d: serial rate %.5f != parallel %.5f", distances[i], serial[i], parallel[i])
		}
	}

	// Trend: d=3 -> d=5 must strictly decrease at this error rate.
	for i := 1; i < len(serial); i++ {
		if !(serial[i] < serial[i-1]) {
			t.Fatalf("sampled logical failure not decreasing at step d=%d->%d: %.5f -> %.5f",
				distances[i-1], distances[i], serial[i-1], serial[i])
		}
	}

	// Sampled estimates should track exact values within Monte-Carlo error.
	for i, d := range distances {
		exact, err := qec.ExactLogicalFailure(d, rates)
		if err != nil {
			t.Fatal(err)
		}
		tol := 5.0 * math.Sqrt(exact*(1-exact)/shots)
		if tol < 0.0004 {
			tol = 0.0004
		}
		if math.Abs(serial[i]-exact) > tol {
			t.Fatalf("d=%d sampled %.5f deviates from exact %.5f beyond tol %.5f",
				d, serial[i], exact, tol)
		}
	}
}
