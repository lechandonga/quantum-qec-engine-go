package qec

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sync"
)

// Rates configures the independent symmetric-ish per-qubit Pauli fault model
// used for statistics: bit flip probability pX, phase flip probability pZ.
// The two draws are independent per qubit, so a qubit may receive both
// (a Y fault).
type Rates struct {
	PX float64
	PZ float64
}

func (r Rates) validate() error {
	if r.PX < 0 || r.PX > 1 || r.PZ < 0 || r.PZ > 1 {
		return fmt.Errorf("%w: pX=%v pZ=%v", ErrInvalidProbability, r.PX, r.PZ)
	}
	return nil
}

// SampleConfig controls one Monte-Carlo experiment. Workers <= 1 runs
// serially; 0 means GOMAXPROCS.
type SampleConfig struct {
	Shots   int
	Workers int
	Seed    int64
}

// SampleReport aggregates decoder outcomes for a batch of fault patterns.
type SampleReport struct {
	Shots               int
	Corrected           int
	PossibleLogicalFlip int
	LogicalFlip         int
	LogicalFailureRate  float64 // (LogicalFlip + PossibleLogicalFlip) / Shots
	DefiniteFlipRate    float64 // LogicalFlip / Shots
}

// Sample draws fault patterns from the independent per-qubit model, decodes
// each one, and aggregates outcomes. A given (code, rates, seed, shots) tuple
// always produces the same report regardless of worker count: each shot uses
// a deterministic child stream keyed by its index, and writes into a
// pre-sized slice at its own position.
func Sample(c *Code, rates Rates, cfg SampleConfig) (*SampleReport, error) {
	if err := rates.validate(); err != nil {
		return nil, err
	}
	if cfg.Shots <= 0 {
		return nil, fmt.Errorf("%w: shots must be positive, got %d", ErrResourceLimit, cfg.Shots)
	}
	workers := cfg.Workers
	if workers <= 0 {
		workers = runtime.GOMAXPROCS(0)
	}
	if workers > cfg.Shots {
		workers = cfg.Shots
	}

	outcomes := make([]DecodeOutcome, cfg.Shots)
	if workers == 1 {
		for i := 0; i < cfg.Shots; i++ {
			o, err := sampleOne(c, rates, childRand(cfg.Seed, i))
			if err != nil {
				return nil, err
			}
			outcomes[i] = o
		}
	} else {
		jobs := make(chan int)
		var wg sync.WaitGroup
		var firstErr error
		var once sync.Once
		for w := 0; w < workers; w++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for i := range jobs {
					o, err := sampleOne(c, rates, childRand(cfg.Seed, i))
					if err != nil {
						once.Do(func() { firstErr = err })
						return
					}
					outcomes[i] = o
				}
			}()
		}
		for i := 0; i < cfg.Shots; i++ {
			if firstErr != nil {
				break
			}
			jobs <- i
		}
		close(jobs)
		wg.Wait()
		if firstErr != nil {
			return nil, firstErr
		}
	}

	rep := &SampleReport{Shots: cfg.Shots}
	for _, o := range outcomes {
		switch o {
		case Corrected:
			rep.Corrected++
		case PossibleLogicalFlip:
			rep.PossibleLogicalFlip++
		case LogicalFlip:
			rep.LogicalFlip++
		}
	}
	rep.LogicalFailureRate = float64(rep.LogicalFlip+rep.PossibleLogicalFlip) / float64(rep.Shots)
	rep.DefiniteFlipRate = float64(rep.LogicalFlip) / float64(rep.Shots)
	return rep, nil
}

// childRand derives a per-shot random source from (seed, shotIndex). It is a
// pure function of its inputs, so scheduling cannot perturb results.
func childRand(seed int64, index int) *rand.Rand {
	var x uint64 = uint64(seed) ^ (uint64(index)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03)
	x ^= x >> 30
	x *= 0xBF58476D1CE4E5B9
	x ^= x >> 27
	x *= 0x94D049BB133111EB
	x ^= x >> 31
	return rand.New(rand.NewSource(int64(x)))
}

func sampleOne(c *Code, rates Rates, rng *rand.Rand) (DecodeOutcome, error) {
	e := NewErrorSet(c)
	// Fixed draw order: qubit q: bit draw, then phase draw.
	for q := 0; q < c.n; q++ {
		if rng.Float64() < rates.PX {
			e.X[q] = true
		}
		if rng.Float64() < rates.PZ {
			e.Z[q] = true
		}
	}
	res, err := c.Decode(e)
	if err != nil {
		return 0, err
	}
	return res.Outcome, nil
}

// ExactLogicalFailure computes the decoder-failure probability without
// Monte-Carlo error, using the product structure of the code: bit-flip
// failure in any one block plus majority failure at the phase level.
//
//	bit level  : repetition block length d, per-qubit error pX
//	phase level: block-parity error probability q=(1-(1-2pZ)^d)/2, then
//	             majority over d blocks
func ExactLogicalFailure(d int, rates Rates) (float64, error) {
	c, err := NewCode(d)
	if err != nil {
		return 0, err
	}
	if err := rates.validate(); err != nil {
		return 0, err
	}
	t := c.CorrectableWeight()

	// Each of the d blocks is an independent length-d repetition code: a
	// block corrects iff it contains at most t bit faults.
	blockCorrect := binomialAtMost(d, rates.PX, t)
	bitLevel := 1 - math.Pow(blockCorrect, float64(d))

	q := (1 - math.Pow(1-2*rates.PZ, float64(d))) / 2
	fPhase := majorityFailure(d, q, t)

	return bitLevel + fPhase - bitLevel*fPhase, nil
}

// binomialAtMost returns P(Binomial(n,p) <= t).
func binomialAtMost(n int, p float64, t int) float64 {
	var s float64
	for k := 0; k <= t; k++ {
		logC := logBinom(n, k)
		s += math.Exp(logC + float64(k)*math.Log(p) + float64(n-k)*math.Log(1-p))
	}
	return s
}

// majorityFailure returns P(Binomial(n,p) > t), the failure probability of
// majority decoding with guaranteed correction up to t errors.
func majorityFailure(n int, p float64, t int) float64 {
	if p == 0 {
		return 0
	}
	if p == 1 {
		return 1
	}
	var f float64
	for k := t + 1; k <= n; k++ {
		logC := logBinom(n, k)
		logP := float64(k) * math.Log(p)
		logQ := float64(n-k) * math.Log(1-p)
		f += math.Exp(logC + logP + logQ)
	}
	return f
}

func logBinom(n, k int) float64 {
	if k < 0 || k > n {
		return math.Inf(-1)
	}
	return lgamma(float64(n)+1) - lgamma(float64(k)+1) - lgamma(float64(n-k)+1)
}

// lgamma wraps math.Lgamma (returns log|Gamma|; Gamma positive for x>0).
func lgamma(x float64) float64 {
	v, _ := math.Lgamma(x)
	return v
}
