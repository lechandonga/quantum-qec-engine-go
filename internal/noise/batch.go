package noise

import (
	"fmt"
	"math/rand"
	"runtime"
	"sync"

	"github.com/lechandonga/quantum-qec-engine-go/internal/qcircuit"
)

// BatchLimits caps batch sampling resources. Zero MaxShots means 10_000_000.
type BatchLimits struct {
	MaxShots int
}

// DefaultBatchLimits is used by RunBatch when no override is supplied.
var DefaultBatchLimits = BatchLimits{MaxShots: 10_000_000}

// BatchRequest configures one batch-sampling run. Circuit is the ideal
// circuit; every shot starts from |0...0>, evolves through it, then receives
// one realisation of Model. MeasureQubits selects which qubits are read
// (nil means all). Workers controls parallelism (0 => GOMAXPROCS, 1 => serial).
type BatchRequest struct {
	Circuit       *qcircuit.Circuit
	Model         Model
	Shots         int
	Workers       int
	Seed          int64
	MeasureQubits []int
	Limits        BatchLimits
}

// ShotResult is one sample: bits indexed by the requested MeasureQubits list
// and Index as the integer with bits[0] the least-significant position.
type ShotResult struct {
	Bits  []int
	Index uint64
}

// BatchResult is the full ordered set of shots plus a convenience histogram.
type BatchResult struct {
	Shots     []ShotResult
	Histogram map[uint64]int
}

// EmpiricalProbability returns the observed frequency of an outcome index.
func (r *BatchResult) EmpiricalProbability(index uint64) float64 {
	if len(r.Shots) == 0 {
		return 0
	}
	return float64(r.Histogram[index]) / float64(len(r.Shots))
}

// streamSeed derives a deterministic per-shot seed from the batch seed and
// the shot index using SplitMix64, so shot i always consumes the same
// pseudo-random stream regardless of scheduling or worker count.
func streamSeed(root int64, index int) int64 {
	var x uint64 = uint64(root) ^ (uint64(index)*0x9E3779B97F4A7C15 + 0xD1B54A32D192ED03)
	x ^= x >> 30
	x *= 0xBF58476D1CE4E5B9
	x ^= x >> 27
	x *= 0x94D049BB133111EB
	x ^= x >> 31
	return int64(x)
}

// RunBatch validates everything up front (rejecting over-limit requests
// before allocating), then executes shots, in parallel when Workers > 1.
// Results are ordered by shot index and are identical for equal Seeds no
// matter how many workers are used.
func RunBatch(req BatchRequest) (*BatchResult, error) {
	if req.Circuit == nil {
		return nil, fmt.Errorf("noise: nil circuit")
	}
	n := req.Circuit.NumQubits()
	if req.Shots <= 0 {
		return nil, fmt.Errorf("%w: shots must be positive, got %d", ErrResourceLimit, req.Shots)
	}
	limits := req.Limits
	if limits.MaxShots == 0 {
		limits = DefaultBatchLimits
	}
	if limits.MaxShots > 0 && req.Shots > limits.MaxShots {
		return nil, fmt.Errorf("%w: %d shots > max %d", ErrResourceLimit, req.Shots, limits.MaxShots)
	}
	if err := req.Model.Validate(n); err != nil {
		return nil, err
	}
	measure := req.MeasureQubits
	if measure == nil {
		measure = make([]int, n)
		for i := range measure {
			measure[i] = i
		}
	}
	for _, q := range measure {
		if q < 0 || q >= n {
			return nil, fmt.Errorf("%w: measure qubit %d for %d-qubit system", ErrInvalidQubit, q, n)
		}
	}

	workers := req.Workers
	if workers <= 0 {
		workers = runtime.GOMAXPROCS(0)
	}
	if workers > req.Shots {
		workers = req.Shots
	}

	results := make([]ShotResult, req.Shots)
	if workers == 1 {
		for i := 0; i < req.Shots; i++ {
			res, err := runOneShot(req.Circuit, req.Model, measure, req.Seed, i)
			if err != nil {
				return nil, err
			}
			results[i] = res
		}
	} else {
		jobs := make(chan int)
		var wg sync.WaitGroup
		var firstErr error
		var errOnce sync.Once
		for w := 0; w < workers; w++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for i := range jobs {
					res, err := runOneShot(req.Circuit, req.Model, measure, req.Seed, i)
					if err != nil {
						errOnce.Do(func() { firstErr = err })
						return
					}
					results[i] = res
				}
			}()
		}
		for i := 0; i < req.Shots; i++ {
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

	hist := make(map[uint64]int, 256)
	for i := range results {
		hist[results[i].Index]++
	}
	return &BatchResult{Shots: results, Histogram: hist}, nil
}

func runOneShot(c *qcircuit.Circuit, m Model, measure []int, seed int64, shotIndex int) (ShotResult, error) {
	rng := rand.New(rand.NewSource(streamSeed(seed, shotIndex)))
	st, err := c.Run()
	if err != nil {
		return ShotResult{}, err
	}
	if !m.IsNoNoise() {
		if err := m.Apply(st, rng); err != nil {
			return ShotResult{}, err
		}
	}
	bits := make([]int, len(measure))
	// Collapse in ascending qubit order so the draw sequence is fixed.
	for k, q := range measure {
		b, err := st.Measure(q, rng)
		if err != nil {
			return ShotResult{}, err
		}
		bits[k] = b
	}
	var idx uint64
	for k, b := range bits {
		if b == 1 {
			idx |= 1 << uint(k)
		}
	}
	return ShotResult{Bits: bits, Index: idx}, nil
}
