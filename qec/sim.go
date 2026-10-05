package qec

import (
	"fmt"
	"math/rand/v2"
	"runtime"
	"sync"
)

// Limits 是 QEC 仿真的资源上限。
type Limits struct {
	MaxDistance int
	MaxShots    int
	MaxWorkers  int
}

// DefaultLimits 为默认上限。
var DefaultLimits = Limits{
	MaxDistance: 7, // d=7 => 49 物理比特（受 PauliError 位掩码约束）
	MaxShots:    1_000_000,
	MaxWorkers:  64,
}

// Simulator 驱动固定错误率下的批量纠错仿真。
type Simulator struct {
	Limits Limits
}

// NewSimulator 创建仿真器，零值字段回填默认上限。
func NewSimulator(l Limits) *Simulator {
	if l.MaxDistance <= 0 {
		l.MaxDistance = DefaultLimits.MaxDistance
	}
	if l.MaxShots <= 0 {
		l.MaxShots = DefaultLimits.MaxShots
	}
	if l.MaxWorkers <= 0 {
		l.MaxWorkers = DefaultLimits.MaxWorkers
	}
	return &Simulator{Limits: l}
}

// SampleError 使用给定随机源，按对称 Pauli 噪声（每比特每种错误 p/3）
// 抽样一个多比特错误。
func (c Code) SampleError(p float64, rng *rand.Rand) (PauliError, error) {
	if p < 0 || p > 1 {
		return PauliError{}, ErrInvalidErrorRate
	}
	if rng == nil {
		return PauliError{}, fmt.Errorf("qec: deterministic random source is required")
	}
	e := PauliError{}
	third := p / 3
	for q := 0; q < c.N; q++ {
		r := rng.Float64()
		switch {
		case r < third:
			e.Xmask |= 1 << q // X
		case r < 2*third:
			e.Xmask |= 1 << q // Y = XZ
			e.Zmask |= 1 << q
		case r < p:
			e.Zmask |= 1 << q // Z
		}
	}
	return e, nil
}

// TrialResult 汇总一次完整「抽样错误 → 解码」。
type TrialResult struct {
	Error   PauliError
	Outcome Outcome
	Weight  int
}

// Run 执行 shots 次纠错试验并统计。结果只取决于 (seed, 试验编号)，
// 与并行度无关；workers<=0 时取 runtime.NumCPU。
func (s *Simulator) Run(distance int, p float64, shots int, seed uint64, workers int) ([]TrialResult, Stats, error) {
	code, err := NewShorCode(distance)
	if err != nil {
		return nil, Stats{}, err
	}
	if p < 0 || p > 1 {
		return nil, Stats{}, fmt.Errorf("%w: p=%g", ErrInvalidErrorRate, p)
	}
	if shots < 0 || shots > s.Limits.MaxShots {
		return nil, Stats{}, fmt.Errorf("%w: shots=%d > max=%d", ErrResourceLimit, shots, s.Limits.MaxShots)
	}
	if distance > s.Limits.MaxDistance {
		return nil, Stats{}, fmt.Errorf("%w: distance=%d > max=%d", ErrResourceLimit, distance, s.Limits.MaxDistance)
	}
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	if workers > s.Limits.MaxWorkers {
		return nil, Stats{}, fmt.Errorf("%w: workers=%d > max=%d", ErrResourceLimit, workers, s.Limits.MaxWorkers)
	}

	results := make([]TrialResult, shots)
	var mu sync.Mutex
	var firstErr error
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			for shot := workerID; shot < shots; shot += workers {
				rng := newChildRNG(seed, shot)
				e, err := code.SampleError(p, rng)
				if err != nil {
					mu.Lock()
					if firstErr == nil {
						firstErr = err
					}
					mu.Unlock()
					return
				}
				o := code.Decode(e)
				results[shot] = TrialResult{
					Error:   e,
					Outcome: o,
					Weight:  e.Weight(code.N),
				}
			}
		}(w)
	}
	wg.Wait()
	if firstErr != nil {
		return nil, Stats{}, firstErr
	}
	return results, Summarize(code, results), nil
}

// Stats 是批量试验的统计量。
type Stats struct {
	Shots            int
	LogicalErrors    int // 纠正后仍发生逻辑 X/Z 翻转的次数
	Uncorrectable    int // 被显式标记为超出保证纠正能力的次数
	Ambiguous        int
	LogicalErrorRate float64
	FlaggedRate      float64 // 逻辑翻转或显式标记不可靠的总比例
}

// Summarize 汇总试验结果。
func Summarize(code Code, results []TrialResult) Stats {
	st := Stats{Shots: len(results)}
	for _, r := range results {
		if r.Outcome.LogicalX || r.Outcome.LogicalZ {
			st.LogicalErrors++
		}
		if !r.Outcome.Correctable {
			st.Uncorrectable++
		}
		if r.Outcome.Ambiguous {
			st.Ambiguous++
		}
	}
	if st.Shots > 0 {
		st.LogicalErrorRate = float64(st.LogicalErrors) / float64(st.Shots)
		st.FlaggedRate = float64(st.LogicalErrors+st.Uncorrectable) / float64(st.Shots)
		if st.FlaggedRate > 1 {
			st.FlaggedRate = 1
		}
	}
	return st
}
