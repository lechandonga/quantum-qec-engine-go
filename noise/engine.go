package noise

import (
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/lechandonga/quantum-qec-engine-go/circuit"
)

// Limits 是可配置的资源上限。任何超出上限的任务必须在开始计算前拒绝，
// 不允许崩掉或静默截断。零值字段表示沿用 DefaultLimits。
type Limits struct {
	MaxQubits  int // 线路最大比特数
	MaxShots   int // 单次批量抽样最大样本数
	MaxWorkers int // 最大并行度
}

// DefaultLimits 为默认资源上限。
var DefaultLimits = Limits{
	MaxQubits:  20, // 精确密度矩阵为 4^n，20 比特时抽样路径仍便宜
	MaxShots:   1_000_000,
	MaxWorkers: 64,
}

// Engine 在线路之上叠加噪声模型，提供精确（密度矩阵）与抽样两种口径。
type Engine struct {
	Limits Limits
}

// NewEngine 创建引擎，零值字段回填默认上限。
func NewEngine(l Limits) *Engine {
	if l.MaxQubits <= 0 {
		l.MaxQubits = DefaultLimits.MaxQubits
	}
	if l.MaxShots <= 0 {
		l.MaxShots = DefaultLimits.MaxShots
	}
	if l.MaxWorkers <= 0 {
		l.MaxWorkers = DefaultLimits.MaxWorkers
	}
	return &Engine{Limits: l}
}

func (e *Engine) check(c *circuit.Circuit, shots int) error {
	if c == nil {
		return ErrInvalidCircuit
	}
	if c.N() > e.Limits.MaxQubits {
		return fmt.Errorf("%w: qubits=%d > max=%d", ErrResourceLimit, c.N(), e.Limits.MaxQubits)
	}
	if shots < 0 || shots > e.Limits.MaxShots {
		return fmt.Errorf("%w: shots=%d > max=%d", ErrResourceLimit, shots, e.Limits.MaxShots)
	}
	return nil
}

// ExactProbabilities 先执行无噪线路，再逐比特施加 Kraus 通道，
// 对密度矩阵做精确演化，返回计算基测量概率与迹（应恒为 1）。
func (e *Engine) ExactProbabilities(c *circuit.Circuit, model Model) ([]float64, float64, error) {
	if err := e.check(c, 0); err != nil {
		return nil, 0, err
	}
	n := c.N()
	if err := model.validate(n); err != nil {
		return nil, 0, err
	}
	st, err := circuit.NewState(n)
	if err != nil {
		return nil, 0, err
	}
	if err := st.Run(c); err != nil {
		return nil, 0, err
	}
	dim := 1 << n
	rho := densityFromPure(st.Amplitudes())

	// 先把门演化施加到密度矩阵上（st 已经是纯态演化结果，这里 rho 即 U|0>）。
	// 噪声通道按 Kraus 和 rho' = (1-p)rho + p P rho P 逐比特展开。
	for _, q := range model.targets(n) {
		rho = applyChannelDensity(rho, dim, n, q, model)
	}
	tr := traceDensity(rho, dim)
	if math.Abs(tr-1) > 1e-9 {
		return nil, tr, fmt.Errorf("noise: trace drift detected, trace=%.12f", tr)
	}
	return densityToProbabilities(rho, dim), tr, nil
}

// ExactWithMidCircuitNoise 执行 prefix，施加噪声，再执行 suffix，
// 返回计算基测量概率与迹。用于验证噪声基变换效应（例如 H-相位-H）。
func (e *Engine) ExactWithMidCircuitNoise(prefix, suffix *circuit.Circuit, model Model) ([]float64, float64, error) {
	if prefix == nil || suffix == nil || prefix.N() != suffix.N() {
		return nil, 0, ErrInvalidCircuit
	}
	if err := e.check(prefix, 0); err != nil {
		return nil, 0, err
	}
	n := prefix.N()
	if err := model.validate(n); err != nil {
		return nil, 0, err
	}
	st, err := circuit.NewState(n)
	if err != nil {
		return nil, 0, err
	}
	if err := st.Run(prefix); err != nil {
		return nil, 0, err
	}
	dim := 1 << n
	rho := densityFromPure(st.Amplitudes())
	for _, q := range model.targets(n) {
		rho = applyChannelDensity(rho, dim, n, q, model)
	}
	// suffix 门作用在密度矩阵上。
	gates := suffix.Gates()
	for _, g := range gates {
		rho = applyCircuitGateDensity(rho, dim, g)
	}
	tr := traceDensity(rho, dim)
	if math.Abs(tr-1) > 1e-9 {
		return nil, tr, fmt.Errorf("noise: trace drift detected, trace=%.12f", tr)
	}
	return densityToProbabilities(rho, dim), tr, nil
}

func applyChannelDensity(rho []complex128, dim, n, q int, model Model) []complex128 {
	terms := model.krausTerms()
	survive := 1.0
	for _, t := range terms {
		survive -= t.prob
	}
	out := make([]complex128, len(rho))
	for i := range out {
		out[i] = complex(survive, 0) * rho[i]
	}
	for _, t := range terms {
		xm, zm := 0, 0
		if t.x {
			xm = 1 << q
		}
		if t.z {
			zm = 1 << q
		}
		// P rho P：对行、列基矢同时作用 P。利用 rho_{ij} -> s_i s_j rho_{i^x,j^x}。
		next := make([]complex128, len(rho))
		for i := 0; i < dim; i++ {
			ip := i ^ xm
			si := zSign(i, zm)
			for j := 0; j < dim; j++ {
				jp := j ^ xm
				sj := zSign(j, zm)
				next[ip*dim+jp] += complex(float64(si*sj), 0) * rho[i*dim+j]
			}
		}
		for i := range out {
			out[i] += complex(t.prob, 0) * next[i]
		}
	}
	return out
}

func zSign(idx, zmask int) int {
	if idx&zmask != 0 {
		return -1
	}
	return 1
}

// ShotResult 是单次抽样结果。Counts 为计算基整数出现次数。
type ShotResult struct {
	Bits  []int // 每个比特的测量值
	Basis int   // 计算基整数
}

// Sample 运行一次「线路 + 噪声 + 测量」。噪声在门演化后、测量前逐比特抽样。
func (e *Engine) Sample(c *circuit.Circuit, model Model, seed uint64, shot int) (ShotResult, error) {
	if err := e.check(c, 0); err != nil {
		return ShotResult{}, err
	}
	n := c.N()
	if err := model.validate(n); err != nil {
		return ShotResult{}, err
	}
	st, err := circuit.NewState(n)
	if err != nil {
		return ShotResult{}, err
	}
	if err := st.Run(c); err != nil {
		return ShotResult{}, err
	}
	rng := newChildRNG(seed, shot)
	amp := st.Amplitudes()
	for _, q := range model.targets(n) {
		x, z := model.sampleError(rng, q)
		xm, zm := 0, 0
		if x {
			xm = 1 << q
		}
		if z {
			zm = 1 << q
		}
		applyPauli(amp, n, xm, zm)
	}
	noisy, err := circuit.NewStateFromAmps(amp)
	if err != nil {
		return ShotResult{}, err
	}
	bits, err := noisy.MeasureAll(rng)
	if err != nil {
		return ShotResult{}, err
	}
	basis := 0
	for q, b := range bits {
		basis |= b << q
	}
	return ShotResult{Bits: bits, Basis: basis}, nil
}

// BatchResult 汇总批量抽样。
type BatchResult struct {
	Shots     int
	Counts    map[int]int
	Empirical []float64 // 每个计算基的经验概率
}

// Batch 并行执行 shots 次抽样。结果只取决于 (seed, shot 编号)，
// 与 workers 数量和调度无关；workers<=0 时取 runtime.NumCPU。
func (e *Engine) Batch(c *circuit.Circuit, model Model, shots int, seed uint64, workers int) (BatchResult, error) {
	if err := e.check(c, shots); err != nil {
		return BatchResult{}, err
	}
	n := c.N()
	if err := model.validate(n); err != nil {
		return BatchResult{}, err
	}
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	if workers > e.Limits.MaxWorkers {
		return BatchResult{}, fmt.Errorf("%w: workers=%d > max=%d", ErrResourceLimit, workers, e.Limits.MaxWorkers)
	}
	if shots == 0 {
		return BatchResult{Shots: 0, Counts: map[int]int{}, Empirical: make([]float64, 1<<n)}, nil
	}

	counts := make([]int, 1<<n)
	var mu sync.Mutex
	var wg sync.WaitGroup
	var firstErr error
	for w := 0; w < workers; w++ {
		start := w
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			local := make([]int, 1<<n)
			for shot := workerID; shot < shots; shot += workers {
				r, err := e.Sample(c, model, seed, shot)
				if err != nil {
					mu.Lock()
					if firstErr == nil {
						firstErr = err
					}
					mu.Unlock()
					return
				}
				local[r.Basis]++
			}
			mu.Lock()
			for k, v := range local {
				counts[k] += v
			}
			mu.Unlock()
		}(start)
	}
	wg.Wait()
	if firstErr != nil {
		return BatchResult{}, firstErr
	}

	out := BatchResult{Shots: shots, Counts: map[int]int{}, Empirical: make([]float64, 1<<n)}
	for basis, v := range counts {
		if v > 0 {
			out.Counts[basis] = v
		}
		out.Empirical[basis] = float64(v) / float64(shots)
	}
	return out, nil
}

// TotalVariation 返回两个概率分布的总变差距离 1/2 Σ|p-q|。
func TotalVariation(p, q []float64) float64 {
	d := 0.0
	for i := range p {
		d += math.Abs(p[i] - q[i])
	}
	return d / 2
}
