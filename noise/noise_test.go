package noise

import (
	"errors"
	"math"
	"testing"

	"github.com/lechandonga/quantum-qec-engine-go/circuit"
)

func hadamardCircuit(t *testing.T) *circuit.Circuit {
	t.Helper()
	c, err := circuit.New(1)
	if err != nil {
		t.Fatal(err)
	}
	if err := c.H(0); err != nil {
		t.Fatal(err)
	}
	return c
}

// TestInvalidModels 覆盖非法通道名、强度与比特下标。
func TestInvalidModels(t *testing.T) {
	e := NewEngine(DefaultLimits)
	c := hadamardCircuit(t)
	cases := []struct {
		name  string
		model Model
		want  error
	}{
		{"unknown channel", Model{Channel: "teleport", Strength: 0.1}, ErrInvalidModel},
		{"bit flip negative", Model{Channel: BitFlip, Strength: -0.01}, ErrInvalidStrength},
		{"bit flip > 1", Model{Channel: BitFlip, Strength: 1.1}, ErrInvalidStrength},
		{"phase > 1", Model{Channel: PhaseFlip, Strength: 1.0001}, ErrInvalidStrength},
		{"depol > 4/3", Model{Channel: Depolarizing, Strength: 1.34}, ErrInvalidStrength},
		{"bad qubit", Model{Channel: BitFlip, Strength: 0.1, Qubits: []int{7}}, circuit.ErrInvalidQubit},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			_, _, err := e.ExactProbabilities(c, tc.model)
			if !errors.Is(err, tc.want) {
				t.Fatalf("err=%v want %v", err, tc.want)
			}
			t.Logf("业务[非法噪声] 用例[%s] 结果[拒绝: %v]", tc.name, err)
		})
	}
}

// TestZeroStrengthIdentity 强度为零时精确结果必须与无噪线路逐位一致，
// 抽样结果也必须完全落在无噪输出上。
func TestZeroStrengthIdentity(t *testing.T) {
	e := NewEngine(DefaultLimits)
	for _, ch := range []Channel{BitFlip, PhaseFlip, Depolarizing} {
		c, _ := circuit.New(2)
		_ = c.H(0)
		_ = c.X(1, 0)
		p, tr, err := e.ExactProbabilities(c, Model{Channel: ch, Strength: 0})
		if err != nil {
			t.Fatal(err)
		}
		// 无噪: H|0> 后 CNOT => Bell 态，|00>,|11> 各 1/2。
		for basis, want := range map[int]float64{0: 0.5, 3: 0.5} {
			if !approxT(p[basis], want, 1e-12) {
				t.Fatalf("%s p=%v", ch, p)
			}
		}
		if !approxT(tr, 1, 1e-12) {
			t.Fatalf("%s trace=%v", ch, tr)
		}
		// 抽样应只出现 0 和 3。
		batch, err := e.Batch(c, Model{Channel: ch, Strength: 0}, 500, 99, 4)
		if err != nil {
			t.Fatal(err)
		}
		for basis := range batch.Counts {
			if basis != 0 && basis != 3 {
				t.Fatalf("%s produced unexpected basis %d", ch, basis)
			}
		}
		t.Logf("业务[零强度退化] 用例[%s] 结果[精确分布与抽样均等于无噪]", ch)
	}
}

func approxT(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// TestExactSmallSolutions 与可手算的小规模解析结果逐项对比。
func TestExactSmallSolutions(t *testing.T) {
	e := NewEngine(DefaultLimits)

	// 比特翻转 p 作用在 |0>（无门）上：P(1)=p。
	idle, _ := circuit.New(1)
	p, tr, err := e.ExactProbabilities(idle, Model{Channel: BitFlip, Strength: 0.2})
	if err != nil {
		t.Fatal(err)
	}
	if !approxT(p[1], 0.2, 1e-12) {
		t.Fatalf("bitflip P(1)=%v", p[1])
	}
	if !approxT(tr, 1, 1e-12) {
		t.Fatalf("trace drift %v", tr)
	}

	// 退极化 p 作用在 |0> 上：X/Y 共概率 2p/3 翻到 |1>。
	pd, _, _ := e.ExactProbabilities(idle, Model{Channel: Depolarizing, Strength: 0.3})
	if !approxT(pd[1], 0.2, 1e-12) {
		t.Fatalf("depol P(1)=%v want 0.2", pd[1])
	}

	// 相位翻转不改变计算基概率。
	pp, _, _ := e.ExactProbabilities(idle, Model{Channel: PhaseFlip, Strength: 0.7})
	if !approxT(pp[0], 1, 1e-12) {
		t.Fatalf("phase flip changed probabilities %v", pp)
	}

	// H - 相位翻转 - H：Z 噪声在 X 基制备下表现为翻转，P(1)=p。
	first, _ := circuit.New(1)
	_ = first.H(0)
	second, _ := circuit.New(1)
	_ = second.H(0)
	ph, tr2, _ := e.ExactWithMidCircuitNoise(first, second, Model{Channel: PhaseFlip, Strength: 0.3})
	if !approxT(ph[1], 0.3, 1e-12) {
		t.Fatalf("H-phase-H P(1)=%v want 0.3", ph[1])
	}
	if !approxT(tr2, 1, 1e-12) {
		t.Fatalf("trace drift %v", tr2)
	}
	t.Logf("业务[精确解析解] 用例[比特翻转/退极化/相位翻转+迹] 结果[通过]")
}

// TestStatisticalConvergence 大样本统计量必须收敛到精确分布。
func TestStatisticalConvergence(t *testing.T) {
	e := NewEngine(DefaultLimits)
	c, _ := circuit.New(1)
	_ = c.H(0)
	model := Model{Channel: BitFlip, Strength: 0.25}
	exact, tr, err := e.ExactProbabilities(c, model)
	if err != nil {
		t.Fatal(err)
	}
	if !approxT(tr, 1, 1e-12) {
		t.Fatalf("trace drift")
	}
	batch, err := e.Batch(c, model, 40000, 202602, 8)
	if err != nil {
		t.Fatal(err)
	}
	d := TotalVariation(exact, batch.Empirical)
	if d > 0.01 {
		t.Fatalf("empirical distribution did not converge, TVD=%v exact=%v emp=%v", d, exact, batch.Empirical)
	}
	t.Logf("业务[统计收敛] 用例[40000 样本 vs 精确分布] 结果[TVD=%.4f, 通过]", d)
}

// TestParallelDeterminism 相同种子下结果必须与并发度无关、逐样本一致。
func TestParallelDeterminism(t *testing.T) {
	e := NewEngine(DefaultLimits)
	c, _ := circuit.New(2)
	_ = c.H(0)
	_ = c.X(1, 0)
	model := Model{Channel: Depolarizing, Strength: 0.08}
	var ref map[int]int
	for _, workers := range []int{1, 2, 3, 8, 16} {
		b, err := e.Batch(c, model, 2000, 424242, workers)
		if err != nil {
			t.Fatal(err)
		}
		if ref == nil {
			ref = b.Counts
			continue
		}
		if len(ref) != len(b.Counts) {
			t.Fatalf("worker=%d counts differ", workers)
		}
		for k, v := range ref {
			if b.Counts[k] != v {
				t.Fatalf("worker=%d basis %d: %d vs %d", workers, k, v, b.Counts[k])
			}
		}
	}
	// 逐样本也必须一致（不同并发分片）。
	s1, err := e.Sample(c, model, 777, 123)
	if err != nil {
		t.Fatal(err)
	}
	s2, err := e.Sample(c, model, 777, 123)
	if err != nil {
		t.Fatal(err)
	}
	if s1.Basis != s2.Basis {
		t.Fatalf("per-shot not reproducible: %d vs %d", s1.Basis, s2.Basis)
	}
	t.Logf("业务[并行可复现] 用例[并发度 1/2/3/8/16 计数完全一致, 逐样本一致] 结果[通过]")
}

// TestResourceLimits 超限输入必须在计算开始前被拒绝。
func TestResourceLimits(t *testing.T) {
	e := NewEngine(Limits{MaxQubits: 4, MaxShots: 100, MaxWorkers: 4})
	big, err := circuit.New(5)
	if err != nil {
		t.Fatal(err)
	}
	if _, _, err := e.ExactProbabilities(big, Model{Channel: BitFlip}); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("qubits limit err=%v", err)
	}
	small := hadamardCircuit(t)
	if _, err := e.Batch(small, Model{Channel: BitFlip}, 101, 1, 1); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("shots limit err=%v", err)
	}
	if _, err := e.Batch(small, Model{Channel: BitFlip}, 10, 1, 5); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("workers limit err=%v", err)
	}
	// 零样本合法。
	b, err := e.Batch(small, Model{Channel: BitFlip}, 0, 1, 1)
	if err != nil || b.Shots != 0 {
		t.Fatalf("zero shots err=%v", err)
	}
	t.Logf("业务[资源上限] 用例[比特数/样本数/并发度超限 + 零样本] 结果[拒绝/通过]")
}
