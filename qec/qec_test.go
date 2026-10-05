package qec

import (
	"errors"
	"testing"

	"github.com/lechandonga/quantum-qec-engine-go/circuit"
)

func TestInvalidCode(t *testing.T) {
	for _, d := range []int{1, 2, 4, 9} {
		if _, err := NewShorCode(d); !errors.Is(err, ErrInvalidCode) {
			t.Fatalf("NewShorCode(%d) err=%v want ErrInvalidCode", d, err)
		}
	}
	if _, err := NewShorCode(3); err != nil {
		t.Fatal(err)
	}
	t.Logf("业务[非法码距] 用例[1/2/4/9 拒绝, 3 接受] 结果[通过]")
}

// recoverLogicalZero 跑完整 编码-注入错误-解码-恢复 流程，
// 返回逻辑 |0> 保真度与解码结论。
func recoverLogicalZero(t *testing.T, code Code, e PauliError) (float64, Outcome) {
	t.Helper()
	en := NewEncoder(code)
	enc, err := en.Circuit()
	if err != nil {
		t.Fatal(err)
	}
	st, err := circuit.NewState(code.N)
	if err != nil {
		t.Fatal(err)
	}
	if err := st.Run(enc); err != nil {
		t.Fatal(err)
	}
	amp := st.Amplitudes()
	if err := en.ApplyError(amp, e); err != nil {
		t.Fatal(err)
	}
	o := code.Decode(e)
	rec, err := en.RecoveryCircuit(o)
	if err != nil {
		t.Fatal(err)
	}
	if err := rec.Apply(amp); err != nil {
		t.Fatal(err)
	}
	p0, err := en.LogicalPurity(amp)
	if err != nil {
		t.Fatal(err)
	}
	return p0, o
}

// TestCorrectableErrors d=3 时所有 27 种单比特 X/Y/Z 错误必须：
// 解码器承诺可纠正、无逻辑翻转标记，且编码态恢复保真度为 1。
func TestCorrectableErrors(t *testing.T) {
	code, err := NewShorCode(3)
	if err != nil {
		t.Fatal(err)
	}
	count := 0
	for q := 0; q < code.N; q++ {
		for _, e := range []PauliError{
			{Xmask: 1 << uint(q)},
			{Zmask: 1 << uint(q)},
			{Xmask: 1 << uint(q), Zmask: 1 << uint(q)},
		} {
			p0, o := recoverLogicalZero(t, code, e)
			if !o.Correctable || o.LogicalX || o.LogicalZ || p0 < 1-1e-9 {
				t.Fatalf("single error %v not corrected: %+v p0=%v", e, o, p0)
			}
			count++
		}
	}
	if count != 27 {
		t.Fatalf("checked %d cases, want 27", count)
	}
	t.Logf("业务[可纠正错误] 用例[d=3 全部 27 种单比特 X/Y/Z] 结果[全部纠正, 保真度=1]")
}

// TestUncorrectableErrors 纯逻辑算符与高权重错误必须被显式标记。
func TestUncorrectableErrors(t *testing.T) {
	code, _ := NewShorCode(3)

	// 纯逻辑 X̄ = 全局 X；纯逻辑 Z̄ = 块首 Z 链。
	logicalX := PauliError{Xmask: 1<<uint(code.N) - 1}
	logicalZ := PauliError{}
	for b := 0; b < code.Nz; b++ {
		logicalZ.Zmask |= 1 << uint(b*code.Nx)
	}
	for name, e := range map[string]PauliError{
		"logical X":     logicalX,
		"logical Z":     logicalZ,
		"2X same block": {Xmask: 1<<0 | 1<<1},
		"2Z two blocks": {Zmask: 1<<0 | 1<<3},
	} {
		o := code.Decode(e)
		if o.Correctable {
			t.Fatalf("%s wrongly marked correctable", name)
		}
		t.Logf("业务[不可纠正判定] 用例[%s 权重=%d] 结果[%s]", name, e.Weight(code.N), o.Message)
	}

	// 逻辑算符必须显式报告对应逻辑翻转。
	ox := code.Decode(logicalX)
	if !ox.LogicalX && !ox.LogicalZ {
		t.Fatalf("logical X not flagged: %+v", ox)
	}
	oz := code.Decode(logicalZ)
	if !oz.LogicalX && !oz.LogicalZ {
		t.Fatalf("logical Z not flagged: %+v", oz)
	}

	// 纯稳定子错误（块内 Z_j Z_{j+1}）应可纠正/等价无害。
	stab := PauliError{Zmask: 1<<1 | 1<<2}
	if o := code.Decode(stab); !o.Correctable || o.LogicalX || o.LogicalZ {
		t.Fatalf("stabilizer error mishandled: %+v", o)
	}
}

// TestReproducibleDecode 相同输入解码结论必须完全一致，
// 且固定种子的抽样试验逐次一致、与并发度无关。
func TestReproducibleDecode(t *testing.T) {
	code, _ := NewShorCode(3)
	e := PauliError{Xmask: 1<<2 | 1<<5, Zmask: 1 << 1}
	a := code.Decode(e)
	b := code.Decode(e)
	if a.Message != b.Message || a.Correctable != b.Correctable ||
		a.LogicalX != b.LogicalX || a.LogicalZ != b.LogicalZ ||
		a.Correction != b.Correction {
		t.Fatalf("decode not deterministic: %+v vs %+v", a, b)
	}
	sim := NewSimulator(DefaultLimits)
	r1, s1, err := sim.Run(3, 0.03, 3000, 2026, 1)
	if err != nil {
		t.Fatal(err)
	}
	r2, s2, err := sim.Run(3, 0.03, 3000, 2026, 8)
	if err != nil {
		t.Fatal(err)
	}
	if s1 != s2 || len(r1) != len(r2) {
		t.Fatalf("parallel simulation not deterministic: %+v vs %+v", s1, s2)
	}
	for i := range r1 {
		if r1[i].Error != r2[i].Error || r1[i].Outcome.Message != r2[i].Outcome.Message {
			t.Fatalf("trial %d differs", i)
		}
	}
	t.Logf("业务[解码可复现] 用例[固定错误解码 + 3000 样本 1/8 并发] 结果[完全一致, 逻辑错误率=%.5f]", s1.LogicalErrorRate)
}

// TestDistanceTrend 固定低错误率下逻辑错误率必须随码距单调下降。
func TestDistanceTrend(t *testing.T) {
	sim := NewSimulator(DefaultLimits)
	const shots = 60000
	const p = 0.02
	rates := map[int]float64{}
	for _, d := range []int{3, 5, 7} {
		_, st, err := sim.Run(d, p, shots, 13579, 8)
		if err != nil {
			t.Fatal(err)
		}
		rates[d] = st.LogicalErrorRate
		t.Logf("业务[码距趋势] 用例[d=%d, p=%.3f, shots=%d] 结果[逻辑错误率=%.5f 标记率=%.5f]",
			d, p, shots, st.LogicalErrorRate, st.FlaggedRate)
	}
	if !(rates[3] > rates[5] && rates[5] >= rates[7]) {
		t.Fatalf("monotonic decrease violated: %v", rates)
	}
}

// TestQECResourceLimits 超上限任务必须开始前拒绝。
func TestQECResourceLimits(t *testing.T) {
	sim := NewSimulator(Limits{MaxDistance: 5, MaxShots: 100, MaxWorkers: 4})
	if _, _, err := sim.Run(7, 0.01, 10, 1, 1); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("distance limit err=%v", err)
	}
	if _, _, err := sim.Run(3, 0.01, 101, 1, 1); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("shots limit err=%v", err)
	}
	if _, _, err := sim.Run(3, 0.01, 10, 1, 8); !errors.Is(err, ErrResourceLimit) {
		t.Fatalf("workers limit err=%v", err)
	}
	if _, _, err := sim.Run(3, 1.2, 10, 1, 1); !errors.Is(err, ErrInvalidErrorRate) {
		t.Fatalf("error rate limit err=%v", err)
	}
	t.Logf("业务[QEC资源上限] 用例[码距/样本/并发/错误率超限] 结果[全部开始前拒绝]")
}

func TestZeroErrorRate(t *testing.T) {
	sim := NewSimulator(DefaultLimits)
	_, st, err := sim.Run(3, 0, 5000, 1, 4)
	if err != nil {
		t.Fatal(err)
	}
	if st.LogicalErrors != 0 || st.Uncorrectable != 0 {
		t.Fatalf("zero error rate produced failures: %+v", st)
	}
	t.Logf("业务[零错误率退化] 用例[d=3, p=0, 5000 样本] 结果[逻辑错误=0, 不可纠正=0]")
}
