package circuit

import (
	"errors"
	"math"
	"math/rand/v2"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// TestInvalidCircuits 覆盖全部非法输入分类。
func TestInvalidCircuits(t *testing.T) {
	c, err := New(3)
	if err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name string
		g    Gate
		want error
	}{
		{"unknown gate", Gate{Name: "BOGUS", Target: 0}, ErrUnknownGate},
		{"target negative", Gate{Name: "X", Target: -1}, ErrInvalidQubit},
		{"target above range", Gate{Name: "X", Target: 3}, ErrInvalidQubit},
		{"control out of range", Gate{Name: "X", Target: 0, Controls: []int{9}}, ErrInvalidQubit},
		{"control equals target", Gate{Name: "X", Target: 1, Controls: []int{1}}, ErrControlConflict},
		{"duplicate controls", Gate{Name: "X", Target: 0, Controls: []int{2, 2}}, ErrControlConflict},
		{"NaN rotation", Gate{Name: "RX", Target: 0, Param: math.NaN()}, ErrInvalidParam},
		{"+Inf rotation", Gate{Name: "RZ", Target: 0, Param: math.Inf(1)}, ErrInvalidParam},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			err := c.Add(tc.g)
			if !errors.Is(err, tc.want) {
				t.Fatalf("Add(%v) err=%v, want %v", tc.g, err, tc.want)
			}
			t.Logf("业务[非法线路] 用例[%s] 结果[拒绝: %v]", tc.name, err)
		})
	}
	if _, err := New(0); !errors.Is(err, ErrInvalidQubit) {
		t.Fatalf("New(0) err=%v", err)
	}
}

// TestGateMath 精确核对固定门与旋转门的作用。
func TestGateMath(t *testing.T) {
	// H|0> => 概率各 1/2。
	c, _ := New(1)
	if err := c.H(0); err != nil {
		t.Fatal(err)
	}
	s, _ := NewState(1)
	if err := s.Run(c); err != nil {
		t.Fatal(err)
	}
	p := s.Probabilities()
	if !approx(p[0], 0.5, 1e-12) || !approx(p[1], 0.5, 1e-12) {
		t.Fatalf("H probabilities = %v", p)
	}
	// RX(pi)|0> = -i|1>
	c2, _ := New(1)
	if err := c2.RX(0, math.Pi); err != nil {
		t.Fatal(err)
	}
	s2, _ := NewState(1)
	_ = s2.Run(c2)
	p1, _ := s2.QubitProbability(0)
	if !approx(p1, 1, 1e-12) {
		t.Fatalf("RX(pi) p1=%v", p1)
	}
	// 受控门：CNOT |10> = |11>
	c3, _ := New(2)
	_ = c3.X(0)
	if err := c3.X(1, 0); err != nil {
		t.Fatal(err)
	}
	s3, _ := NewState(2)
	_ = s3.Run(c3)
	a, _ := s3.Amplitude(0b11)
	if approx(absC(a), 1, 1e-12) == false {
		t.Fatalf("CNOT amplitude on |11> = %v", a)
	}
	t.Logf("业务[门数学] 用例[H/RX/CNOT] 结果[通过]")
}

func absC(z complex128) float64 { return math.Hypot(real(z), imag(z)) }

// TestExpectation 核对 Pauli 张量积期望值。
func TestExpectation(t *testing.T) {
	c, _ := New(1)
	_ = c.H(0)
	s, _ := NewState(1)
	_ = s.Run(c)
	z, err := s.Expectation([]Pauli{PauliZ})
	if err != nil || !approx(z, 0, 1e-12) {
		t.Fatalf("<Z> on |+> = %v, err=%v", z, err)
	}
	x, _ := s.Expectation([]Pauli{PauliX})
	if !approx(x, 1, 1e-12) {
		t.Fatalf("<X> on |+> = %v", x)
	}
	// |1> 上 <Y>=0, <Z>=-1
	s1, _ := NewState(1)
	_ = s1.Run(mustCircuit(t, 1, func(cc *Circuit) { _ = cc.X(0) }))
	z1, _ := s1.Expectation([]Pauli{PauliZ})
	if !approx(z1, -1, 1e-12) {
		t.Fatalf("<Z> on |1> = %v", z1)
	}
	if _, err := s1.Expectation([]Pauli{PauliZ, PauliZ}); !errors.Is(err, ErrInvalidObservable) {
		t.Fatalf("wrong observable err=%v", err)
	}
	t.Logf("业务[期望值] 用例[<X>,<Z>,非法可观测量] 结果[通过]")
}

func mustCircuit(t *testing.T, n int, build func(*Circuit)) *Circuit {
	t.Helper()
	c, err := New(n)
	if err != nil {
		t.Fatal(err)
	}
	build(c)
	return c
}

// TestMeasurementCollapse 验证塌缩语义与可复现性。
func TestMeasurementCollapse(t *testing.T) {
	// 对 |+> 固定种子测量，结果序列必须逐位可复现。
	seq := func(seed uint64) []int {
		out := make([]int, 0, 20)
		rng := rand.New(rand.NewPCG(seed, seed^0x9e3779b97f4a7c15))
		for i := 0; i < 20; i++ {
			c, _ := New(1)
			_ = c.H(0)
			s, _ := NewState(1)
			_ = s.Run(c)
			r, err := s.Measure(0, rng)
			if err != nil {
				t.Fatal(err)
			}
			out = append(out, r)
		}
		return out
	}
	a := seq(12345)
	b := seq(12345)
	c := seq(12346)
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("same seed not reproducible at %d: %v vs %v", i, a, b)
		}
	}
	sameAsOther := true
	for i := range a {
		if a[i] != c[i] {
			sameAsOther = false
		}
	}
	if sameAsOther {
		t.Fatalf("different seeds produced identical sequences: %v", a)
	}
	// 塌缩后态必须是计算基本征态，且重复测量保持不变。
	circ, _ := New(1)
	_ = circ.H(0)
	s, _ := NewState(1)
	_ = s.Run(circ)
	rng := rand.New(rand.NewPCG(1, 2))
	r, _ := s.Measure(0, rng)
	p1, _ := s.QubitProbability(0)
	if !((r == 1 && approx(p1, 1, 1e-12)) || (r == 0 && approx(p1, 0, 1e-12))) {
		t.Fatalf("state did not collapse, r=%d p1=%v", r, p1)
	}
	r2, _ := s.Measure(0, rng)
	if r2 != r {
		t.Fatalf("repeated measurement changed result %d -> %d", r, r2)
	}
	// 统计频率应接近 1/2。
	ones := 0
	for _, v := range a {
		ones += v
	}
	if fr := float64(ones) / float64(len(a)); math.Abs(fr-0.5) > 0.25 {
		t.Fatalf("|+> measurement frequency drifted: %v", fr)
	}
	t.Logf("业务[测量塌缩] 用例[同种子复现/异种子不同/塌缩/重复测量] 结果[通过] 序列=%v", a)
}

func TestMeasureRequiresRNG(t *testing.T) {
	s, _ := NewState(1)
	if _, err := s.Measure(0, nil); !errors.Is(err, ErrNoRNG) {
		t.Fatalf("nil rng err=%v", err)
	}
}
