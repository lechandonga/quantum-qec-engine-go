package circuit

import (
	"math"
	"math/rand/v2"
)

// State 是 n 比特纯态的概率幅向量，长度 2^n。
type State struct {
	n   int
	amp []complex128
}

// NewState 创建全部处于 |0> 的初始态。
func NewState(n int) (*State, error) {
	if n <= 0 || n > 63 {
		return nil, ErrInvalidQubit
	}
	amp := make([]complex128, 1<<n)
	amp[0] = 1
	return &State{n: n, amp: amp}, nil
}

// NewStateFromAmps 用给定概率幅构造态，校验维度与归一化（容差 1e-9）。
func NewStateFromAmps(amp []complex128) (*State, error) {
	if len(amp) < 2 || len(amp)&(len(amp)-1) != 0 {
		return nil, ErrInvalidState
	}
	n := int(math.Log2(float64(len(amp))))
	norm := 0.0
	for _, a := range amp {
		norm += real(a)*real(a) + imag(a)*imag(a)
	}
	if math.Abs(norm-1) > 1e-9 {
		return nil, ErrInvalidState
	}
	cp := make([]complex128, len(amp))
	copy(cp, amp)
	return &State{n: n, amp: cp}, nil
}

// N 返回比特数。
func (s *State) N() int { return s.n }

// Amplitudes 返回概率幅的副本。
func (s *State) Amplitudes() []complex128 {
	out := make([]complex128, len(s.amp))
	copy(out, s.amp)
	return out
}

// Amplitude 返回指定计算基矢的概率幅；basis 越界返回错误。
func (s *State) Amplitude(basis int) (complex128, error) {
	if basis < 0 || basis >= len(s.amp) {
		return 0, ErrInvalidQubit
	}
	return s.amp[basis], nil
}

// Run 在线路当前态上演化线路。
func (s *State) Run(c *Circuit) error {
	if c.n != s.n {
		return ErrInvalidQubit
	}
	return c.Apply(s.amp)
}

// Probabilities 返回各计算基矢的测量概率（副本）。
func (s *State) Probabilities() []float64 {
	p := make([]float64, len(s.amp))
	for i, a := range s.amp {
		p[i] = real(a)*real(a) + imag(a)*imag(a)
	}
	return p
}

// QubitProbability 返回指定比特测得 1 的概率。
func (s *State) QubitProbability(q int) (float64, error) {
	if q < 0 || q >= s.n {
		return 0, ErrInvalidQubit
	}
	bit := 1 << q
	sum := 0.0
	for i, a := range s.amp {
		if i&bit != 0 {
			sum += real(a)*real(a) + imag(a)*imag(a)
		}
	}
	return sum, nil
}

// Measurement 是一次测量的结果。
type Measurement struct {
	Qubit  int
	Result int
}

// Measure 对指定比特做计算基测量并塌缩：返回 0/1，且将与结果矛盾的
// 分量清零、对剩余分量按测量后条件概率重新归一化。
// 随机源必须由调用方提供，保证同种子完全可复现。
func (s *State) Measure(q int, rng *rand.Rand) (int, error) {
	if q < 0 || q >= s.n {
		return 0, ErrInvalidQubit
	}
	if rng == nil {
		return 0, ErrNoRNG
	}
	p1, err := s.QubitProbability(q)
	if err != nil {
		return 0, err
	}
	res := 0
	if rng.Float64() < p1 {
		res = 1
	}
	s.collapse(q, res)
	return res, nil
}

// MeasureAll 按比特 0..n-1 依次测量，全部塌缩，返回每个比特的结果。
func (s *State) MeasureAll(rng *rand.Rand) ([]int, error) {
	if rng == nil {
		return nil, ErrNoRNG
	}
	out := make([]int, s.n)
	for q := 0; q < s.n; q++ {
		r, err := s.Measure(q, rng)
		if err != nil {
			return nil, err
		}
		out[q] = r
	}
	return out, nil
}

func (s *State) collapse(q, result int) {
	bit := 1 << q
	keep := 0.0
	for i, a := range s.amp {
		is1 := i&bit != 0
		if (result == 1) != is1 {
			s.amp[i] = 0
		} else {
			keep += real(a)*real(a) + imag(a)*imag(a)
		}
	}
	if keep == 0 {
		return
	}
	scale := complex(1/math.Sqrt(keep), 0)
	for i := range s.amp {
		s.amp[i] *= scale
	}
}

// Pauli 表示单比特 Pauli 算符。
type Pauli int

const (
	PauliI Pauli = iota
	PauliX
	PauliY
	PauliZ
)

// Expectation 计算可观测量 <psi| O |psi>，O 为各比特 Pauli 算符的张量积
// （term[q] 为比特 q 上的算符，PauliI 表示单位算符）。
func (s *State) Expectation(term []Pauli) (float64, error) {
	if len(term) != s.n {
		return 0, ErrInvalidObservable
	}
	for _, p := range term {
		if p < PauliI || p > PauliZ {
			return 0, ErrInvalidObservable
		}
	}
	var sum complex128
	for idx, a := range s.amp {
		j := idx
		ph := complex(1.0, 0.0)
		for q, p := range term {
			if p == PauliI {
				continue
			}
			bit := 1 << q
			switch p {
			case PauliX:
				j ^= bit
			case PauliY:
				j ^= bit
				// Y|0> = i|1>，Y|1> = -i|0>
				if idx&bit == 0 {
					ph *= complex(0, 1)
				} else {
					ph *= complex(0, -1)
				}
			case PauliZ:
				if idx&bit != 0 {
					ph = -ph
				}
			}
		}
		b := s.amp[j]
		sum += ph * complex(real(b)*real(a)+imag(b)*imag(a), real(b)*imag(a)-imag(b)*real(a))
	}
	return real(sum), nil
}
