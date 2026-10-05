package noise

import (
	"errors"
	"fmt"
	"math"

	"github.com/lechandonga/quantum-qec-engine-go/circuit"
)

// 可区分的错误类型。
var (
	ErrInvalidStrength = errors.New("noise: channel strength must be within [0,1] (depolarizing within [0, 4/3])")
	ErrInvalidModel    = errors.New("noise: unknown channel name")
	ErrResourceLimit   = errors.New("noise: requested workload exceeds configured resource limits")
	ErrInvalidCircuit  = errors.New("noise: invalid circuit")
)

// Channel 是作用在单比特上的 Pauli 噪声通道名称。
type Channel string

const (
	BitFlip      Channel = "bit_flip"     // X 错误，概率 p
	PhaseFlip    Channel = "phase_flip"   // Z 错误，概率 p
	Depolarizing Channel = "depolarizing" // 等概率 X/Y/Z，单比特约定每类概率 p/3
)

// Model 描述一个均匀噪声模型：线路结束后对指定比特（为空表示全部比特）
// 独立施加同一通道。
type Model struct {
	Channel  Channel
	Strength float64
	Qubits   []int // 空切片表示全部比特
}

func (m Model) validate(n int) error {
	if m.Channel != BitFlip && m.Channel != PhaseFlip && m.Channel != Depolarizing {
		return fmt.Errorf("%w: %q", ErrInvalidModel, m.Channel)
	}
	limit := 1.0
	if m.Channel == Depolarizing {
		limit = 4.0 / 3.0
	}
	if math.IsNaN(m.Strength) || m.Strength < 0 || m.Strength > limit+1e-12 {
		return fmt.Errorf("%w: %s strength=%g", ErrInvalidStrength, m.Channel, m.Strength)
	}
	for _, q := range m.Qubits {
		if q < 0 || q >= n {
			return circuit.ErrInvalidQubit
		}
	}
	return nil
}

func (m Model) targets(n int) []int {
	if len(m.Qubits) > 0 {
		out := append([]int(nil), m.Qubits...)
		// 去重，避免同一比特施加两次噪声。
		seen := map[int]struct{}{}
		dedup := out[:0]
		for _, q := range out {
			if _, ok := seen[q]; ok {
				continue
			}
			seen[q] = struct{}{}
			dedup = append(dedup, q)
		}
		return dedup
	}
	out := make([]int, n)
	for i := range out {
		out[i] = i
	}
	return out
}

// sampleError 在给定比特上按通道抽样一个 Pauli 错误。
// 返回值为 xmask/zmask 分量（单比特，故为 0 或 1<<q 的掩码分量）。
func (m Model) sampleError(rng interface {
	Float64() float64
}, q int) (x, z bool) {
	r := rng.Float64()
	switch m.Channel {
	case BitFlip:
		return r < m.Strength, false
	case PhaseFlip:
		return false, r < m.Strength
	case Depolarizing:
		p := m.Strength / 3
		switch {
		case r < p:
			return true, false // X
		case r < 2*p:
			return true, true // Y
		case r < 3*p:
			return false, true // Z
		}
	}
	return false, false
}

// krausTerms 返回通道的非平凡 Kraus 项（Pauli + 概率），
// 供密度矩阵精确演化使用。
type krausTerm struct {
	prob float64
	x, z bool
}

func (m Model) krausTerms() []krausTerm {
	switch m.Channel {
	case BitFlip:
		return []krausTerm{{m.Strength, true, false}}
	case PhaseFlip:
		return []krausTerm{{m.Strength, false, true}}
	case Depolarizing:
		p := m.Strength / 3
		return []krausTerm{{p, true, false}, {p, true, true}, {p, false, true}}
	default:
		return nil
	}
}
