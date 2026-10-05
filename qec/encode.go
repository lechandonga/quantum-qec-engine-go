package qec

import (
	"fmt"

	"github.com/lechandonga/quantum-qec-engine-go/circuit"
)

// Encoder 构建广义 Shor 码的编码电路与恢复工具，用于编码态级验证。
//
// 约定物理比特编号：块 b 内位置 j 的比特下标为 b*Nx + j。
// |0>_L 与 |1>_L 的编码：
//  1. 数据比特 0 作为逻辑输入；先做块间相位重复（块首间 CNOT 扇出）；
//  2. 每个块首施加 H；
//  3. 块首向块内其余比特 CNOT 扇出（比特翻转重复）。
//
// 逻辑态 |ψ> = a|0> + b|1> 放在比特 0 即可直接编码。
type Encoder struct {
	Code Code
}

// NewEncoder 创建编码器。
func NewEncoder(code Code) *Encoder { return &Encoder{Code: code} }

// Circuit 返回编码线路（共 N 个比特，逻辑输入在比特 0）。
func (en *Encoder) Circuit() (*circuit.Circuit, error) {
	c, err := circuit.New(en.Code.N)
	if err != nil {
		return nil, err
	}
	// 块间 CNOT：块首 0 -> b*Nx。
	for b := 1; b < en.Code.Nz; b++ {
		if err := c.X(b*en.Code.Nx, 0); err != nil {
			return nil, err
		}
	}
	// 块首 H。
	for b := 0; b < en.Code.Nz; b++ {
		if err := c.H(b * en.Code.Nx); err != nil {
			return nil, err
		}
	}
	// 块内 CNOT 扇出。
	for b := 0; b < en.Code.Nz; b++ {
		head := b * en.Code.Nx
		for j := 1; j < en.Code.Nx; j++ {
			if err := c.X(head+j, head); err != nil {
				return nil, err
			}
		}
	}
	return c, nil
}

// LogicalZeroState 返回编码 |0>_L 的理想概率幅（全部校验位约束下的
// 均匀叠加态），供保真度比对。
func (en *Encoder) LogicalZeroState() ([]complex128, error) {
	c, err := en.Circuit()
	if err != nil {
		return nil, err
	}
	st, err := circuit.NewState(en.Code.N)
	if err != nil {
		return nil, err
	}
	if err := st.Run(c); err != nil {
		return nil, err
	}
	return st.Amplitudes(), nil
}

// ApplyError 将抽样/给定的 Pauli 错误作用到概率幅上。
func (en *Encoder) ApplyError(amp []complex128, e PauliError) error {
	if len(amp) != 1<<en.Code.N {
		return fmt.Errorf("qec: amplitude dimension mismatch")
	}
	applyPauliAmp(amp, en.Code.N, e.Xmask, e.Zmask)
	return nil
}

// RecoveryCircuit 根据解码结论构建恢复线路：施加纠正 Pauli，
// 再按编码的逆序解码，比特 0 输出恢复后的逻辑态。
func (en *Encoder) RecoveryCircuit(o Outcome) (*circuit.Circuit, error) {
	c, err := circuit.New(en.Code.N)
	if err != nil {
		return nil, err
	}
	// 先施加纠正（X 与 Z 分量；同位即 Y）。
	for q := 0; q < en.Code.N; q++ {
		b := 1 << q
		if o.Correction.Xmask&b != 0 {
			if err := c.X(q); err != nil {
				return nil, err
			}
		}
		if o.Correction.Zmask&b != 0 {
			if err := c.Z(q); err != nil {
				return nil, err
			}
		}
	}
	// 逆编码：块内 CNOT（CNOT 自逆）、块首 H、块间 CNOT。
	for b := en.Code.Nz - 1; b >= 0; b-- {
		head := b * en.Code.Nx
		for j := en.Code.Nx - 1; j >= 1; j-- {
			if err := c.X(head+j, head); err != nil {
				return nil, err
			}
		}
	}
	for b := en.Code.Nz - 1; b >= 0; b-- {
		if err := c.H(b * en.Code.Nx); err != nil {
			return nil, err
		}
	}
	for b := en.Code.Nz - 1; b >= 1; b-- {
		if err := c.X(b*en.Code.Nx, 0); err != nil {
			return nil, err
		}
	}
	return c, nil
}

// LogicalPurity 返回恢复后逻辑比特（比特 0）处于 |0> 的概率，
// 对 |0>_L 编码-恢复实验即成功保真度。
func (en *Encoder) LogicalPurity(amp []complex128) (float64, error) {
	st, err := circuit.NewStateFromAmps(amp)
	if err != nil {
		return 0, err
	}
	p1, err := st.QubitProbability(0)
	if err != nil {
		return 0, err
	}
	return 1 - p1, nil
}

// applyPauliAmp 对概率幅直接施加 Pauli（全局相位无关）。
func applyPauliAmp(amp []complex128, n, xmask, zmask int) {
	if xmask == 0 {
		for idx := range amp {
			if idx&zmask != 0 {
				amp[idx] = -amp[idx]
			}
		}
		return
	}
	out := make([]complex128, len(amp))
	for idx, a := range amp {
		j := idx ^ xmask
		s := 1
		for q := 0; q < n; q++ {
			b := 1 << q
			if xmask&b != 0 && zmask&b != 0 && idx&b == 0 {
				s = -s
			}
		}
		out[j] += complex(float64(s), 0) * a
	}
	copy(amp, out)
}
