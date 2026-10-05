// Package qec 实现基于广义 Shor 码 [[d^2,1,d]] 的稳定子量子纠错。
//
// 编码结构：n_z 个相位重复块，每块含 n_x 个比特翻转重复比特。
// 本包约定取奇数 n_x = n_z = d（码距），物理比特数 N = d^2。
// 物理错误采用独立同分布的 X/Y/Z 对称 Pauli 噪声（每比特每种错误
// 概率 p/3，与 noise.Depolarizing 一致）。
package qec

import (
	"errors"
	"fmt"
)

var (
	// ErrInvalidCode 表示码参数非法。
	ErrInvalidCode = errors.New("qec: invalid code parameters")
	// ErrInvalidErrorRate 表示物理错误率越界。
	ErrInvalidErrorRate = errors.New("qec: error rate must be within [0,1]")
	// ErrResourceLimit 表示任务超出资源上限。
	ErrResourceLimit = errors.New("qec: workload exceeds configured resource limits")
)

// Code 是一个广义 Shor 码实例。
type Code struct {
	D  int // 码距（奇数）
	N  int // 物理比特数 d*d
	Nx int // 块内比特翻转重复长度
	Nz int // 相位重复块数
}

// NewShorCode 构造码距为 d 的广义 Shor 码。d 必须为不小于 3 的奇数。
func NewShorCode(d int) (Code, error) {
	if d < 3 || d%2 == 0 {
		return Code{}, fmt.Errorf("%w: distance must be an odd integer >= 3, got %d", ErrInvalidCode, d)
	}
	// PauliError 以 Go int 位掩码表示错误，N=d*d 必须能放入掩码；
	// d<=7（N=49）时既保留码距扩展空间，又远低于掩码上限。
	if d > 7 {
		return Code{}, fmt.Errorf("%w: distance %d exceeds supported maximum 7 (bit-mask limit)", ErrInvalidCode, d)
	}
	return Code{D: d, N: d * d, Nx: d, Nz: d}, nil
}

// blockOf 返回物理比特 q 所属的相位块编号。
func (c Code) blockOf(q int) int { return q / c.Nx }

// withinBlock 返回块内位置。
func (c Code) withinBlock(q int) int { return q % c.Nx }
