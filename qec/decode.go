package qec

import "fmt"

// PauliError 用两个位掩码描述多比特 Pauli 错误：
// X 掩码与 Z 掩码同时置位的比特表示 Y。
type PauliError struct {
	Xmask int
	Zmask int
}

// Weight 返回错误的 Pauli 权重。
func (e PauliError) Weight(n int) int {
	w := 0
	for q := 0; q < n; q++ {
		b := 1 << q
		if e.Xmask&b != 0 || e.Zmask&b != 0 {
			w++
		}
	}
	return w
}

// Syndrome 是稳定子综合征。
type Syndrome struct {
	// XSyndrome[b][j] 为第 b 块内第 j 条 Z 型稳定子（比较比特 j 与
	// j+1 上的 X 奇偶）的测量值。
	XSyndrome [][]int
	// ZSyndrome[j] 为第 j 条 X 型稳定子（比较第 j 块与 j+1 块全体
	// 比特上的 Z 奇偶）的测量值。
	ZSyndrome []int
}

// ExtractSyndrome 直接由 Pauli 错误确定性地计算综合征。
func (c Code) ExtractSyndrome(e PauliError) Syndrome {
	syn := Syndrome{
		XSyndrome: make([][]int, c.Nz),
		ZSyndrome: make([]int, c.Nz-1),
	}
	for b := 0; b < c.Nz; b++ {
		syn.XSyndrome[b] = make([]int, c.Nx-1)
		for j := 0; j < c.Nx-1; j++ {
			q0 := b*c.Nx + j
			q1 := q0 + 1
			syn.XSyndrome[b][j] = bitAt(e.Xmask, q0) ^ bitAt(e.Xmask, q1)
		}
	}
	for j := 0; j < c.Nz-1; j++ {
		v := 0
		for k := 0; k < c.Nx; k++ {
			q0 := j*c.Nx + k
			q1 := q0 + c.Nx
			v ^= bitAt(e.Zmask, q0) ^ bitAt(e.Zmask, q1)
		}
		syn.ZSyndrome[j] = v
	}
	return syn
}

func bitAt(mask, q int) int { return (mask >> q) & 1 }

// Outcome 是一次解码结论。
type Outcome struct {
	Code       Code
	Syndrome   Syndrome
	Correction PauliError
	// Correctable 为 false 表示错误超出码的保证纠正半径 t=(d-1)/2，
	// 或落在最小权重错误不唯一的平局边界。
	Correctable bool
	// Ambiguous 表示综合征对应的最小权重错误不唯一。
	Ambiguous bool
	// LogicalX / LogicalZ 报告所选恢复下残余错误的逻辑分量。
	LogicalX bool
	LogicalZ bool
	Message  string
}

// repetitionPatternRaw 以 pattern[0]=0 锚定，逐位递推恢复重复码模式。
// 该模式与其整体补对应同一综合征，但可能相差一个逻辑算符。
func repetitionPatternRaw(syn []int) []int {
	L := len(syn) + 1
	p := make([]int, L)
	for j := 0; j < len(syn); j++ {
		p[j+1] = p[j] ^ syn[j]
	}
	return p
}

func complementPattern(p []int) []int {
	out := make([]int, len(p))
	for i, v := range p {
		out[i] = v ^ 1
	}
	return out
}

func patternWeight(p []int) int {
	w := 0
	for _, v := range p {
		w += v
	}
	return w
}

func repOnes(bits []int) int {
	n := 0
	for _, b := range bits {
		n += b
	}
	return n
}

// Decode 完成广义 Shor 码的多数判决解码。
//
// 判定规则（t=(d-1)/2）：
//   - X 分量：逐块按 Nx 长度重复码恢复；块内 X/Y 错误数 > t 即越界。
//   - Z 分量：先取每块 Z 奇偶作为块符号，再按 Nz 长度重复码恢复；
//     异常块数 > t 即越界。
//   - 综合征只确定错误“形状”，模式整体取补可能恰好等价于逻辑算符。
//     因此对 Z 恢复的两个互补方向分别计算残余的逻辑分量：
//     只有存在一个方向使残余完全位于稳定子空间，才承诺纠正成功；
//     否则显式标记 LogicalX/LogicalZ。
//   - 越界或平局时 Correctable=false，绝不静默给出成功结论。
func (c Code) Decode(e PauliError) Outcome {
	syn := c.ExtractSyndrome(e)
	out := Outcome{Code: c, Syndrome: syn, Correctable: true}
	t := (c.D - 1) / 2

	// --- X 分量恢复（逐块锚定模式） ----------------------------
	corr := PauliError{}
	xExceed := false
	for b := 0; b < c.Nz; b++ {
		actual := make([]int, c.Nx)
		for j := 0; j < c.Nx; j++ {
			actual[j] = bitAt(e.Xmask, b*c.Nx+j)
		}
		if n := repOnes(actual); n > t {
			xExceed = true
		}
		raw := repetitionPatternRaw(syn.XSyndrome[b])
		// 块内补自由度相差整块 X（对 X 重复码是逻辑算符），
		// 物理恢复选最小权重方向。
		pattern := raw
		if patternWeight(raw) > patternWeight(complementPattern(raw)) {
			pattern = complementPattern(raw)
		}
		for j, v := range pattern {
			if v == 1 {
				corr.Xmask |= 1 << uint(b*c.Nx+j)
			}
		}
	}

	// --- Z 分量恢复：块符号 + 两个互补方向评估 ------------------
	blockSigns := make([]int, c.Nz)
	for b := 0; b < c.Nz; b++ {
		s := 0
		for j := 0; j < c.Nx; j++ {
			s ^= bitAt(e.Zmask, b*c.Nx+j)
		}
		blockSigns[b] = s
	}
	zOnes := repOnes(blockSigns)
	zExceed := zOnes > t
	if zOnes == c.Nz-zOnes {
		out.Ambiguous = true
	}

	zA := repetitionPatternRaw(syn.ZSyndrome)
	zB := complementPattern(zA)
	corrA := corr
	corrA.Zmask = zBitsFromPattern(c, zA)
	corrB := corr
	corrB.Zmask = zBitsFromPattern(c, zB)

	lxA, lzA := c.logicalClass(PauliError{Xmask: e.Xmask ^ corrA.Xmask, Zmask: e.Zmask ^ corrA.Zmask})
	lxB, lzB := c.logicalClass(PauliError{Xmask: e.Xmask ^ corrB.Xmask, Zmask: e.Zmask ^ corrB.Zmask})

	// 两个互补恢复方向中，残余完全落入稳定子空间（逻辑分量为 0）
	// 的方向才能承诺成功。纯逻辑算符错误（平凡综合征）在两个方向上
	// 都会留下逻辑分量，从而被显式标记。
	// 特殊情形：若错误自身（仅扣掉 X 恢复）经 GF2 约化后非零、
	// 但综合征平凡（zA 全 0），说明错误就是裸逻辑算符，
	// 任何使其残余为 0 的恢复都等于逻辑翻转本身。
	trivialZSyn := patternWeight(zA) == 0
	selfLx, selfLz, pureLogical := false, false, false
	if trivialZSyn {
		reduced := cosetReduction(c, e, corr)
		if reduced != 0 {
			selfLx, selfLz = c.logicalClass(PauliError{Xmask: e.Xmask ^ corr.Xmask, Zmask: e.Zmask})
			pureLogical = selfLx || selfLz
		}
	}
	aValid := !lxA && !lzA && !pureLogical
	bValid := !lxB && !lzB && !pureLogical

	logical := true
	switch {
	case aValid:
		corr = corrA
		logical = false
	case bValid:
		corr = corrB
		logical = false
	default:
		// 两个方向都有逻辑分量：取较小权重恢复并报告该方向的翻转。
		if trivialZSyn {
			out.LogicalX, out.LogicalZ = selfLx, selfLz
		} else {
			out.LogicalX, out.LogicalZ = lxA, lzA
		}
		if patternWeight(zA) <= patternWeight(zB) {
			corr = corrA
			out.LogicalX, out.LogicalZ = lxA, lzA
		} else {
			corr = corrB
			out.LogicalX, out.LogicalZ = lxB, lzB
		}
	}

	out.Correction = corr
	if xExceed || zExceed {
		out.Correctable = false
	}
	if logical {
		out.Correctable = false
	}

	switch {
	case out.LogicalX || out.LogicalZ:
		out.Message = fmt.Sprintf("possible logical flip (X=%v, Z=%v)", out.LogicalX, out.LogicalZ)
	case out.Ambiguous:
		out.Message = "decoder boundary (non-unique minimum weight error)"
	case !out.Correctable:
		out.Message = "error weight exceeds guaranteed correction radius"
	default:
		out.Message = "corrected"
	}
	return out
}

// zBitsFromPattern 按块符号模式在每个异常块的块首放置一个 Z 恢复
// （块内单比特 Z 对块间相位码等价，块首为最小权重位置）。
func zBitsFromPattern(c Code, blockPattern []int) int {
	mask := 0
	for b, v := range blockPattern {
		if v == 1 {
			mask |= 1 << uint(b*c.Nx)
		}
	}
	return mask
}
