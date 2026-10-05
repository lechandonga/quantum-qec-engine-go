package qec

// stabilizerRows 返回 Shor 码独立稳定子生成元的二进制行。
// 行宽 2N：位 0..N-1 为 X 掩码，位 N..2N-1 为 Z 掩码。
// 受 PauliError 的 int 位掩码限制，N 必须 <= 约 60（见 NewShorCode）。
func (c Code) stabilizerRows() []uint {
	n := c.N
	rows := make([]uint, 0, 2*n)
	// 块内 Z 型稳定子 Z_{b,j} Z_{b,j+1}（定位 X 错误）。
	for b := 0; b < c.Nz; b++ {
		for j := 0; j < c.Nx-1; j++ {
			q0 := b*c.Nx + j
			q1 := q0 + 1
			rows = append(rows, uint(1<<uint(n+q0)|1<<uint(n+q1)))
		}
	}
	// 块间 X 型稳定子：相邻两块全体 X 的乘积（Nz-1 个独立生成元，
	// 定位块间相位错误）。全局 X（逻辑 X̄）在此空间之外。
	for b := 0; b < c.Nz-1; b++ {
		var row uint
		for k := 0; k < c.Nx; k++ {
			row |= 1 << uint(b*c.Nx+k)
			row |= 1 << uint((b+1)*c.Nx+k)
		}
		rows = append(rows, row)
	}
	return rows
}

// gf2Reduce 对 rows 做高斯消元得到行阶梯，并用它把 vec 约化为
// 稳定子陪集代表。返回 (约化结果, 行阶梯)；约化结果为 0 表示 vec
// 完全属于稳定子空间。
// gf2Reduce 对 rows 做 GF(2) 高斯-若当消元，返回简化行阶梯及其 pivot 列。
func gf2Basis(rows []uint, width int) ([]uint, []int) {
	basis := append([]uint(nil), rows...)
	pivots := make([]int, 0, len(rows))
	rank := 0
	for col := 0; col < width; col++ {
		mask := uint(1) << uint(col)
		pivot := -1
		for r := rank; r < len(basis); r++ {
			if basis[r]&mask != 0 {
				pivot = r
				break
			}
		}
		if pivot < 0 {
			continue
		}
		basis[rank], basis[pivot] = basis[pivot], basis[rank]
		for r := 0; r < len(basis); r++ {
			if r != rank && basis[r]&mask != 0 {
				basis[r] ^= basis[rank]
			}
		}
		pivots = append(pivots, col)
		rank++
	}
	return basis, pivots
}

func gf2Reduce(rows []uint, vec uint, width int) (uint, []uint) {
	basis, pivots := gf2Basis(rows, width)
	v := vec
	// 用记录的 pivot 列消去 vec（RREF 保证每个 pivot 列只有一行置位）。
	for i, row := range basis {
		if row == 0 || i >= len(pivots) {
			continue
		}
		if v&(uint(1)<<uint(pivots[i])) != 0 {
			v ^= row
		}
	}
	return v, basis
}

// logicalClass 计算给定 Pauli 在逻辑算符上的分量（陪集代表）。
// 先把错误约化到稳定子陪集代表 v=(x,z)，再与逻辑算符代表做辛内积：
//
//	<e, L_Z> = x·z_{LZ} ⊕ z·x_{LZ}（LZ 只有 Z 分量 => x·z_{LZ}）；
//	<e, L_X> = x·z_{LX} ⊕ z·x_{LX}（LX 只有 X 分量 => z·x_{LX}）。
//
// 内积为 1 表示残余错误反交换于该逻辑算符，即发生对应逻辑翻转。
//
//	逻辑 X̄ 代表链：每个块块内位置 1 的 X（跨全部 d 块）；
//	逻辑 Z̄ 代表链：块 0 全部比特的 Z。
func (c Code) logicalClass(r PauliError) (lx, lz bool) {
	n := c.N
	vec := uint(r.Xmask) | uint(r.Zmask)<<uint(n)
	v, _ := gf2Reduce(c.stabilizerRows(), vec, 2*n)
	xRep := int(v & (uint(1)<<uint(n) - 1))
	zRep := int(v >> uint(n))

	// 逻辑 Z̄：跨越全部 d 个块、每块块首一个 Z 的权重 d 链。
	// 与每条块间 X 稳定子重叠 2 点（对易），与全局 X 反对易。
	lzChain := 0
	for b := 0; b < c.Nz; b++ {
		lzChain |= 1 << uint(b*c.Nx)
	}
	// 逻辑 X̄：全局 X（N 个物理比特的 X 乘积），与所有稳定子对易。
	lxChain := 0
	for q := 0; q < c.N; q++ {
		lxChain |= 1 << uint(q)
	}
	// <P, L_Z> = x·z_LZ；<P, L_X> = z·x_LX。
	// 稳定子与两个逻辑算符都对易，故无需先做陪集约化，
	// 直接对原始 Pauli 计算内积即可，结果在同一稳定子陪集内恒定。
	// 返回值 (antiLx, antiLz)：
	// antiLx=true 表示反交换逻辑 X 可观测量（即发生逻辑 Z 翻转）；
	// antiLz=true 表示反交换逻辑 Z 可观测量（即发生逻辑 X 翻转）。
	antiLx := onesCount(uint(zRep&lxChain))%2 == 1
	antiLz := onesCount(uint(xRep&lzChain))%2 == 1
	return antiLx, antiLz
}

func onesCount(x uint) int {
	c := 0
	for x != 0 {
		x &= x - 1
		c++
	}
	return c
}
