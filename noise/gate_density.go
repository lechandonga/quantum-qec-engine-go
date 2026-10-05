package noise

import "github.com/lechandonga/quantum-qec-engine-go/circuit"

// applyCircuitGateDensity 将一条（含受控的）门作用到密度矩阵 rho -> U rho U†。
// 对每一对 (target=0/1) 且控制位满足的基矢 lo/hi，rho 中行、列指标位于
// {lo,hi} 的 2x2 子块按 U 做同构变换；其余子块仅做整体重排，
// 这里采用遍历全部 2x2 子块的方式逐门更新。
func applyCircuitGateDensity(rho []complex128, dim int, g circuit.Gate) []complex128 {
	m, err := circuit.GateMatrix(g)
	if err != nil {
		panic(err)
	}
	tbit := 1 << g.Target
	cmask := 0
	for _, q := range g.Controls {
		cmask |= 1 << q
	}

	// 选出“控制条件满足且 target 固定为 0”的基矢代表 lo，
	// 每个代表对应一对 (lo, hi=lo|tbit)。
	pairs := make([][2]int, 0, dim/2)
	for idx := 0; idx < dim; idx++ {
		if idx&tbit != 0 {
			continue
		}
		if cmask != 0 && idx&cmask != cmask {
			continue
		}
		pairs = append(pairs, [2]int{idx, idx | tbit})
	}

	// 先对行作用 U（rho <- U rho）：列固定，取该列在 pairs 上的 2 维子向量。
	tmp := make([]complex128, len(rho))
	copy(tmp, rho)
	for _, p := range pairs {
		for j := 0; j < dim; j++ {
			a0 := tmp[p[0]*dim+j]
			a1 := tmp[p[1]*dim+j]
			rho[p[0]*dim+j] = m[0]*a0 + m[2]*a1
			rho[p[1]*dim+j] = m[1]*a0 + m[3]*a1
		}
	}
	// 再对列作用 U*（rho <- rho U†）：行固定。
	copy(tmp, rho)
	for i := 0; i < dim; i++ {
		for _, p := range pairs {
			a0 := tmp[i*dim+p[0]]
			a1 := tmp[i*dim+p[1]]
			rho[i*dim+p[0]] = a0*conjM(m[0]) + a1*conjM(m[1])
			rho[i*dim+p[1]] = a0*conjM(m[2]) + a1*conjM(m[3])
		}
	}
	return rho
}
