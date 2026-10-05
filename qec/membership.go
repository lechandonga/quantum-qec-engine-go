package qec

// cosetReduction 返回 e（已叠加 X 恢复、不含 Z 恢复）的 GF2 稳定子
// 陪集约化结果；返回 0 表示完全属于稳定子空间。
func cosetReduction(c Code, e, xCorr PauliError) uint {
	n := c.N
	r := PauliError{Xmask: e.Xmask ^ xCorr.Xmask, Zmask: e.Zmask}
	vec := uint(r.Xmask) | uint(r.Zmask)<<uint(n)
	v, _ := gf2Reduce(c.stabilizerRows(), vec, 2*n)
	return v
}
