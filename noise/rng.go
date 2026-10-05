package noise

import "math/rand/v2"

// splitmix64 是确定性的 64 位混合函数。主种子经过它派生出每一个样本
// 的独立子种子：子种子序列只取决于主种子和样本编号，与工作协程数量、
// 调度顺序完全无关，从而保证并行抽样的逐位可复现性。
func splitmix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

func newChildRNG(seed uint64, shot int) *rand.Rand {
	s := splitmix64(seed ^ splitmix64(uint64(shot)+0xd1b54a32d192ed03))
	return rand.New(rand.NewPCG(s, splitmix64(s^0x243f6a8885a308d3)))
}
