package qec

import "math/rand/v2"

// splitmix64 与 noise 包保持同一确定性派生约定：
// 第 shot 次试验的随机种子只依赖 (主种子, shot)，与协程调度无关。
func splitmix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

func newChildRNG(seed uint64, shot int) *rand.Rand {
	s := splitmix64(seed ^ splitmix64(uint64(shot)+0x128ad98b5f3c2e61))
	return rand.New(rand.NewPCG(s, splitmix64(s^0x3c6ef372fe94f82b)))
}
