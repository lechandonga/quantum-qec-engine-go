package noise

import "math"

// applyPauli 原地将 Pauli 算符 P（X/Y/Z/I 的逐比特张量积）作用在
// 2^n 维向量上。xmask 为 X 或 Y 的比特位，zmask 为 Z 或 Y 的比特位。
// 全局相位无关（本包中 Kraus 算符最多只差符号/虚单位的全局因子）。
func applyPauli(amp []complex128, n int, xmask, zmask int) {
	if xmask == 0 && zmask == 0 {
		return
	}
	ySign := func(idx int) int {
		s := 1
		for q := 0; q < n; q++ {
			bit := 1 << q
			if zmask&bit != 0 && xmask&bit != 0 && idx&bit == 0 {
				s = -s
			}
		}
		return s
	}
	if xmask == 0 {
		for idx := range amp {
			if zmask&idx != 0 {
				amp[idx] = -amp[idx]
			}
		}
		return
	}
	out := make([]complex128, len(amp))
	for idx, a := range amp {
		j := idx ^ xmask
		s := ySign(idx)
		out[j] += complex(float64(s), 0) * a
	}
	copy(amp, out)
}

func conjM(z complex128) complex128 {
	return complex(real(z), -imag(z))
}

// traceDensity 返回密度矩阵的迹。
func traceDensity(rho []complex128, dim int) float64 {
	t := 0.0
	for i := 0; i < dim; i++ {
		t += real(rho[i*dim+i])
	}
	return t
}

func densityFromPure(amp []complex128) []complex128 {
	dim := len(amp)
	rho := make([]complex128, dim*dim)
	for i := 0; i < dim; i++ {
		for j := 0; j < dim; j++ {
			rho[i*dim+j] = amp[i] * conjM(amp[j])
		}
	}
	return rho
}

func densityToProbabilities(rho []complex128, dim int) []float64 {
	p := make([]float64, dim)
	for i := 0; i < dim; i++ {
		v := real(rho[i*dim+i])
		if math.Abs(v) < 1e-15 {
			v = 0
		}
		p[i] = v
	}
	return p
}
