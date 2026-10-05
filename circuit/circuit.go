package circuit

// 本文件只依赖标准库的内建运算。

// Circuit 是有序的量子门序列，全部作用在 N 个比特上。
type Circuit struct {
	n     int
	gates []Gate
}

// New 创建一个 n 比特的空线路。n 必须为正。
func New(n int) (*Circuit, error) {
	if n <= 0 {
		return nil, ErrInvalidQubit
	}
	if n > 63 {
		return nil, ErrInvalidQubit
	}
	return &Circuit{n: n}, nil
}

// N 返回比特数。
func (c *Circuit) N() int { return c.n }

// Len 返回门数量。
func (c *Circuit) Len() int { return len(c.gates) }

// Gates 返回门序列的副本。
func (c *Circuit) Gates() []Gate {
	out := make([]Gate, len(c.gates))
	copy(out, c.gates)
	return out
}

func (c *Circuit) validate(g Gate) error {
	if !IsFixedGate(g.Name) && !IsRotationGate(g.Name) {
		return ErrUnknownGate
	}
	if IsRotationGate(g.Name) && (isBadFloat(g.Param)) {
		return ErrInvalidParam
	}
	if g.Target < 0 || g.Target >= c.n {
		return ErrInvalidQubit
	}
	seen := make(map[int]struct{}, len(g.Controls)+1)
	seen[g.Target] = struct{}{}
	for _, q := range g.Controls {
		if q < 0 || q >= c.n {
			return ErrInvalidQubit
		}
		if _, dup := seen[q]; dup {
			return ErrControlConflict
		}
		seen[q] = struct{}{}
	}
	return nil
}

func isBadFloat(x float64) bool { return x != x || x > 1e308 || x < -1e308 }

// Add 追加一条门指令并在追加前完成全部校验。
func (c *Circuit) Add(g Gate) error {
	if err := c.validate(g); err != nil {
		return err
	}
	g.Controls = append([]int(nil), g.Controls...)
	c.gates = append(c.gates, g)
	return nil
}

// MustAdd 与 Add 相同，但校验失败时 panic，适合测试与字面量构造。
func (c *Circuit) MustAdd(g Gate) {
	if err := c.Add(g); err != nil {
		panic(err)
	}
}

// 便捷构造方法 ---------------------------------------------------------------

func (c *Circuit) X(target int, controls ...int) error {
	return c.Add(Gate{Name: "X", Target: target, Controls: controls})
}
func (c *Circuit) Y(target int, controls ...int) error {
	return c.Add(Gate{Name: "Y", Target: target, Controls: controls})
}
func (c *Circuit) Z(target int, controls ...int) error {
	return c.Add(Gate{Name: "Z", Target: target, Controls: controls})
}
func (c *Circuit) H(target int) error { return c.Add(Gate{Name: "H", Target: target}) }
func (c *Circuit) S(target int) error { return c.Add(Gate{Name: "S", Target: target}) }
func (c *Circuit) T(target int) error { return c.Add(Gate{Name: "T", Target: target}) }

func (c *Circuit) RX(target int, theta float64, controls ...int) error {
	return c.Add(Gate{Name: "RX", Target: target, Param: theta, Controls: controls})
}
func (c *Circuit) RY(target int, theta float64, controls ...int) error {
	return c.Add(Gate{Name: "RY", Target: target, Param: theta, Controls: controls})
}
func (c *Circuit) RZ(target int, theta float64, controls ...int) error {
	return c.Add(Gate{Name: "RZ", Target: target, Param: theta, Controls: controls})
}
func (c *Circuit) Phase(target int, theta float64) error {
	return c.Add(Gate{Name: "PHASE", Target: target, Param: theta})
}

// Apply 将整条线路作用到纯态向量上（原地修改）。
func (c *Circuit) Apply(amp []complex128) error {
	if len(amp) != 1<<c.n {
		return ErrInvalidState
	}
	for _, g := range c.gates {
		applyGate(amp, c.n, g)
	}
	return nil
}

func applyGate(amp []complex128, n int, g Gate) {
	m := gateMatrix(g)
	tbit := 1 << g.Target
	cmask := 0
	for _, q := range g.Controls {
		cmask |= 1 << q
	}
	for idx := 0; idx < len(amp); idx++ {
		if idx&tbit != 0 {
			continue // 每对只处理 target=0 的基矢
		}
		if cmask != 0 && idx&cmask != cmask {
			continue
		}
		i0 := idx
		i1 := idx | tbit
		a0 := amp[i0]
		a1 := amp[i1]
		amp[i0] = m[0]*a0 + m[2]*a1
		amp[i1] = m[1]*a0 + m[3]*a1
	}
}
