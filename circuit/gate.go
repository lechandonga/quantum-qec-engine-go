package circuit

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Gate 描述一条量子门指令。Controls 为控制位集合（全部为 |1> 时生效），
// Target 为作用比特；Param 供带参数旋转门使用（弧度）。
type Gate struct {
	Name     string
	Target   int
	Controls []int
	Param    float64
}

func (g Gate) String() string {
	if len(g.Controls) > 0 {
		return fmt.Sprintf("C%v-%s(target=%d, param=%v)", g.Controls, g.Name, g.Target, g.Param)
	}
	return fmt.Sprintf("%s(target=%d, param=%v)", g.Name, g.Target, g.Param)
}

// 固定的 2x2 单比特门矩阵（列主序下标 m[row+2*col]）。
var fixedGates = map[string][4]complex128{
	"I": {1, 0, 0, 1},
	"X": {0, 1, 1, 0},
	"Y": {0, complex(0, -1), complex(0, 1), 0},
	"Z": {1, 0, 0, -1},
	"H": {1 / math.Sqrt2, 1 / math.Sqrt2, 1 / math.Sqrt2, -1 / math.Sqrt2},
	"S": {1, 0, 0, complex(0, 1)},
	"T": {1, 0, 0, cmplx.Exp(complex(0, math.Pi/4))},
}

// IsFixedGate 报告名称是否为无参固定门。
func IsFixedGate(name string) bool {
	_, ok := fixedGates[name]
	return ok
}

// IsRotationGate 报告名称是否为带参数旋转门。
func IsRotationGate(name string) bool {
	switch name {
	case "RX", "RY", "RZ", "PHASE":
		return true
	default:
		return false
	}
}

// rotationMatrix 返回绕坐标轴旋转或相位门的矩阵。角度为弧度。
func rotationMatrix(name string, theta float64) [4]complex128 {
	c := complex(math.Cos(theta/2), 0)
	s := complex(math.Sin(theta/2), 0)
	i := complex(0, 1)
	switch name {
	case "RX":
		return [4]complex128{c, -i * s, -i * s, c}
	case "RY":
		return [4]complex128{c, -s, s, c}
	case "RZ":
		em := cmplx.Exp(-i * complex(theta/2, 0))
		ep := cmplx.Exp(i * complex(theta/2, 0))
		return [4]complex128{em, 0, 0, ep}
	case "PHASE":
		return [4]complex128{1, 0, 0, cmplx.Exp(i * complex(theta, 0))}
	default:
		panic("circuit: not a rotation gate")
	}
}

func gateMatrix(g Gate) [4]complex128 {
	if m, ok := fixedGates[g.Name]; ok {
		return m
	}
	if IsRotationGate(g.Name) {
		if math.IsNaN(g.Param) || math.IsInf(g.Param, 0) {
			panic(ErrInvalidParam)
		}
		return rotationMatrix(g.Name, g.Param)
	}
	panic(ErrUnknownGate)
}

// GateMatrix 导出单比特门的 2x2 复矩阵（列主序），供密度矩阵演化复用。
func GateMatrix(g Gate) ([4]complex128, error) {
	if !IsFixedGate(g.Name) && !IsRotationGate(g.Name) {
		return [4]complex128{}, ErrUnknownGate
	}
	if IsRotationGate(g.Name) && (math.IsNaN(g.Param) || math.IsInf(g.Param, 0)) {
		return [4]complex128{}, ErrInvalidParam
	}
	return gateMatrix(g), nil
}
