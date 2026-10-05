package qcircuit

import (
	"math"
	"math/cmplx"
)

// Gate is an instruction in a circuit. Qubits holds target qubits,
// Controls holds control qubits (gate fires only when all controls are 1).
// Param is the rotation angle in radians for parameterised gates.
type Gate struct {
	Name     string
	Qubits   []int
	Controls []int
	Param    float64
}

type gateSpec struct {
	targets  int  // required number of target qubits
	controls int  // required number of controls (-1 means any number allowed)
	rotParam bool // gate takes an angle parameter
}

// gateTable maps the canonical upper-case gate name to its specification.
var gateTable = map[string]gateSpec{
	"I":    {targets: 1, controls: -1},
	"X":    {targets: 1, controls: -1},
	"Y":    {targets: 1, controls: -1},
	"Z":    {targets: 1, controls: -1},
	"H":    {targets: 1, controls: -1},
	"S":    {targets: 1, controls: 0},
	"SDG":  {targets: 1, controls: 0},
	"T":    {targets: 1, controls: 0},
	"TDG":  {targets: 1, controls: 0},
	"RX":   {targets: 1, controls: -1, rotParam: true},
	"RY":   {targets: 1, controls: -1, rotParam: true},
	"RZ":   {targets: 1, controls: -1, rotParam: true},
	"CNOT": {targets: 1, controls: 1}, // CX alias
	"CX":   {targets: 1, controls: 1},
	"CY":   {targets: 1, controls: 1},
	"CZ":   {targets: 1, controls: 1},
	"CH":   {targets: 1, controls: 1},
	"CRX":  {targets: 1, controls: 1, rotParam: true},
	"CRY":  {targets: 1, controls: 1, rotParam: true},
	"CRZ":  {targets: 1, controls: 1, rotParam: true},
	"SWAP": {targets: 2, controls: -1},
}

func (g Gate) validate(n int) error {
	spec, ok := gateTable[g.Name]
	if !ok {
		return ErrUnknownGate
	}
	if spec.rotParam && (math.IsNaN(g.Param) || math.IsInf(g.Param, 0)) {
		return ErrInvalidParameter
	}
	if len(g.Qubits) != spec.targets {
		return ErrControlConflict
	}
	switch {
	case spec.controls == -1:
		// arbitrary (including zero) controls allowed
	case len(g.Controls) != spec.controls:
		return ErrControlConflict
	}
	seen := make(map[int]bool, len(g.Qubits)+len(g.Controls))
	for _, q := range g.Qubits {
		if q < 0 || q >= n {
			return ErrQubitOutOfRange
		}
		if seen[q] {
			return ErrControlConflict
		}
		seen[q] = true
	}
	for _, q := range g.Controls {
		if q < 0 || q >= n {
			return ErrQubitOutOfRange
		}
		if seen[q] {
			return ErrControlConflict
		}
		seen[q] = true
	}
	return nil
}

// singleMatrix returns the 2x2 unitary (rows/cols indexed by local bit)
// for a single-qubit gate.
func singleMatrix(name string, theta float64) [4]complex128 {
	c := func(x float64) complex128 { return complex(x, 0) }
	switch name {
	case "I":
		return [4]complex128{1, 0, 0, 1}
	case "X":
		return [4]complex128{0, 1, 1, 0}
	case "Y":
		return [4]complex128{0, complex128(-1i), complex128(1i), 0}
	case "Z":
		return [4]complex128{1, 0, 0, -1}
	case "H":
		v := complex(1/math.Sqrt2, 0)
		return [4]complex128{v, v, v, -v}
	case "S":
		return [4]complex128{1, 0, 0, complex128(1i)}
	case "SDG":
		return [4]complex128{1, 0, 0, complex128(-1i)}
	case "T":
		return [4]complex128{1, 0, 0, cmplx.Exp(1i * math.Pi / 4)}
	case "TDG":
		return [4]complex128{1, 0, 0, cmplx.Exp(-1i * math.Pi / 4)}
	case "RX":
		cos := c(math.Cos(theta / 2))
		sin := complex(0, -math.Sin(theta/2))
		return [4]complex128{cos, sin, sin, cos}
	case "RY":
		cos := c(math.Cos(theta / 2))
		sin := c(math.Sin(theta / 2))
		return [4]complex128{cos, -sin, sin, cos}
	case "RZ":
		em := cmplx.Exp(complex(0, -theta/2))
		ep := cmplx.Exp(complex(0, theta/2))
		return [4]complex128{em, 0, 0, ep}
	}
	return [4]complex128{}
}
