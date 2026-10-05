package qcircuit

import (
	"fmt"
)

// Limits caps the resources a circuit (and its simulation) may consume.
// A zero field means "unlimited".
type Limits struct {
	MaxQubits int
	MaxGates  int
}

// DefaultLimits rejects obviously excessive inputs before allocation.
var DefaultLimits = Limits{MaxQubits: 24, MaxGates: 200_000}

// Circuit is an ordered list of gates over a fixed number of qubits.
type Circuit struct {
	n      int
	gates  []Gate
	limits Limits
}

// NewCircuit creates an n-qubit circuit. n must be positive and within limits.
func NewCircuit(n int) (*Circuit, error) {
	return NewCircuitWithLimits(n, DefaultLimits)
}

// NewCircuitWithLimits is NewCircuit with an explicit resource budget.
func NewCircuitWithLimits(n int, limits Limits) (*Circuit, error) {
	if n <= 0 {
		return nil, fmt.Errorf("%w: qubit count must be positive, got %d", ErrQubitOutOfRange, n)
	}
	if limits.MaxQubits > 0 && n > limits.MaxQubits {
		return nil, fmt.Errorf("%w: %d qubits > max %d", ErrResourceLimit, n, limits.MaxQubits)
	}
	if n > 63 {
		return nil, fmt.Errorf("%w: state-vector simulator supports at most 63 qubits", ErrResourceLimit)
	}
	return &Circuit{n: n, limits: limits}, nil
}

// NumQubits reports the circuit width.
func (c *Circuit) NumQubits() int { return c.n }

// NumGates reports the number of appended gates.
func (c *Circuit) NumGates() int { return len(c.gates) }

// Add appends a generic gate and validates it immediately.
func (c *Circuit) Add(name string, targets, controls []int, param float64) error {
	if c.limits.MaxGates > 0 && len(c.gates) >= c.limits.MaxGates {
		return fmt.Errorf("%w: gate count > max %d", ErrResourceLimit, c.limits.MaxGates)
	}
	t := append([]int(nil), targets...)
	ctl := append([]int(nil), controls...)
	g := Gate{Name: name, Qubits: t, Controls: ctl, Param: param}
	if err := g.validate(c.n); err != nil {
		return fmt.Errorf("add %s: %w", name, err)
	}
	c.gates = append(c.gates, g)
	return nil
}

// Gates returns a defensive copy of the gate list.
func (c *Circuit) Gates() []Gate {
	out := make([]Gate, len(c.gates))
	copy(out, c.gates)
	return out
}

// X appends a (optionally controlled) Pauli-X gate.
func (c *Circuit) X(q int, controls ...int) error { return c.Add("X", []int{q}, controls, 0) }

// Y appends a (optionally controlled) Pauli-Y gate.
func (c *Circuit) Y(q int, controls ...int) error { return c.Add("Y", []int{q}, controls, 0) }

// Z appends a (optionally controlled) Pauli-Z gate.
func (c *Circuit) Z(q int, controls ...int) error { return c.Add("Z", []int{q}, controls, 0) }

// H appends a (optionally controlled) Hadamard gate.
func (c *Circuit) H(q int, controls ...int) error { return c.Add("H", []int{q}, controls, 0) }

// S appends the phase gate S.
func (c *Circuit) S(q int) error { return c.Add("S", []int{q}, nil, 0) }

// Sdg appends the adjoint phase gate.
func (c *Circuit) Sdg(q int) error { return c.Add("SDG", []int{q}, nil, 0) }

// T appends the T gate.
func (c *Circuit) T(q int) error { return c.Add("T", []int{q}, nil, 0) }

// Tdg appends the adjoint T gate.
func (c *Circuit) Tdg(q int) error { return c.Add("TDG", []int{q}, nil, 0) }

// RX appends an X-axis rotation by theta radians.
func (c *Circuit) RX(q int, theta float64) error {
	return c.Add("RX", []int{q}, nil, theta)
}

// RY appends a Y-axis rotation by theta radians.
func (c *Circuit) RY(q int, theta float64) error {
	return c.Add("RY", []int{q}, nil, theta)
}

// RZ appends a Z-axis rotation by theta radians.
func (c *Circuit) RZ(q int, theta float64) error {
	return c.Add("RZ", []int{q}, nil, theta)
}

// CNOT appends a controlled-X (CX) gate.
func (c *Circuit) CNOT(control, target int) error {
	return c.Add("CNOT", []int{target}, []int{control}, 0)
}

// CZ appends a controlled-Z gate.
func (c *Circuit) CZ(control, target int) error {
	return c.Add("CZ", []int{target}, []int{control}, 0)
}

// SWAP appends a swap gate.
func (c *Circuit) SWAP(a, b int) error {
	return c.Add("SWAP", []int{a, b}, nil, 0)
}
