// Package qcircuit implements a pure-Go state-vector quantum circuit simulator.
//
// Bit ordering convention: qubit q corresponds to bit position 1<<q in the
// basis index, so basis state |q1 q0> with q0 the least-significant qubit is
// stored at index (q1<<1)|q0.
package qcircuit

import "errors"

// Distinguishable sentinel errors. Callers can use errors.Is to branch on the
// exact validation failure.
var (
	// ErrQubitOutOfRange is returned when a gate or measurement references a
	// qubit index outside [0, NumQubits()).
	ErrQubitOutOfRange = errors.New("qcircuit: qubit index out of range")
	// ErrUnknownGate is returned when the gate name is not registered.
	ErrUnknownGate = errors.New("qcircuit: unknown gate")
	// ErrControlConflict is returned when controls overlap each other or the
	// target set, or when the arity of controls/targets is wrong.
	ErrControlConflict = errors.New("qcircuit: control/target conflict or bad arity")
	// ErrInvalidParameter is returned when a rotation parameter is NaN/+/-Inf
	// or an otherwise invalid parameter is supplied.
	ErrInvalidParameter = errors.New("qcircuit: invalid gate parameter")
	// ErrInvalidState is returned when a state vector has the wrong length or
	// is not normalised within tolerance.
	ErrInvalidState = errors.New("qcircuit: invalid state vector")
	// ErrInvalidPauli is returned when an expectation Pauli string is malformed.
	ErrInvalidPauli = errors.New("qcircuit: invalid Pauli string")
	// ErrResourceLimit is returned before any heavy allocation when the
	// configured resource limits would be exceeded.
	ErrResourceLimit = errors.New("qcircuit: resource limit exceeded")
)
