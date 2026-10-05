package qec

import (
	"fmt"
)

// ErrorSet is a Pauli fault pattern on the physical qubits: X[q] means an
// X (bit) fault, Z[q] means a Z (phase) fault, both means a Y fault.
type ErrorSet struct {
	X []bool
	Z []bool
}

// NewErrorSet returns an empty fault pattern for code c.
func NewErrorSet(c *Code) *ErrorSet {
	return &ErrorSet{X: make([]bool, c.n), Z: make([]bool, c.n)}
}

// Weight returns the number of physical qubits carrying a nontrivial fault.
func (e *ErrorSet) Weight() int {
	w := 0
	for q := range e.X {
		if e.X[q] || e.Z[q] {
			w++
		}
	}
	return w
}

// Syndrome holds the extracted check outcomes.
//
// BitSyndrome[b][i] is the Z-check Z_i Z_{i+1} inside block b
// (i in 0..d-2): the parity of bit-flip (X) faults at positions i,i+1.
//
// PhaseInner[b][i] is the X-check X_i X_{i+1} inside block b
// (i in 0..d-2): the parity of phase-flip (Z) faults at positions i,i+1.
//
// PhaseOuter[i] is the X-check (all X in block i)(all X in block i+1)
// (i in 0..d-2): parity of the phase-fault parity vectors of blocks i,i+1.
type Syndrome struct {
	BitSyndrome [][]bool // [block][inner edge] X-fault localisation
	PhaseInner  [][]bool // [block][inner edge] Z-fault localisation
	PhaseOuter  []bool   // [block-edge] block-parity defects
}

// Extract measures every stabiliser check against the (known, for analysis)
// Pauli fault pattern. It is the ideal-syndrome-extraction model: check
// outcomes are computed exactly, with no extra measurement faults.
func (c *Code) Extract(e *ErrorSet) (*Syndrome, error) {
	if e == nil || len(e.X) != c.n || len(e.Z) != c.n {
		return nil, fmt.Errorf("qec: error set does not match code size %d", c.n)
	}
	syn := &Syndrome{
		BitSyndrome: make([][]bool, c.blocks),
		PhaseInner:  make([][]bool, c.blocks),
		PhaseOuter:  make([]bool, c.blocks-1),
	}
	blockPhaseParity := make([]bool, c.blocks)
	for b := 0; b < c.blocks; b++ {
		syn.BitSyndrome[b] = make([]bool, c.size-1)
		syn.PhaseInner[b] = make([]bool, c.size-1)
		for i := 0; i < c.size-1; i++ {
			q0 := c.Qubit(b, i)
			q1 := c.Qubit(b, i+1)
			syn.BitSyndrome[b][i] = e.X[q0] != e.X[q1]
			syn.PhaseInner[b][i] = e.Z[q0] != e.Z[q1]
		}
		var parity bool
		for i := 0; i < c.size; i++ {
			parity = parity != e.Z[c.Qubit(b, i)]
		}
		blockPhaseParity[b] = parity
	}
	for i := 0; i < c.blocks-1; i++ {
		syn.PhaseOuter[i] = blockPhaseParity[i] != blockPhaseParity[i+1]
	}
	return syn, nil
}
