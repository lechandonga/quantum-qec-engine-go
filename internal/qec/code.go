// Package qec implements stabiliser quantum error correction for the
// generalised Shor code (a CSS concatenation of a bit-flip repetition code
// inside a phase-flip repetition code), syndrome extraction, minimum-weight
// decoding, and logical-error statistics.
//
// For odd code distance d the code uses n=d*d physical qubits arranged as d
// blocks of d qubits:
//
//   - Z-checks: adjacent Z Z pairs inside each block (d-1 per block, d(d-1)
//     total), which localise X/bit errors to positions.
//   - X-checks: the all-X operator over one block times the all-X operator
//     over the next block (d-1 total), which localise Z/phase errors to
//     blocks.
//
// The logical operators are X_L = X on every qubit of any one block and
// Z_L = Z on any single physical qubit. The decoder corrects up to
// (d-1)/2 bit errors per block and (d-1)/2 errored blocks at the phase
// level; heavier errors are explicitly flagged rather than reported as
// silent successes.
package qec

import (
	"errors"
	"fmt"
)

var (
	// ErrInvalidDistance is returned when the code distance is not an odd
	// positive integer >= 3.
	ErrInvalidDistance = errors.New("qec: code distance must be an odd integer >= 3")
	// ErrInvalidProbability is returned when error rates are outside [0,1].
	ErrInvalidProbability = errors.New("qec: error probability out of range")
	// ErrResourceLimit is returned before any work starts when the code size
	// or sample count exceeds the configured budget.
	ErrResourceLimit = errors.New("qec: resource limit exceeded")
)

// Code describes one fixed generalised Shor-code instance.
type Code struct {
	d      int // odd code distance
	blocks int // number of blocks (== d)
	size   int // qubits per block (== d)
	n      int // total physical qubits == d*d
}

// NewCode constructs a distance-d generalised Shor code.
func NewCode(d int) (*Code, error) {
	if d < 3 || d%2 == 0 {
		return nil, fmt.Errorf("%w: got %d", ErrInvalidDistance, d)
	}
	if d*d <= 0 {
		return nil, fmt.Errorf("%w: code size overflow", ErrResourceLimit)
	}
	return &Code{d: d, blocks: d, size: d, n: d * d}, nil
}

// Distance returns the code distance.
func (c *Code) Distance() int { return c.d }

// NumQubits returns the physical qubit count d*d.
func (c *Code) NumQubits() int { return c.n }

// Qubit returns the flat index of (block, positionInBlock).
func (c *Code) Qubit(block, pos int) int { return block*c.size + pos }

// BlockOf returns the block containing physical qubit q.
func (c *Code) BlockOf(q int) int { return q / c.size }

// PosInBlock returns the within-block position of physical qubit q.
func (c *Code) PosInBlock(q int) int { return q % c.size }

// CorrectableWeight is (d-1)/2.
func (c *Code) CorrectableWeight() int { return (c.d - 1) / 2 }
