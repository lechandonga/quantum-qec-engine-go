package noise

import (
	"fmt"
	"math"
)

// Density is a small exact density matrix in the computational basis, stored
// row-major as rho[i*dim+j]. It exists for exact, tolerance-based checks of
// the stochastic sampler and is not intended for large systems.
type Density struct {
	n   int
	dim int
	rho []complex128
}

// NewDensity rejects systems larger than the exact-reference budget (2^n x
// 2^n complex entries); n=5 is the supported default.
func NewDensity(n int) (*Density, error) {
	const maxExactQubits = 5
	if n < 1 || n > maxExactQubits {
		return nil, fmt.Errorf("%w: exact density reference supports 1..%d qubits, got %d",
			ErrResourceLimit, maxExactQubits, n)
	}
	dim := 1 << n
	return &Density{n: n, dim: dim, rho: make([]complex128, dim*dim)}, nil
}

// FromPure builds rho = |psi><psi| from an amplitude vector.
func FromPure(amp []complex128) (*Density, error) {
	n := 0
	for dim := 1; dim < len(amp); dim <<= 1 {
		n++
	}
	if len(amp) != 1<<n {
		return nil, fmt.Errorf("%w: amplitude length %d is not a power of two", ErrResourceLimit, len(amp))
	}
	d, err := NewDensity(n)
	if err != nil {
		return nil, err
	}
	for i := 0; i < d.dim; i++ {
		for j := 0; j < d.dim; j++ {
			d.rho[i*d.dim+j] = amp[i] * conjC(amp[j])
		}
	}
	return d, nil
}

func conjC(z complex128) complex128 { return complex(real(z), -imag(z)) }

// Trace returns Tr(rho).
func (d *Density) Trace() complex128 {
	var tr complex128
	for i := 0; i < d.dim; i++ {
		tr += d.rho[i*d.dim+i]
	}
	return tr
}

// Probabilities returns the exact computational-basis diagonal.
func (d *Density) Probabilities() []float64 {
	out := make([]float64, d.dim)
	var total float64
	for i := 0; i < d.dim; i++ {
		out[i] = real(d.rho[i*d.dim+i])
		total += out[i]
	}
	if math.Abs(total-1) < 1e-12 && total > 0 {
		for i := range out {
			out[i] /= total
		}
	}
	return out
}

// PauliExpectation returns Tr(P rho) for an I/X/Y/Z tensor product.
func (d *Density) PauliExpectation(pauli []byte) (float64, error) {
	if len(pauli) != d.n {
		return 0, fmt.Errorf("noise: want %d Pauli letters, got %d", d.n, len(pauli))
	}
	for _, ch := range pauli {
		switch ch {
		case 'I', 'i', 'X', 'x', 'Y', 'y', 'Z', 'z':
		default:
			return 0, fmt.Errorf("noise: invalid Pauli letter %q", ch)
		}
	}
	var sum complex128
	for j := 0; j < d.dim; j++ {
		k, phase := pauliApply(pauli, j)
		// <j|P|k> is nonzero only when k is the image; Tr(P rho)=sum_j P_{k,j} rho_{j,k}
		sum += phase * d.rho[j*d.dim+k]
	}
	return real(sum), nil
}

// pauliApply maps basis state j under the Pauli string and returns
// (image index, matrix coefficient).
func pauliApply(pauli []byte, j int) (int, complex128) {
	var phase complex128 = 1
	for q, ch := range pauli {
		bit := (j >> uint(q)) & 1
		switch ch {
		case 'X', 'x':
			j ^= 1 << uint(q)
		case 'Y', 'y':
			j ^= 1 << uint(q)
			if bit == 0 {
				phase *= 1i // Y|0> = i|1>
			} else {
				phase *= -1i // Y|1> = -i|0>
			}
		case 'Z', 'z':
			if bit == 1 {
				phase *= -1
			}
		}
	}
	return j, phase
}

// krausTerm applies E rho E† for a single Pauli Kraus operator on qubit q
// (letter I/X/Y/Z) scaled by sqrt(p).
func (d *Density) krausTerm(q int, letter byte, scale complex128, out []complex128) {
	dim := d.dim
	for a := 0; a < dim; a++ {
		ai, ac := singlePauli(letter, q, a)
		for b := 0; b < dim; b++ {
			bi, bc := singlePauli(letter, q, b)
			out[ai*dim+bi] += scale * ac * d.rho[a*dim+b] * scale * conjC(bc)
		}
	}
}

// singlePauli applies one single-qubit Pauli on qubit q to basis index j.
func singlePauli(letter byte, q, j int) (int, complex128) {
	bit := (j >> uint(q)) & 1
	switch letter {
	case 'I', 'i':
		return j, 1
	case 'X', 'x':
		return j ^ (1 << uint(q)), 1
	case 'Y', 'y':
		j ^= 1 << uint(q)
		if bit == 0 {
			return j, 1i
		}
		return j, -1i
	case 'Z', 'z':
		if bit == 1 {
			return j, -1
		}
		return j, 1
	}
	return j, 1
}

// applyChannel replaces rho by sum_k E_k rho E_k† for the Kraus letters given
// with probabilities. letters[i] is applied with probability probs[i].
func (d *Density) applyChannel(q int, letters []byte, probs []float64) {
	out := make([]complex128, len(d.rho))
	for i, p := range probs {
		if p == 0 {
			continue
		}
		d.krausTerm(q, letters[i], complex(math.Sqrt(p), 0), out)
	}
	d.rho = out
}

// ApplyExact applies the full noise model exactly to the density matrix.
func (d *Density) ApplyExact(m Model) error {
	if err := m.Validate(d.n); err != nil {
		return err
	}
	for _, q := range m.targetQubits(d.n) {
		if m.BitFlip > 0 {
			d.applyChannel(q, []byte{'I', 'X'}, []float64{1 - m.BitFlip, m.BitFlip})
		}
		if m.PhaseFlip > 0 {
			d.applyChannel(q, []byte{'I', 'Z'}, []float64{1 - m.PhaseFlip, m.PhaseFlip})
		}
		if m.Depolarising > 0 {
			p := m.Depolarising
			d.applyChannel(q, []byte{'I', 'X', 'Y', 'Z'},
				[]float64{1 - p, p / 3, p / 3, p / 3})
		}
	}
	return nil
}

// FidelityWithPure returns <psi|rho|psi> for a pure reference state.
func (d *Density) FidelityWithPure(amp []complex128) (float64, error) {
	if len(amp) != d.dim {
		return 0, fmt.Errorf("noise: amplitude length %d != dimension %d", len(amp), d.dim)
	}
	var sum complex128
	for i := 0; i < d.dim; i++ {
		for j := 0; j < d.dim; j++ {
			sum += conjC(amp[i]) * d.rho[i*d.dim+j] * amp[j]
		}
	}
	return real(sum), nil
}
