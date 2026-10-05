package qec

// DecodeOutcome reports the decoder's conclusion for one fault pattern.
type DecodeOutcome int

const (
	// Corrected means a unique minimum-weight correction removed all logical
	// action; the logical state is recovered.
	Corrected DecodeOutcome = iota
	// PossibleLogicalFlip means the syndrome had degenerate minimum-weight
	// explanations, or after correction a residual logical action remains
	// plausible; the logical state may have flipped.
	PossibleLogicalFlip
	// LogicalFlip means the post-correction residual provably contains a
	// logical operator (the error weight was beyond the guaranteed
	// correction radius).
	LogicalFlip
)

func (o DecodeOutcome) String() string {
	switch o {
	case Corrected:
		return "corrected"
	case PossibleLogicalFlip:
		return "possible_logical_flip"
	case LogicalFlip:
		return "logical_flip"
	}
	return "unknown"
}

// DecodeResult holds the correction (as a Pauli fault pattern), the decoder
// verdict, and the residual fault pattern correction applied to the injected
// errors (modulo stabilisers as detected by logical-parity tests).
type DecodeResult struct {
	Correction        *ErrorSet
	Outcome           DecodeOutcome
	Residual          *ErrorSet
	BitLogicalFlip    bool
	PhaseLogicalFlip  bool
	AmbiguousSyndrome bool
}

// Decode extracts the syndrome of e, finds minimum-weight corrections,
// applies them, and classifies the residual. Correction is deterministic for
// a fixed syndrome and code; ties are reported as PossibleLogicalFlip.
func (c *Code) Decode(e *ErrorSet) (*DecodeResult, error) {
	syn, err := c.Extract(e)
	if err != nil {
		return nil, err
	}
	corr := NewErrorSet(c)
	ambiguous := false

	// Bit level: independent repetition-code instance per block.
	for b := 0; b < c.blocks; b++ {
		chain, amb, err := decodeRepetition(syn.BitSyndrome[b])
		if err != nil {
			return nil, err
		}
		if amb {
			ambiguous = true
		}
		for i, flip := range chain {
			if flip {
				corr.X[c.Qubit(b, i)] = true
			}
		}
	}

	// Phase level is two nested repetition codes:
	//  1. each block is matched internally on PhaseInner edges, localising
	//     Z faults to positions (same line matching as the bit level);
	//  2. the resulting per-block residual parity feeds the outer block
	//     repetition code on PhaseOuter edges; a block correction adds one
	//     Z at the position chosen by the inner chain.
	innerChains := make([][]bool, c.blocks)
	for b := 0; b < c.blocks; b++ {
		chain, amb, err := decodeRepetition(syn.PhaseInner[b])
		if err != nil {
			return nil, err
		}
		if amb {
			ambiguous = true
		}
		innerChains[b] = chain
	}
	outer, amb, err := decodeRepetition(syn.PhaseOuter)
	if err != nil {
		return nil, err
	}
	if amb {
		ambiguous = true
	}
	for b, chain := range innerChains {
		pos := -1
		for i, flip := range chain {
			if flip {
				corr.Z[c.Qubit(b, i)] = true
				if pos < 0 {
					pos = i
				}
			}
		}
		if outer[b] {
			// Outer block parity correction: one Z at the inner chain's
			// representative (position 0 if the inner chain is empty).
			if pos < 0 {
				pos = 0
			}
			corr.Z[c.Qubit(b, pos)] = !corr.Z[c.Qubit(b, pos)]
		}
	}

	residual := NewErrorSet(c)
	for q := 0; q < c.n; q++ {
		residual.X[q] = e.X[q] != corr.X[q]
		residual.Z[q] = e.Z[q] != corr.Z[q]
	}

	bitFlip := c.hasBitLogical(residual)
	phaseFlip := c.hasPhaseLogical(residual)

	outcome := Corrected
	switch {
	case bitFlip || phaseFlip:
		outcome = LogicalFlip
	case ambiguous:
		outcome = PossibleLogicalFlip
	default:
		overweight := c.overweightAtLevels(e)
		// Syndrome-free faults that are heavier than the guaranteed
		// correction radius, or faults whose minimum-weight explanation
		// succeeded while the injected pattern exceeded the radius at
		// both code levels, cannot be certified harmless: report them
		// honestly instead of claiming success.
		if c.syndromeFreeOverweight(syn, e) || overweight {
			outcome = PossibleLogicalFlip
		}
	}
	return &DecodeResult{
		Correction:        corr,
		Outcome:           outcome,
		Residual:          residual,
		BitLogicalFlip:    bitFlip,
		PhaseLogicalFlip:  phaseFlip,
		AmbiguousSyndrome: ambiguous,
	}, nil
}

// overweightAtLevels reports whether the injected fault exceeds the
// guaranteed correction radius at either level: more than t X faults in any
// one block (bit level) or more than t blocks with odd Z parity (phase
// level). Such faults may be corrected by luck (their syndrome can coincide
// with a light error's) but success cannot be guaranteed, so the result is
// reported as a possible logical flip.
func (c *Code) overweightAtLevels(e *ErrorSet) bool {
	t := c.CorrectableWeight()
	for b := 0; b < c.blocks; b++ {
		w := 0
		for i := 0; i < c.size; i++ {
			if e.X[c.Qubit(b, i)] {
				w++
			}
		}
		if w > t {
			return true
		}
	}
	badBlocks := 0
	for b := 0; b < c.blocks; b++ {
		var parity bool
		for i := 0; i < c.size; i++ {
			parity = parity != e.Z[c.Qubit(b, i)]
		}
		if parity {
			badBlocks++
		}
	}
	return badBlocks > t
}

// syndromeFreeOverweight reports a zero-syndrome fault whose physical weight
// exceeds the per-level correction guarantee, so the decoder cannot rule out a
// minimum-weight logical action.
func (c *Code) syndromeFreeOverweight(syn *Syndrome, e *ErrorSet) bool {
	zero := true
	for b := range syn.BitSyndrome {
		for _, s := range syn.BitSyndrome[b] {
			if s {
				zero = false
			}
		}
		for _, s := range syn.PhaseInner[b] {
			if s {
				zero = false
			}
		}
	}
	for _, s := range syn.PhaseOuter {
		if s {
			zero = false
		}
	}
	if !zero {
		return false
	}
	// Any syndrome-free fault above the guarantee can be logical. At the
	// bit level the guarantee is per block; at the phase level per block
	// parity. Check both raw levels.
	t := c.CorrectableWeight()
	for b := 0; b < c.blocks; b++ {
		w := 0
		for i := 0; i < c.size; i++ {
			if e.X[c.Qubit(b, i)] {
				w++
			}
		}
		if w > t {
			return true
		}
	}
	badBlocks := 0
	for b := 0; b < c.blocks; b++ {
		var parity bool
		for i := 0; i < c.size; i++ {
			parity = parity != e.Z[c.Qubit(b, i)]
		}
		if parity {
			badBlocks++
		}
	}
	return badBlocks > t
}

// hasBitLogical reports whether the residual X error contains the block-level
// logical operator X_L. Within a block the Z-check stabilisers are adjacent
// pairs, so any residual with zero edge-syndrome is either even-weight
// (stabiliser-equivalent to identity) or odd-weight; the latter is X_L.
func (c *Code) hasBitLogical(r *ErrorSet) bool {
	for b := 0; b < c.blocks; b++ {
		first := r.X[c.Qubit(b, 0)]
		uniform := true
		for i := 0; i < c.size; i++ {
			if r.X[c.Qubit(b, i)] != first {
				uniform = false
			}
		}
		// Post-decoding residual in a block has zero edge syndrome, so it
		// is uniform; the odd-weight uniform residual is X_L.
		if uniform && first {
			return true
		}
	}
	return false
}

// hasPhaseLogical reports whether the residual Z error contains Z_L. The
// phase-level code sees only per-block parity, so a block with odd residual
// parity and zero inter-block syndrome is a logical Z action.
func (c *Code) hasPhaseLogical(r *ErrorSet) bool {
	for b := 0; b < c.blocks; b++ {
		first := r.Z[c.Qubit(b, 0)]
		uniform := true
		var parity bool
		for i := 0; i < c.size; i++ {
			if r.Z[c.Qubit(b, i)] != first {
				uniform = false
			}
			parity = parity != r.Z[c.Qubit(b, i)]
		}
		// Uniform odd-weight residual in a block is Z_L; non-uniform
		// residuals have inner edge syndrome and are returned as trivial
		// or as PossibleLogicalFlip via the overweight zero-syndrome rule.
		if uniform && parity {
			return true
		}
	}
	return false
}

// decodeRepetition solves the minimum-weight perfect matching on a line for
// the length-n repetition code, given its n-1 edge syndromes.
//
// Flipping site p toggles edges p-1 and p (edge i joins sites i and i+1);
// pairing defect edges a<b therefore flips sites a..b-1. Odd numbers of
// defects pair the leftover with one of the two virtual boundary defects
// (left of edge 0, cost a+1; right of edge n-1, cost n-1-a).
//
// On a line the minimum matching is found by the interval DP
//
//	dp[i] = min cost to match the first i defect points (i even),
//
// where point i-1 is paired with some point j, with the cost of the path
// between them. The same DP counted with multiplicities reports ties.
func decodeRepetition(edges []bool) ([]bool, bool, error) {
	n := len(edges) + 1
	defects := make([]int, 0, len(edges))
	for i, s := range edges {
		if s {
			defects = append(defects, i+1)
		}
	}
	k := len(defects)
	// Points to match: the real defects plus, when k is odd, the cheaper
	// virtual boundary. Both choices are evaluated by the DP below.
	// match returns the minimum total path length, a list of matched
	// point-index pairs, and the number of distinct minimum-cost matchings.
	// Distances on the line are d(a,b) = pts[b]-pts[a] for a<b.
	pairAdjacent := func(defects []int) [][2]int {
		var pairs [][2]int
		for i := 0; i+1 < len(defects); i += 2 {
			pairs = append(pairs, [2]int{defects[i], defects[i+1]})
		}
		return pairs
	}
	pathCost := func(pairs [][2]int) int {
		total := 0
		for _, pr := range pairs {
			a, b := pr[0], pr[1]
			switch {
			case a == -1:
				total += b + 1
			case b == n:
				total += n - 1 - a
			default:
				total += b - a
			}
		}
		return total
	}

	var chosen [][2]int
	ambiguous := false
	switch {
	case k%2 == 0:
		// No boundary: pairing adjacent defect edges is uniquely minimum
		// on a line, and no ties exist between distinct defect pairs.
		chosen = pairAdjacent(defects)
	default:
		// One leftover defect pairs with the nearer virtual boundary; ties
		// give two distinct minimum-weight corrections.
		a := defects[k-1]
		leftPairs := append([][2]int{{-1, a}}, pairAdjacent(defects[:k-1])...)
		rightPairs := append(pairAdjacent(defects[:k-1]), [2]int{a, n})
		costL, costR := pathCost(leftPairs), pathCost(rightPairs)
		switch {
		case costL < costR:
			chosen = leftPairs
		case costR < costL:
			chosen = rightPairs
		default:
			chosen = leftPairs
			ambiguous = true
		}
	}

	chain := make([]bool, n)
	for _, pr := range chosen {
		a, b := pr[0], pr[1]
		// Virtual boundary -1 with edge defect b: flip sites 0..b.
		if a == -1 {
			for p := 0; p <= b; p++ {
				chain[p] = !chain[p]
			}
			continue
		}
		// Virtual boundary n with edge defect a: flip sites a+1..n-1.
		if b == n {
			for p := a + 1; p < n; p++ {
				chain[p] = !chain[p]
			}
			continue
		}
		for p := a; p < b; p++ {
			chain[p] = !chain[p]
		}
	}
	return chain, ambiguous, nil
}
