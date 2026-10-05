// Command qec-demo exercises the three packages end to end, with no flags
// required: it builds a Bell circuit, samples it under depolarising noise,
// and reports stabiliser-code logical-failure statistics versus the exact
// values for several code distances.
package main

import (
	"fmt"

	"github.com/lechandonga/quantum-qec-engine-go/internal/noise"
	"github.com/lechandonga/quantum-qec-engine-go/internal/qcircuit"
	"github.com/lechandonga/quantum-qec-engine-go/internal/qec"
)

func main() {
	// 1. Circuit construction and noisy sampling.
	c, err := qcircuit.NewCircuit(2)
	must(err)
	must(c.H(0))
	must(c.CNOT(0, 1))

	res, err := noise.RunBatch(noise.BatchRequest{
		Circuit: c,
		Model:   noise.Model{Depolarising: 0.02},
		Shots:   20000,
		Seed:    20261006,
		Workers: 4,
	})
	must(err)
	fmt.Println("noisy Bell sampling (depolarising p=0.02, 20000 shots):")
	for outcome := uint64(0); outcome < 4; outcome++ {
		fmt.Printf("  outcome %02b frequency %.4f\n", outcome, res.EmpiricalProbability(outcome))
	}

	// 2. Stabiliser code logical-failure trend.
	rates := qec.Rates{PX: 0.001, PZ: 0.001}
	fmt.Println("\nlogical failure vs code distance (pX=pZ=0.001):")
	for _, d := range []int{3, 5, 7, 9} {
		code, err := qec.NewCode(d)
		must(err)
		rep, err := qec.Sample(code, rates, qec.SampleConfig{
			Shots: 40000, Seed: 20261006, Workers: 4,
		})
		must(err)
		exact, err := qec.ExactLogicalFailure(d, rates)
		must(err)
		fmt.Printf("  d=%d n=%d sampled=%.6f exact=%.3e\n",
			d, code.NumQubits(), rep.LogicalFailureRate, exact)
	}
}

func must(err error) {
	if err != nil {
		panic(err)
	}
}
