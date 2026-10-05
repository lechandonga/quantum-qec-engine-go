package circuit

import "errors"

// 可区分的输入校验错误。调用方可用 errors.Is 精确判断失败原因。
var (
	ErrInvalidQubit      = errors.New("circuit: qubit index out of range")
	ErrUnknownGate       = errors.New("circuit: unknown gate name")
	ErrControlConflict   = errors.New("circuit: control qubits must be distinct and exclude the target")
	ErrInvalidParam      = errors.New("circuit: invalid gate parameter")
	ErrInvalidState      = errors.New("circuit: invalid state")
	ErrInvalidObservable = errors.New("circuit: invalid pauli observable")
	ErrNoRNG             = errors.New("circuit: deterministic random source is required")
)
