# quantum-qec-engine-go

一个完全本地运行的量子线路模拟与量子纠错引擎。仅依赖 Go 标准库，不修改
`go.mod` 的 module 声明与 Go 版本，不访问任何外部服务。

## 目录结构

| 包 | 职责 |
| --- | --- |
| `internal/qcircuit` | 态矢量线路模拟：门、演化、概率幅、测量塌缩、Pauli 期望值 |
| `internal/noise` | Pauli 噪声信道、密度矩阵精确参照实现、确定性并行批量抽样 |
| `internal/qec` | 广义 Shor 稳定子码、综合征提取、最小权重解码、逻辑错误统计 |
| `cmd/qec-demo` | 端到端示例：含噪 Bell 采样 + 码距逻辑错误率趋势 |

## 线路模型约定（internal/qcircuit）

- 比特序：比特 `q` 对应计算基索引中的位 `1<<q`，即 `|q1 q0>`
  （`q0` 为最低位）存放在索引 `(q1<<1)|q0`。
- 初始态为 `|0...0>`，`State.Amplitudes/Probabilities/ProbabilityOf`
  只读返回拷贝；`Measure/MeasureAll` 使用调用方提供的 `*rand.Rand`
  就地塌缩并归一化，符合投影测量语义。
- 支持的门：
  - 单比特：`I X Y Z H S SDG T TDG`，带参数旋转 `RX RY RZ`（参数为弧度）；
  - 受控：`CNOT/CX CY CZ CH CRX CRY CRZ`，并且 `X/Y/Z/H/RX/RY/RZ` 本身
    也允许任意数量的控制位；
  - 双比特：`SWAP`（允许受控形式）。
- Pauli 期望值 `PauliExpectation([]byte)`，按比特位置传 `I/X/Y/Z`。
- 错误均为可区分的哨兵错误，可用 `errors.Is` 判断：
  `ErrQubitOutOfRange`、`ErrUnknownGate`、`ErrControlConflict`（含控制位
  与目标位冲突/控制位重复/门元数错误）、`ErrInvalidParameter`（NaN/Inf）、
  `ErrInvalidState`、`ErrInvalidPauli`、`ErrResourceLimit`。
- 资源上限：`DefaultLimits{MaxQubits:24, MaxGates:200000}`，也可用
  `NewCircuitWithLimits` 自定义；超限在建线路或加门之前即返回错误，绝不静默
  截断，态矢量上限为 63 比特。

## 噪声模型约定（internal/noise）

`Model{BitFlip, PhaseFlip, Depolarising, Qubits}` 描述逐比特独立的 Pauli
信道，三个信道按顺序逐比特复合；`Qubits == nil` 表示作用于全部比特：

- 比特翻转：`(1-p)ρ + p XρX`
- 相位翻转：`(1-p)ρ + p ZρZ`
- 退极化：`(1-p)ρ + p/3 (XρX + YρY + ZρZ)`

规则：

- 所有强度必须在 `[0,1]`（NaN 也拒绝，返回 `ErrInvalidStrength`），
  比特越界返回 `ErrInvalidQubit`。
- 强度全为 0 时逐样本与无噪结果**逐位完全一致**（测试中用 `!=` 比较）。
- 小规模精确结果由 `Density`（5 比特及以下密度矩阵，Kraus 精确演化，
  返回 `ErrResourceLimit` 拒绝更大系统）给出，用于与抽样统计对照；每一步
  校验迹不漂移。
- `RunBatch` 并行抽样：每个 shot 使用 `SplitMix64(seed, shotIndex)` 派生出
  独立的随机流，因此相同 `(seed, shots)` 下，无论 `Workers` 为 1 还是
  并行、无论调度如何，逐 shot 结果完全一致；结果按 shot 索引写入预分配切片，
  通过 `-race` 检测。
- `BatchLimits.MaxShots`（默认 10,000,000）与参数合法性在任何分配/抽样
  之前检查，超限直接返回 `ErrResourceLimit`。

统计收敛性用固定种子的大批量抽样验证（例如比特翻转 p=0.2，2 万样本落在
理论值 ±0.02 内）。

## 纠错与容错判定（internal/qec）

实现广义 Shor（CSS）码：码距 `d` 为奇数且 `d>=3`，物理比特数 `n=d*d`，
排成 `d` 个大小为 `d` 的块；保证纠正半径 `t=(d-1)/2`。

- 稳定子检查（综合征）：
  - 比特级：每块内相邻 `Z_i Z_{i+1}`（X 错误定位）；
  - 相位内层：每块内相邻 `X_i X_{i+1}`（Z 错误在块内的定位）；
  - 相位外层：相邻块的整块全 X 算符（块间 Z 奇偶定位）。
- 解码器：两级一维最小权重匹配（带虚拟边界与简并计数的区间 DP），
  综合征→最小权重纠正→残余分析，整个过程确定、可复现。
- 判定结果：
  - `Corrected`：唯一最小权重纠正后无逻辑作用；
  - `LogicalFlip`：残余包含逻辑算符（`X_L` 为块内全 X，`Z_L` 为
    单物理比特 Z），逻辑态确定被翻转；
  - `PossibleLogicalFlip`：最小权重解简并，或注入模式在任一码层超过
    `t`（即使其综合征恰好对应一个轻错误，也不保证纠正成功）。
  绝不把不可保证成功的错误静默报告为成功。
- 综合征提取是理想模型（检查本身无故障）。

`Rates{PX,PZ}` 为逐比特独立的 X/Z 概率（两次独立抽样，同时命中即为 Y）。
`Sample` 与噪声批量抽样使用同样的按 shot 索引派生随机流方案，固定种子、
固定 shots 时串行/并行逐样本一致；`SampleReport` 分别统计纠正、可能
翻转、确定翻转。

`ExactLogicalFailure(d, rates)` 用乘积结构给出去模拟误差的逻辑失败概率：
每块内层长度 `d` 的多数判决在 `d` 个块上复合（比特级），加上块奇偶
概率 `q=(1-(1-2*PZ)^d)/2` 上的长度 `d` 多数判决（相位级），取并集。
在固定低错误率（默认文档/测试用 `pX=pZ=0.001`）下，精确逻辑错误率随
`d=3,5,7,9,...` 单调下降；测试同时以固定种子抽样验证 d=3→d=5 的
下降趋势及与精确值的容差一致。

## 可复现性与容差

- 所有随机性都来自显式传入的 `*math/rand.Rand`；库内部不使用全局随机源。
- 固定种子（测试种子：线路 `42`、噪声 `7/11/12345`、纠错 `20261006`）。
- 默认容差：精确密度矩阵对照 `1e-12`、迹漂移 `1e-12`、态归一化
  `1e-9`、抽样统计 `±0.02`（2 万样本量级）、抽样对精确值使用
  `5*sqrt(p(1-p)/N)` 的蒙特卡洛容差。

## 本地验证方法

```bash
go build ./...
go vet ./...
go test ./... -v          # 子测试名称展示每个业务用例与通过情况
CGO_ENABLED=1 go test -race ./...   # 竞态检测（需要 cgo）
go run ./cmd/qec-demo     # 端到端演示
```

测试按业务组织，`-v` 日志可直接看到覆盖的非法线路、测量与采样语义、
噪声退化/精确对照/统计收敛/并行确定性，以及可纠/不可纠错误判定和码距
趋势；所有随机用例固定种子，重复运行结果稳定。
