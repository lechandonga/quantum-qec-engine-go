# quantum-qec-engine-go

完全本地运行的量子线路模拟与量子纠错引擎。只依赖 Go 标准库，
不接入任何外部服务，线路演化、噪声抽样与解码结论均可在本地复现验证。

## 包结构

- `circuit`：量子线路的构造、纯态概率幅演化、测量塌缩与 Pauli 期望值。
- `noise`：Pauli 噪声通道、密度矩阵精确演化、确定性并行批量抽样。
- `qec`：广义 Shor 稳定子码 `[[d², 1, d]]` 的综合征提取、多数判决解码、
  逻辑翻转标记与编码态级恢复验证。

## 线路模型约定

- 比特编号 `0 .. n-1`，计算基整数的第 `q` 位对应比特 `q`。
- 概率幅向量长度 `2^n`，支持的比特数上限为 63。
- 固定门：`I X Y Z H S T`；带参数旋转门（角度为弧度）：
  `RX(θ) RY(θ) RZ(θ) PHASE(θ)`。
- 任意门都可带若干控制位（全部为 `|1>` 时生效），控制位之间、
  控制位与目标位必须互不相同。
- 门在追加到线路前完成全部校验，非法输入以可 `errors.Is` 区分的错误拒绝：
  - `circuit.ErrInvalidQubit`：比特下标越界；
  - `circuit.ErrUnknownGate`：未知门名；
  - `circuit.ErrControlConflict`：控制位重复或与目标冲突；
  - `circuit.ErrInvalidParam`：旋转角为 `NaN`/`Inf`；
  - `circuit.ErrInvalidState` / `ErrInvalidObservable`：态或观测量维度非法；
  - `circuit.ErrNoRNG`：测量未提供确定性随机源。

## 测量与可复现性

- `State.Measure` 按 Born 概率抽样，并把与结果矛盾的分量清零、重新归一化，
  严格符合测量塌缩语义；对同一比特重复测量结果不变。
- 随机源必须由调用方传入 `*rand.Rand`（标准库 `math/rand/v2`）。
- 噪声批量抽样与 QEC 试验都用 splitmix64 从 `(主种子, 样本编号)`
  派生每个样本的独立子种子，子种子序列与协程数量、调度顺序无关：
  **同一主种子下，无论并发度是多少，逐样本结果与聚合计数完全一致。**

```go
c, _ := circuit.New(2)
_ = c.H(0)
_ = c.X(1, 0)               // CNOT，构造 Bell 态
s, _ := circuit.NewState(2)
_ = s.Run(c)
p := s.Probabilities()      // [0.5, 0, 0, 0.5]
```

## 噪声模型约定

| 通道 | 名称 | 非平凡 Kraus 算符 |
| --- | --- | --- |
| 比特翻转 | `noise.BitFlip` | `sqrt(p) X` |
| 相位翻转 | `noise.PhaseFlip` | `sqrt(p) Z` |
| 退极化 | `noise.Depolarizing` | `sqrt(p/3) X/Y/Z` |

- 噪声默认在无噪线路演化之后、测量之前，对每个比特独立施加；
  可用 `Model.Qubits` 指定比特子集。
- 强度合法区间：比特/相位翻转 `[0,1]`，退极化 `[0,4/3]`
  （完全混合态对应 `p=3/4`，`4/3` 为等概率 X/Y/Z 参数化的数学上限）。
- `p=0` 时精确分布、密度矩阵迹与抽样结果与无噪线路完全一致。
- 两种计算口径：
  - `Engine.ExactProbabilities` / `ExactWithMidCircuitNoise`：
    密度矩阵 Kraus 精确演化，返回计算基概率与迹，迹漂移超过容差会报错；
  - `Engine.Sample` / `Engine.Batch`：态向量 + Pauli 错误抽样，
    `Batch` 自动并行并保证种子级确定性。
- `noise.TotalVariation` 给出经验分布与精确分布的总变差距离，
  用于统计收敛校验。
- 资源上限 `noise.Limits{MaxQubits, MaxShots, MaxWorkers}`：
  超限任务在开始计算前返回 `noise.ErrResourceLimit`，不会崩溃或截断。

## 纠错与容错判定规则

- 码：广义 Shor 码，码距 `d`（奇数，`3 <= d <= 7`），
  `N=d²` 个物理比特，`d` 个相位重复块、每块 `d` 个比特翻转重复比特。
- 物理错误：对称 Pauli 噪声，每比特独立以概率 `p/3` 发生 X、Y、Z。
- 综合征：
  - 块内 Z 型稳定子给出相邻比特 X 奇偶，定位 X/Y 错误；
  - 块间 X 型稳定子给出相邻块全体 Z 奇偶，定位相位错误。
- 解码：两层重复码多数判决，保证纠正半径 `t=(d-1)/2`。
  - 子码错误数 `> t`（或平局边界）时 `Outcome.Correctable=false`；
  - 对两个互补恢复方向分别计算残余错误的稳定子陪集逻辑分量，
    只有存在使残余完全落入稳定子空间的方向才承诺成功；
  - 裸逻辑算符（平凡综合征的逻辑 X̄/Z̄）与高权重错误会显式置
    `LogicalX`/`LogicalZ`，消息为 `possible logical flip ...`，
    绝不静默返回成功。
- 逻辑判定用 GF(2) 高斯-若当消元把残余投影到稳定子陪集，
  再与逻辑算符代表做辛内积：
  - 逻辑 `X̄ = X^⊗N`（全局 X）；
  - 逻辑 `Z̄ =` 跨全部 `d` 个块、每块块首一个 Z 的权重 `d` 链。
- `Encoder` 提供编码电路与恢复电路，可在真实态向量上验证
  「编码 → 注入错误 → 综合征/解码 → 恢复 → 逻辑保真度」全流程。
- QEC 仿真同样有资源上限 `qec.Limits` 与确定性并行保证，
  超限返回 `qec.ErrResourceLimit` / `ErrInvalidCode` / `ErrInvalidErrorRate`。

## 容差设置

- 态归一化校验：`1e-9`。
- 密度矩阵迹漂移：超过 `1e-9` 报错。
- 精确解析解对比：默认 `1e-12`。
- 统计收敛（40000 样本）：总变差距离阈值 `0.01`。
- 编码态恢复保真度：测试中要求 `> 1 - 1e-9`。

## 本地验证

无需联网、无需第三方依赖：

```bash
go test ./...                      # 全部单测
go test -v ./...                   # 查看每个业务/用例与中文结果日志
CGO_ENABLED=1 go test -race ./...  # 竞态检测（并行抽样必须通过）
go vet ./...
gofmt -l .
```

固定随机种子的用例（测量、统计收敛、并行确定性、解码趋势）在任何机器、
任何并发度下都应稳定通过。
