# UMI 分子共识去噪 API

基于 **Python 3.11 + FastAPI + SQLite** 的 UMI（Unique Molecular Identifier）分子共识
去噪服务。对同一 UMI 分子的多条测序 reads 按 **Phred 质量加权投票**生成一致序列，
剔除 PCR/测序错误；对无法可靠判定的簇（多代表命中、低覆盖、低质量优势不足）整簇淘汰。

- `POST /consensus`：提交一批 reads，创建一个**不可变**的版本化结果（同 `request_id`
  同载荷幂等；载荷冲突返回 `409`；同一 `read_id` 在同请求或跨请求携带不同内容均
  返回 `409`）。
- `GET /consensus/{version}`：按版本号读取创建时的**不可变快照**。
- `GET /healthz`：健康检查。

---

## 1. 算法说明

### 1.1 UMI 分簇（在每个 `locus` 内独立进行）

1. 统计每个不同 UMI 的频次，按 **频次降序、字典序升序** 排序。
2. 依序选择簇代表（representative）：当前 UMI 与所有已选代表的汉明距离都 **> 1**
   时，它成为新代表；否则视为非代表。
3. 每个非代表 UMI 归入与其汉明距离 **≤ 1** 的**所有**代表：
   - 只命中 1 个代表 → 加入该簇；
   - 命中 ≥ 2 个代表 → 合并这些簇，并将**整个合并簇标记为 `ambiguous`**，
     该簇直接淘汰，不进行投票（原因码 `ambiguous_multi_representative`）。
4. UMI 长度不同视为不匹配（汉明距离无穷大）。

### 1.2 Phred 加权逐位投票

对簇内每个位置（0-based），对碱基 `b ∈ {A,C,G,T}` 累加权重：

```
W(b) = Σ 10^(Q_i / 10)        # 覆盖该位置的每条 read 在该碱基上的贡献
```

满足**全部**条件才给出碱基调用，否则该位置记为 `N`：

- 覆盖该位置的 read 数 `coverage ≥ 3`；
- 最高权重 `W_top` 严格占优：`W_top ≥ 2 × (其余碱基权重之和)`（并对浮点相等做
  平票保护：前两名权重相等时直接判 `N`）。

只要共识中出现任意 `N`，该簇整簇淘汰（原因码 `consensus_contains_n`，响应给出每个
N 位置的覆盖数、各碱基权重与是否平票）。

### 1.3 幂等与冲突

- 幂等键是 `request_id`。载荷指纹为规范化（read 行排序、去重后）JSON 的 SHA-256，
  因此 **read 顺序不影响**幂等性。
- 同一 `request_id`：
  - 相同指纹 → 返回已存版本，HTTP 头 `X-Idempotent-Replay: true`（不会产生新版本）；
  - 不同指纹 → **409 Conflict**。
- 同一请求内 `read_id` 相同但内容（locus/umi/sequence/qualities）不同 → **409**。
- **跨请求的 `read_id` 是全局的**：服务在 `read_registry` 表中持久记录每个已提交
  `read_id` 的规范内容（locus/umi/sequence/qualities 的 SHA-256 指纹，随首次接受的
  版本在同一条 `BEGIN IMMEDIATE` 事务中写入）。
  - 用**另一个 `request_id`** 提交同名 read，但 locus / umi / sequence / qualities
    **任一不同** → **409 Conflict**；不会分配新版本，也不会改动任何旧快照或已记录的
    read 内容。
  - 不同 `request_id` 提交**内容完全相同**的同名 read → 允许，正常创建新版本。
  - 并发的冲突提交由单条写事务串行化：**只有一种内容能落库**，其余请求收到 409。
- 完全相同的重复行允许提交，自动去重并计入 `stats.duplicate_read_rows`。

#### 旧版数据库升级

升级到本版本时，若库中还没有 `read_registry` 表，服务会自动建表并扫描已有版本快照：

- **旧版本照常可读**，`GET /consensus/{version}` 返回与升级前完全一致的快照；
- 旧快照的 `clusters[].read_ids` 只保留了 read_id，**无法还原原始的
  locus/umi/sequence/qualities**。这些 id 会以内容为 NULL 的「墓碑」记录登记：
  后续**任何复用该 read_id 的提交都返回 409**（即使新内容看起来合理），服务绝不
  伪造内容记录；
- 使用全新 read_id 的提交不受影响，正常创建版本。

### 1.4 不可变版本与读取一致性

- 版本号由 SQLite `AUTOINCREMENT` 单调分配，结果落盘后永不更新、删除；
- `GET /consensus/{version}` 永远返回创建那一刻的快照；
- 写入在单条 `BEGIN IMMEDIATE` 事务中完成（版本分配 + `requests`/`versions`/
  `read_registry` 三表写入与跨请求 read 内容比对全部原子化），冲突时整体回滚：
  不留空版本号、不登记任何 read；SQLite 配置 WAL 与 30s busy timeout，并发提交
  同样安全。

### 1.5 请求校验（422）

- `sequence` / `umi` 只允许 `A/C/G/T`；
- `qualities` 为 0～40 的整数，且与 `sequence` **等长**；
- `qualities` 也接受字符串：数字/逗号分隔形式（如 `"30, 20, 10"`）按整数字面量解析，
  其余可打印字符串按 **FASTQ Phred+33 ASCII** 解析（如 `"????"` = `[30,30,30,30]`）。

---

## 2. 快速开始

### 方式 A：Docker Compose（推荐，一键启动）

```bash
docker compose up --build
```

启动后：

- 服务地址：http://localhost:8000
- Swagger UI：http://localhost:8000/docs ；OpenAPI：http://localhost:8000/openapi.json
- 健康检查：http://localhost:8000/healthz
- SQLite 数据保存在命名卷 `umi-data`（容器内 `/data/umi.db`）；升级镜像后旧库
  会自动迁移（见 [1.3 旧版数据库升级](#13-幂等与冲突)），旧版本可继续读取。

快速验证跨请求 `read_id` 冲突防护：

```bash
curl -sS -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' --data @examples/request.json          # 200 v1
curl -sS -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' --data @examples/cross_request_reuse.json    # 200 v2
curl -i  -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' --data @examples/cross_request_conflict.json # 409
```

停止：

```bash
docker compose down          # 保留数据
docker compose down -v       # 同时删除数据卷
```

### 方式 B：本地直接运行

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
UMI_DB_PATH=./data/umi.db .venv/bin/uvicorn app.main:app --reload --port 8000
```

- 交互式 API 文档（Swagger UI）：http://localhost:8000/docs
- OpenAPI：http://localhost:8000/openapi.json
- 健康检查：http://localhost:8000/healthz

运行测试：

```bash
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m pytest -q
```

---

## 3. 实际请求 / 响应示例

示例文件见 [`examples/request.json`](examples/request.json)，包含 4 个 locus、
5 个簇，分别演示：正常共识、权重不占优导致 N、多代表命中 ambiguous、覆盖不足淘汰。

### 请求

```bash
curl -sS -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/request.json
```

请求体：

```json
{
  "request_id": "demo-run-001",
  "sample_id": "sample-2026-1001",
  "reads": [
    {"locus": "IGH-A", "read_id": "a01", "umi": "AAAAA", "sequence": "ACGT", "qualities": [40, 40, 40, 40]},
    {"locus": "IGH-A", "read_id": "a02", "umi": "AAAAA", "sequence": "ACGT", "qualities": [40, 39, 40, 38]},
    {"locus": "IGH-A", "read_id": "a03", "umi": "AAAAA", "sequence": "ACGT", "qualities": [38, 40, 39, 40]},
    {"locus": "IGH-A", "read_id": "a04", "umi": "AAAAT", "sequence": "ACGT", "qualities": [35, 36, 35, 34]},

    {"locus": "IGH-A", "read_id": "a05", "umi": "CCCCC", "sequence": "TTTT", "qualities": [40, 40, 40, 40]},
    {"locus": "IGH-A", "read_id": "a06", "umi": "CCCCC", "sequence": "TATT", "qualities": [40, 40, 40, 40]},
    {"locus": "IGH-A", "read_id": "a07", "umi": "CCCCC", "sequence": "TATT", "qualities": [30, 30, 30, 30]},

    {"locus": "IGH-B", "read_id": "b01", "umi": "GGGGG", "sequence": "ACGT", "qualities": [30, 30, 30, 30]},
    {"locus": "IGH-B", "read_id": "b02", "umi": "GGGGG", "sequence": "ACGT", "qualities": [30, 30, 30, 30]},
    {"locus": "IGH-B", "read_id": "b03", "umi": "GGGGT", "sequence": "ACGT", "qualities": [30, 30, 30, 30]},
    {"locus": "IGH-B", "read_id": "b04", "umi": "GGGAT", "sequence": "ACGT", "qualities": [30, 30, 30, 30]},

    {"locus": "IGH-C", "read_id": "c01", "umi": "TTTTA", "sequence": "GATT", "qualities": [40, 40, 40, 40]},
    {"locus": "IGH-C", "read_id": "c02", "umi": "TTTTA", "sequence": "GATT", "qualities": [40, 40, 40, 40]},

    {"locus": "IGH-D", "read_id": "d01", "umi": "AGAGT", "sequence": "CGTA", "qualities": [37, 37, 37, 37]},
    {"locus": "IGH-D", "read_id": "d02", "umi": "AGAGT", "sequence": "CGTA", "qualities": [36, 37, 35, 37]},
    {"locus": "IGH-D", "read_id": "d03", "umi": "AGAGC", "sequence": "CGTA", "qualities": [37, 36, 37, 37]}
  ]
}
```

### 200 响应（真实运行输出，`v1`）

```
HTTP/1.1 200 OK
content-type: application/json
x-idempotent-replay: false
x-version: 1
```

```json
{
  "version": 1,
  "request_id": "demo-run-001",
  "sample_id": "sample-2026-1001",
  "consensus": [
    {"locus": "IGH-A", "representative_umi": "AAAAA", "sequence": "ACGT", "depth": 4, "ambiguous": false},
    {"locus": "IGH-D", "representative_umi": "AGAGT", "sequence": "CGTA", "depth": 3, "ambiguous": false}
  ],
  "clusters": [
    {"locus": "IGH-A", "representative_umi": "AAAAA", "members": ["AAAAA", "AAAAT"], "ambiguous": false, "depth": 4, "read_ids": ["a01", "a02", "a03", "a04"]},
    {"locus": "IGH-A", "representative_umi": "CCCCC", "members": ["CCCCC"], "ambiguous": false, "depth": 3, "read_ids": ["a05", "a06", "a07"]},
    {"locus": "IGH-B", "representative_umi": "GGGGG", "members": ["GGGAT", "GGGGG", "GGGGT"], "ambiguous": true, "depth": 4, "read_ids": ["b01", "b02", "b03", "b04"]},
    {"locus": "IGH-C", "representative_umi": "TTTTA", "members": ["TTTTA"], "ambiguous": false, "depth": 2, "read_ids": ["c01", "c02"]},
    {"locus": "IGH-D", "representative_umi": "AGAGT", "members": ["AGAGC", "AGAGT"], "ambiguous": false, "depth": 3, "read_ids": ["d01", "d02", "d03"]}
  ],
  "eliminated": [
    {
      "locus": "IGH-A",
      "representative_umi": "CCCCC",
      "reason": "1 position(s) failed coverage or Phred-weight dominance and were called N",
      "reason_code": "consensus_contains_n",
      "ambiguous": false,
      "depth": 3,
      "n_positions": [1],
      "details": {
        "positions": [
          {
            "position": 1,
            "coverage": 3,
            "tie": false,
            "weights": {"A": 11000.0, "C": 0.0, "G": 0.0, "T": 10000.0}
          }
        ]
      }
    },
    {
      "locus": "IGH-B",
      "representative_umi": "GGGGG",
      "reason": "UMI hit more than one cluster representative (Hamming distance <= 1); whole cluster is ambiguous",
      "reason_code": "ambiguous_multi_representative",
      "ambiguous": true,
      "depth": 4,
      "n_positions": [],
      "details": {"representatives": ["GGGAT", "GGGGG"]}
    },
    {
      "locus": "IGH-C",
      "representative_umi": "TTTTA",
      "reason": "4 position(s) failed coverage or Phred-weight dominance and were called N",
      "reason_code": "consensus_contains_n",
      "ambiguous": false,
      "depth": 2,
      "n_positions": [0, 1, 2, 3],
      "details": {
        "positions": [
          {"position": 0, "coverage": 2, "tie": false, "weights": {"A": 0.0, "C": 0.0, "G": 20000.0, "T": 0.0}},
          {"position": 1, "coverage": 2, "tie": false, "weights": {"A": 20000.0, "C": 0.0, "G": 0.0, "T": 0.0}},
          {"position": 2, "coverage": 2, "tie": false, "weights": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 20000.0}},
          {"position": 3, "coverage": 2, "tie": false, "weights": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 20000.0}}
        ]
      }
    }
  ],
  "stats": {
    "loci": 4,
    "raw_reads": 16,
    "distinct_reads": 16,
    "duplicate_read_rows": 0,
    "clusters": 5,
    "consensus": 2,
    "eliminated": 3,
    "ambiguous_clusters": 1
  }
}
```

关键解读：

- `IGH-A / CCCCC`：第 2 位（0-based pos=1）A 权重 11000（两条 `TATT`，Q40+Q30 =
  10000+1000）、T 权重 10000（`TTTT`，Q40），`11000 < 2 × 10000`，优势不足 → `N` → 淘汰；
- `IGH-B`：`GGGGT` 与代表 `GGGGG`、`GGGAT` 的汉明距离均为 1（两位点间距为 2
  因而都是代表），整簇合并并标记 `ambiguous` → 淘汰；
- `IGH-C`：深度仅 2，每个位置覆盖都 `< 3` → 全部 `N` → 淘汰。

### 幂等重放（200，不产生新版本）

```bash
curl -i -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/request.json
```

```
HTTP/1.1 200 OK
x-idempotent-replay: true
x-version: 1
```

响应体与首次完全一致（read 顺序打乱也视为同一载荷）。

### 同 request_id 不同载荷 → 409

```bash
# examples/conflict.json：仅把 a01 的第一个质量值从 40 改成 25，request_id 不变
curl -sS -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' --data @examples/conflict.json
```

```json
{"detail": {"error": "request_id='demo-run-001' already submitted with a different payload", "request_id": "demo-run-001"}}
```

### 同 read_id 内容冲突 → 409

同一请求内：

```json
{"detail": {"error": "read_id content conflict", "conflicts": [{"read_id": "dup", "reason": "same read_id with different content"}]}}
```

### 跨请求 read_id 内容冲突 → 409

`read_id` 全局唯一地标识一条 read 的内容。先用 `examples/request.json`
（`request_id=demo-run-001`，其中 `a02.sequence = ACGT`）提交后，再换一个
`request_id` 提交**内容相同**的 read 是允许的（产生新版本）：

```bash
# examples/cross_request_reuse.json：request_id=demo-run-002，
# 复用 a01..a04 且内容与首次提交完全一致
curl -i -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/cross_request_reuse.json
```

```
HTTP/1.1 200 OK
x-idempotent-replay: false
x-version: 2
```

但如果另一个请求里同名 read 的 locus / umi / sequence / qualities 任一不同
（下例把 `a02` 的 `sequence` 从 `ACGT` 改成 `ACGG`）：

```bash
# examples/cross_request_conflict.json：request_id=demo-run-003，
# read_id=a02 与 demo-run-001 中的内容不一致
curl -i -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/cross_request_conflict.json
```

```
HTTP/1.1 409 Conflict
```

```json
{"detail": {"error": "read_id='a02' was already submitted with different content (locus/umi/sequence/qualities)"}}
```

该请求整体回滚：不产生新版本（`/healthz` 的 `latest_version` 不变）、旧快照
不变，请求里其余全新的 read_id 也不会被登记。把 `request_id` 换成与已存内容
一致的载荷后可重新提交成功。

> 升级旧库后，若某 `read_id` 只在旧版本快照中出现过、原始内容已无法还原，
> 复用它会得到：
> `{"detail": {"error": "read_id='...' appeared in a pre-upgrade result whose original content cannot be recovered; the read_id cannot be reused"}}`。

### 读取不可变版本

```bash
curl -sS http://localhost:8000/consensus/1 | jq -r '.consensus[].sequence'
# ACGT
# CGTA

curl -i http://localhost:8000/consensus/999
# HTTP/1.1 404 Not Found
# {"detail":{"error":"version 999 not found"}}
```

---

## 4. 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 路由、启动、错误处理
│   ├── schemas.py       # Pydantic 请求/响应模型与校验
│   ├── consensus.py     # 分簇、加权投票、载荷/read 内容指纹（纯函数，可单测）
│   └── db.py            # SQLite 存储、版本分配、幂等与跨请求 read 冲突事务、旧库迁移
├── tests/               # pytest：算法单测 + API 端到端（39 个用例）
├── examples/request.json
├── examples/conflict.json                 # 同 request_id 异载荷
├── examples/cross_request_reuse.json      # 跨请求同内容复用 → 新版本
├── examples/cross_request_conflict.json   # 跨请求同名 read 异内容 → 409
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── requirements-dev.txt
```
