# BioBirdDNA Node feedback → AutoResearch review contract

## Mục đích

`biobirddna-feedback` biến một normalized production-review envelope đã hoàn tất
thành artifact feedback bất biến để reviewer cân nhắc tạo proposal cải thiện.
Adapter không phải AutoResearch experiment runner và không được tự sửa code/config.

## Ranh giới

- Không gọi `record_experiment`, promotion hoặc rollback của TAICHEL AutoResearch.
- Không dùng health gate 3.0 của Page làm điểm BioBirdDNA.
- Không lấy run-lock, không tạo action/queue, không gửi network request.
- Không ghi vào TAICHEL engine hay proposal ledger.
- Không tạo, approve hoặc apply proposal.
- Output luôn có `proposal.status:"pending_human_review"`,
  `auto_apply:false`, `target_files:[]`.

Reviewer người có thể dùng artifact đã xác minh để viết một entry riêng vào
`docs/PROPOSALS-next-improvements.md`. Bước đó nằm ngoài adapter và không tự động.

## Input production envelope tối thiểu

```json
{
  "schema_version": "1.0.0",
  "run_id": "CM-RUN-...",
  "operation": "cross",
  "engine_version": "0.3.0",
  "input_sha256": "64-hex",
  "output_sha256": "64-hex",
  "package": {
    "name": "@biobirddna/genetics",
    "version": "0.3.0",
    "integrity_sha256": "64-hex"
  },
  "result": {
    "output": {
      "ok": true,
      "operation": "cross"
    },
    "evidence": {
      "evidence_grade": "G0"
    },
    "limits": {
      "decision_support_only": true,
      "gene_status": "not_assessed"
    },
    "missing_data": []
  },
  "human_review_status": "pending"
}
```

`result.output` có thể giữ raw output của một trong sáu public operation để audit.
Runner phải bọc output đó cùng `result.evidence` và normalized review limits; không
đưa raw dispatch result trực tiếp vào adapter như thể nó đã là review envelope.
`engine_version`, package version và `gene_status` phải khớp chính xác `0.3.0`,
`0.3.0` và `not_assessed`; whitespace không được tự trim để chữa input sai.

Adapter phải nhận expected package integrity qua CLI/config độc lập và so sánh
constant-time với envelope. Hash input/output là provenance do production runner
tạo; adapter không tự bịa lại khi thiếu.

## Output `biobirddna_feedback_v1`

```json
{
  "schema_version": "biobirddna_feedback_v1",
  "feedback_id": "sha256-stable-id",
  "source": {
    "package": "@biobirddna/genetics",
    "package_version": "0.3.0",
    "package_integrity_sha256": "64-hex"
  },
  "run": {
    "run_id": "CM-RUN-...",
    "operation": "cross",
    "input_sha256": "64-hex",
    "output_sha256": "64-hex",
    "engine_version": "0.3.0",
    "schema_version": "1.0.0"
  },
  "evidence": {
    "evidence_grade": "not_assessed",
    "missing_data": [],
    "decision_support_only": true,
    "human_review_status": "pending",
    "limits": {}
  },
  "signals": [],
  "proposal": {
    "status": "pending_human_review",
    "auto_apply": false,
    "target_files": []
  }
}
```

`feedback_id` là SHA-256 của canonical JSON gồm source provenance, run provenance,
evidence và signals; không gồm timestamp hoặc đường dẫn máy.

## Fail-closed

Từ chối:

- JSON không hợp lệ, root không phải object hoặc vượt size limit;
- operation ngoài public allowlist;
- thiếu ID/hash/version hoặc hash không phải 64 hex;
- package integrity không khớp expected value;
- `decision_support_only !== true`;
- normalized review `gene_status !== "not_assessed"`;
- `human_review_status !== "pending"`;
- yêu cầu mutation/apply/approve/patch/modify, `auto_apply:true` hoặc
  `target_files` khác rỗng;
- cố tái dùng `health_gate_score`/AutoResearch score làm score BioBirdDNA;
- ngôn ngữ xác nhận gene dưới G5;
- output không có limits/evidence.

Key guard được canonicalize qua snake/camel/Pascal case. Vì vậy các dạng như
`mutationIntent`, `geneClaim`, `variant`, `autoResearchScore` hoặc
`AutoResearch.HealthGate.score` không được dùng để lách gate. Metric nghiệp vụ
trung tính như phenotype `score`, mate-ranking `score` và `fraction` vẫn hợp lệ
khi không nằm trong ngữ cảnh AutoResearch health gate hay mutation/claim.

Nếu evidence thiếu nhưng envelope vẫn hợp lệ, chỉ phát signal
`INSUFFICIENT_DATA`; không tạo score và không đề nghị auto-change. Các metric
nghiệp vụ hợp lệ trong `result.output` như pigment `scores`, mate-ranking `score`
hoặc `fraction` chỉ là raw provenance: adapter không sao chép chúng thành feedback
score.

## Verification

- cùng input canonical tạo cùng `feedback_id`;
- thay một hash làm ID đổi;
- integrity mismatch và mutation request fail closed;
- test monkeypatch network/file mutation để chứng minh adapter chỉ đọc input và
  ghi output được caller chỉ định;
- normalized wrapper chứa real output của cả sáu operation được chấp nhận nhưng
  claim xác nhận gene dưới G5 và health-gate score vẫn bị chặn;
- npm tarball fresh-install smoke phải chạy cả genetics CLI và feedback CLI.
