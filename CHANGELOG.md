# Changelog

## 0.3.0 — 2026-08-08 (release hardening)

- Cross-platform release: golden verification so sánh cây JSON với tolerance số
  (rel 1e-9 / abs 1e-12); byte-exact chỉ yêu cầu cho double-run cùng máy. Gate
  `validate_package.py` pass trên Windows / Linux (đã kiểm chứng cả hai).
- Node engine an toàn với key nguy hiểm: mọi accumulator/lookup nhận key từ
  input dùng null-prototype + Object.hasOwn; label `__proto__` / `constructor`
  / `prototype` cho kết quả đúng hệt Python (regression test kèm theo).
- Hàng rào evidence enforce bằng code: service áp claim-ceiling theo operation,
  mọi envelope mang `limits.claim_language` qua `safe_claim_language`.
- Kết quả infer thêm `warnings` (cả hai runtime): cảnh báo khi thứ hạng model
  có thể do prior cấu-hình quyết định thay vì dữ liệu.
- Contract `female_probability`: reject giá trị khác 0.5 khi có model Z-linked
  tham gia (engine + standalone, cả hai runtime).
- `alternative_hypotheses` của pigment chỉ chứa nhãn score > 0.
- Số học di động hai runtime: port fdlibm log/exp (+ log2, logGamma, logsumexp,
  entropy, quantile) bit-identical Python↔Node; tổng float dùng Neumaier trên
  dãy đã sort (bit-identical bất kể phiên bản Python — builtins.sum đổi thuật
  toán từ CPython 3.12) cho brier/means/entropy.
- Windows: subprocess/console ép UTF-8 (sửa lỗi cp1252 với message tiếng Việt),
  test installer .sh tự skip khi thiếu POSIX sh.
- Verifier/installer/manifest-builder miễn nhiễm cache dev-time
  (.pytest_cache/.git); vẫn chặn *.pyc lạ như một tính năng an ninh.
- CLI tách lỗi input (exit 2, slug + message_vi) khỏi lỗi nội bộ (exit 70,
  internal_error:<Tên>) — mirror ở cả Node bin. TEST_REPORT tách phần kiểm
  được từ chính gói khỏi bằng chứng hạ tầng ngoài gói.

## 0.2.0 — 2026-08-01 (public release candidate; not published)

- Hoàn thành GEN-016 release-candidate boundary: package npm có scoped public
  metadata, zero runtime dependency và packed-install smoke test cho cả hai CLI.
- Chuẩn bị repository metadata cho `vuho60-byte/biobirddna-genetics` nhưng không
  chứa credential, private operator state hoặc bằng chứng vận hành riêng tư.
- Không có `npm publish`, npm login, GitHub push hoặc xác nhận registry nào trong
  thay đổi này. Public authenticated ChatGPT review, native Ubuntu/macOS installer
  evidence và OAuth updater renewal vẫn là cổng external phải hoàn tất riêng.

## 0.2.0 — 2026-07-30

- Hoàn thành GEN-010 với contract breeding-record chuẩn hóa: giữ rõ reciprocal
  direction, phân biệt missing với observed zero, ghi từng giai đoạn thất thoát
  và đối soát stage count khi coverage đầy đủ.
- Hoàn thành GEN-011 với simulation-recovery deterministic cho AR/ZR, uncertainty,
  tie-safe recovery, calibration theo multi-model design và absolute-fit bootstrap
  có refit.
- Giữ đúng một `M0_NONGENETIC` và ít nhất một model đối chứng trong inference/
  recovery; posterior chỉ là hỗ trợ tương đối trong candidate set.
- Từ chối JSON `NaN`/`Infinity`, prior không hữu hạn, genotype khác `NA` ở
  generator non-genetic và semantic text có whitespace ngoài; trạng thái
  log-evidence bất khả thi dùng JSON `null`.
- Khóa label phenotype không được trùng, count inference tối đa `1_000_000`,
  hypothesis/pairing weights phải là số hữu hạn không âm và mọi pair score phải
  hữu hạn; input sai kiểu hoặc field lạ ở cross/rank fail closed thay vì coercion.
- Hard-block `same_individual_pair` bất kể ngưỡng cận huyết hoặc cờ quan hệ được
  cấu hình, nên cùng một cá thể không thể xuất hiện ở cả hai vai trò ghép cặp.
- Đồng nhất seed SHA-256 53-bit, LCG và thứ tự bucket cho simulation Python/Node;
  cùng request v0.2.0 tạo cùng output recovery ở cả hai runtime.
- Khóa semantic snapshot của các built-in model ID kể cả khi generating model bị
  loại khỏi candidate set; generator misspec tùy biến phải dùng ID riêng.
- Thêm package ESM `@biobirddna/genetics` cho Node.js >=20, zero dependencies,
  public API và CLI parity cho sáu operation; regression gate so output Node với
  Python trên golden và adversarial fixtures.
- Thêm `biobirddna-feedback` tạo artifact SHA-256 bất biến, proposal-only,
  fail-closed; không gọi, không sửa và không tự apply AutoResearch.
- Feedback guard nhận diện cả camel/Pascal/underscore mutation intent,
  gene/variant claim đã xác nhận và AutoResearch health-score trá hình, nhưng
  vẫn giữ các metric nghiệp vụ như phenotype score/fraction.
- Thêm schema, template, examples, ADR-005, tài liệu AutoResearch feedback và
  packed-install smoke test cho cả hai Node CLI.
- Phát hành dưới trạng thái `UNLICENSED`, chỉ tạo local tarball; không publish npm.
- Thêm GEN-015 authenticated stateless OAuth MCP adapter với đúng bảy read-only
  tool, exact Edge/Node parity, allowed-user và exact-audience gates; không thêm
  database, bucket, queue, model-run persistence hoặc generic dispatch.
- Thêm fail-closed Node updater/watcher local-to-Supabase, sensitive-change pin,
  immutable release snapshot, live verification và rollback.
- Thêm private ChatGPT Web MCP App và universal personal plugin packaging; plugin
  chỉ tham chiếu registered connection, bundle skill đã kiểm checksum và không
  chứa token/secret. Public directory submission vẫn ngoài phạm vi.
- Hardening acceptance Q0/Q1 của GEN-060: pigment Python modular, standalone và
  Node/Edge giờ hard-abstain thành `UNK` + `needs_recapture:true` trước scoring,
  nên feature mạnh không thể vượt cổng chất lượng ảnh; giới hạn phenotype-only và
  `gene_status:not_assessed` không đổi. GEN-060 vẫn chưa hoàn tất.

## 0.1.1 — 2026-07-30

- Hoàn thành GEN-001: khóa năm golden outputs, kiểm tra mọi phân bố xác suất và byte stability.
- Buộc JSON output dùng LF để ổn định byte trên Windows, Linux và macOS.
- Sửa Node.js standalone entrypoint cho đường dẫn Windows và đường dẫn có khoảng trắng.
- Thêm kiểm tra checksum/inventory/version và installer staging trước khi atomic rename.
- Chặn NaN/Infinity, distribution rỗng và manifest sai kiểu bằng lỗi có cấu trúc.
- Thêm installer lock, no-replace rename và kiểm checksum lại sau final rename.
- Shell installer dùng `renameat2(RENAME_NOREPLACE)` trên Linux hoặc `renamex_np(RENAME_EXCL)` trên macOS và fail-loud nếu không hỗ trợ.
- Ngăn CLI tạo bytecode trong cây skill và chuẩn hóa lệnh chạy với `python -B`.
- Giữ nguyên evidence boundary, `M0_NONGENETIC`, cổng cận huyết/phúc lợi và public JSON contract.

## 0.1.0 — 2026-07-29

- Tạo khung skill `biobirddna-genetics` theo cấu trúc `SKILL.md` + `agents/openai.yaml` + `references/` + `scripts/`.
- Thêm lõi Python standard-library cho năm mô hình một locus: lặn thường, trội thường, trội không hoàn toàn, lặn liên kết Z và trội liên kết Z.
- Thêm mô hình nền `non_genetic` để tránh ép mọi biến đổi màu vào một gene.
- Thêm mô phỏng kiểu gene cha mẹ chắc chắn hoặc phân bố xác suất kiểu gene chưa chắc chắn.
- Thêm bộ so sánh mô hình bằng likelihood/Bayes thăm dò.
- Thêm ma trận quan hệ phả hệ, hệ số cận huyết đời con và chặn cặp cha–con/anh chị em ruột.
- Thêm bộ xếp hạng cặp theo xác suất mục tiêu, giá trị thông tin, cận huyết, sức khỏe và bất định.
- Thêm bộ rule-based nhận định kiểu hình sắc tố, luôn giữ `gene_status:not_assessed`.
- Thêm JSON Schema, dữ liệu mẫu, unit test, CLI modular và bản code độc lập một file.
- Bổ sung bản Node.js một file cho cross simulation và tích hợp n8n/worker JavaScript.
- Bổ sung blueprint hợp nhất, handoff cho coding model, backlog có acceptance criteria và bốn ADR kiến trúc.
- Bổ sung draft PostgreSQL/Supabase schema và OpenAPI cho giai đoạn platform.
- Bổ sung manifest, file index, test report và output demo đã chạy.
