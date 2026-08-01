# Test report — BioBirdDNA Genetics Skill v0.2.0

**Ngày kiểm tra:** 01/08/2026  
**Runtime chính:** Python 3.14.3; package yêu cầu Python >=3.11.  
**Runtime Node:** Node.js v24.13.1, npm 11.12.1; package yêu cầu Node.js >=20.  
**Runtime dependencies:** không.

## Kết quả phát hành

```text
110/110 Python tests passed
60/60 Node/npm tests passed
554/554 cloud, updater, parity, OAuth and plugin tests passed
7/7 golden workflows matched
51 probability distributions checked
same-input JSON is byte-stable
exact Python/Node simulation-recovery parity passed
npm packed-install smoke passed for both CLIs
npm global install and both CLI forward tests passed
Codex skill validator passed
123-file release tree is clean
SHA256SUMS.txt: 121 entries
PACKAGE_CHECKSUMS.json: 122 entries
```

`npm pack --dry-run --json` xác nhận đúng 8 file phát hành. Package là
`UNLICENSED`, không có runtime dependency và có metadata scoped public. Đây là
release candidate local, không phải bằng chứng npm registry đã publish.

## RED → GREEN evidence

- GEN-010 bắt đầu bằng test thất bại cho reciprocal direction, missing khác
  observed zero, stage loss, sex provenance và reconciliation đầy đủ.
- GEN-011 bắt đầu bằng test thất bại cho deterministic recovery AR/ZR,
  small-sample uncertainty, calibration và absolute-fit misspecification.
- Parity test ban đầu phát hiện seed/order khác nhau; Python và Node nay dùng cùng
  SHA-256 seed 53-bit, LCG, UTF-8 bucket order và bootstrap seed material.
- Contract fuzz phát hiện coercion, field lạ, label collision, count quá lớn và
  numerical edge cases; mọi trường hợp được khóa bằng regression test.
- Review cuối phát hiện self-pair có thể vượt gate khi ngưỡng F là 1 và feedback
  có thể nhận `patch`/`modify` hoặc claim “confirmed as causal”/“locus
  identified”. Test đỏ tái hiện từng bypass trước khi sửa; reviewer chạy lại và
  xác nhận PASS.
- Không xóa test, không nới evidence boundary để đạt GREEN.
- GEN-015 governance bắt đầu bằng test đỏ do thiếu deterministic release inventory
  builder và backlog issue; plugin packaging bắt đầu bằng test đỏ do thiếu
  configurator. Review-RED sau đó tái hiện rollback xóa nhầm target, PATH/env quá
  rộng, placeholder app ID, manifest drift và validator nằm ngoài transaction;
  toàn bộ đều có regression test trước khi sửa.
- Q0/Q1 hard-abstention bắt đầu bằng RED trong Python, npm và parity: feature cực
  mạnh vẫn trả `MEL-BRN`/`MEL-ALB` vì `UNK` chỉ được cộng điểm mềm. GREEN chuyển
  cổng chất lượng thành hard gate trước scoring trên modular Python, standalone,
  npm/Node và vendored Edge; kết quả bắt buộc `UNK`, `needs_recapture:true`, giữ
  phenotype-only limits và không có gene/genotype inference.
- Security review plugin bắt đầu bằng RED cho file `.env`/payload sinh học ngoài
  allowlist, junction/reparse root, alias đường dẫn Windows chỉ khác hoa thường và
  mutation vào external target. GREEN khóa exact package closure, realpath
  containment, same-volume promotion và rollback ở cả ba cổng official validator;
  focused suite đạt 35/35 và hai review cuối không còn P0/P1/P2.
- GEN-016 bắt đầu bằng RED cho real packed tarball: package cài được nhưng còn
  `private:true` (`true !== false`). GREEN xóa cờ private, thêm public scoped
  repository metadata, pack/install vào consumer sạch, chạy cả hai CLI, so exact
  cross golden và giữ feedback `pending_human_review`, `auto_apply:false`,
  `target_files:[]`. Test không dùng token, npm login, publish hoặc registry read.

## Coverage chính

- Sáu public operation có Python/Node parity: `cross`, `infer`, `pigment`,
  `rank-mates`, `normalize-records`, `simulate-recovery`.
- GEN-010 biểu diễn pairing, clutch, mọi offspring/loss stage, phenotype timeline
  và nguồn xác định giới tính; reciprocal cross vẫn phân biệt được.
- GEN-011 phục hồi model AR/ZR trên synthetic design, hiển thị uncertainty, giữ
  tie deterministic và cảnh báo candidate set misspecified bằng absolute-fit
  bootstrap có refit.
- Inference luôn yêu cầu đúng một `M0_NONGENETIC` cùng ít nhất một model đối chứng;
  zero-likelihood configuration không làm mất toàn bộ model.
- Posterior chỉ là hỗ trợ tương đối trong candidate set, không phải “xác suất
  gene”.
- Ảnh chỉ dùng cho phenotype screen; Q0/Q1 luôn abstain thành `UNK` +
  `needs_recapture:true` dù feature mạnh, còn `gene_status` vẫn `not_assessed`.
- `same_individual_pair`, full sibling, direct ancestor, cận huyết và welfare là
  hard gate; điểm phenotype không thể làm cặp bị chặn trở thành eligible.
- JSON không hữu hạn, input cycle, field lạ, kiểu coercible, phenotype label
  trùng, trọng số âm/không hữu hạn và work budget quá giới hạn đều fail closed.
- `biobirddna-feedback` tạo artifact SHA-256 deterministic, immutable,
  `pending_human_review`, `auto_apply:false`, `target_files:[]`.
- Feedback không gọi network, child process, AutoResearch mutation/promotion hoặc
  proposal ledger; AutoResearch health gate 3.0 không được dùng làm điểm
  BioBirdDNA.
- Packed-install smoke chạy cả `biobirddna-genetics` và
  `biobirddna-feedback` ngoài source tree và trên đường dẫn có khoảng trắng.
- Installer skill kiểm checksum ở staging và sau atomic no-replace rename.
- OAuth MCP expose đúng bảy read-only tool, kiểm exact audience và allowed user
  trước dispatch, không có anonymous/issuer-only fallback hoặc generic dispatch.
- Edge/source golden parity, changed deploy, unchanged no-op và rollback đều có
  live evidence; dedicated Supabase project không có application table, bucket,
  queue hoặc model-run persistence.
- Private ChatGPT Web App đã refresh đúng bảy action và gọi phenotype tool thật;
  universal personal plugin/Codex cache khớp byte, fresh Codex session gọi cùng
  registered OAuth MCP và giữ `gene_status:"not_assessed"`.

## Giới hạn còn lại

- Chưa có property-based test trên dải seed rộng; exact simulation parity hiện
  được khóa bằng seeded regression fixture và implementation audit.
- Bộ lọc claim/action free-text của feedback là lexical; structured fields vẫn là
  đường an toàn hơn cho phiên bản sau.
- Chưa có genotype-specific survival/observation model (GEN-012), multi-locus,
  database/API integration hoặc external-dataset validation.
- GEN-060 vẫn pending: patch này chỉ khóa acceptance Q0/Q1 của rule engine hiện
  tại, chưa xây calibrated image feature adapter hoặc external-device validation.
- Nhánh installer Linux/macOS được review tĩnh nhưng không chạy trên native host
  trong phiên Windows này.
- Supabase OAuth 2.1 vẫn ở trạng thái beta; release fail closed nếu exact-resource
  binding hoặc refresh behavior thay đổi.
- Candidate này không publish package lên npm registry và không submit public
  plugin directory; cả hai cần ủy quyền riêng cùng evidence native installers,
  public authenticated ChatGPT review và secure updater renewal.
- Synthetic recovery và model-relative posterior không xác nhận gene, variant
  hay genotype sinh học của cá thể thật.
