# p06_atavism — brief X8: xác suất một "gene cổ" hiện lại

Bối cảnh đầy đủ: `research/briefs/X8-xac-suat-hien-lai-gen-co.md`. Tham khảo phong cách CLI và cách ghi
TSV từ `pipeline/p05_admixture/dilution_sim.py` (không import module đó — mỗi thư mục trong
`pipeline/` tự chứa, đúng quy ước hiện có của dự án).

## Câu hỏi khoa học

Chủ dự án hỏi hai câu, và brief X8 đóng khung để cả hai trả lời được bằng số:

1. Xác suất một gene cổ **hiện lại** là bao nhiêu?
2. Có **bằng** xác suất xuất hiện một **đột biến mới** không?

Mô hình ở mức **phân tử**: 1 locus = danh sách phần tử đặt tên A, B, C, ... Z; phần tử **A** đã bị đột
biến làm hỏng chức năng, các phần tử còn lại còn nguyên ("gốc"). Câu trả lời phụ thuộc hoàn toàn vào
**k** — có bao nhiêu phần tử phải quay ngược — không phụ thuộc danh tính/thứ tự cụ thể của từng phần
tử (mô hình coi các phần tử `dot_bien` là đồng nhất/hoán đổi được, dùng chung 1 cặp `mu`/`kappa`/lựa
chọn đường — đơn giản hóa, xem "Giới hạn").

## Chạy — `reversion_prob.py`, bốn lệnh con

```
python -B pipeline/p06_atavism/reversion_prob.py locus --elements A,B,C --mutated A \
    --mu 4.6e-9 --kappa 2.0 --n-diploid 1e6 --generation-years 2 --out locus.txt

python -B pipeline/p06_atavism/reversion_prob.py compare --locus-length 1000 \
    --mu 4.6e-9 --kappa 2.0 --mutated-count 1 --out compare.tsv

python -B pipeline/p06_atavism/reversion_prob.py sweep --mu-list 1e-9,4.6e-9,1e-8 \
    --k-list 1,2,3,5 --n-list 1e4,1e6,1e8 --out sweep.tsv

python -B pipeline/p06_atavism/reversion_prob.py report --sweep sweep.tsv --out report.md
```

Chỉ Python 3.11 stdlib (`argparse`, `math`). Không numpy/scipy, không docker, không tải dữ liệu mạng,
không import module khác trong `pipeline/`. Mọi lệnh con có `--out` (hoặc bắt buộc, hoặc mặc định `-`
= stdout), theo đúng quy ước TSV/text streaming của `dilution_sim.py` (`smart_open_read`/
`smart_open_write`, `_fmt_cell` làm tròn float về 6 chữ số có nghĩa).

### Công thức gốc (từ brief X8, không tự suy diễn)

```
p_rev = mu * p_specific
p_specific = kappa/(kappa+2)   nếu dùng đường CHUYỂN ĐỔI (transition)
p_specific = 1/(kappa+2)       nếu dùng MỘT đường ĐẢO CHUYỂN (transversion) cụ thể
```

`mu` = tốc độ đột biến/vị trí/thế hệ; `kappa` = tỉ lệ tốc độ chuyển đổi/đảo chuyển. Mỗi base có 1 đường
chuyển đổi + 2 đường đảo chuyển, nên `p_specific(chuyển đổi) + 2×p_specific(đảo chuyển) = 1` (kiểm
trong test `test_transition_and_transversion_sum_to_one_for_various_kappa`).

**Mặc định dùng đường CHUYỂN ĐỔI** (`--transversion` để đổi sang đường còn lại) — đây là trường hợp
**thuận lợi nhất** (khi `kappa > 1`, `p_specific` lớn hơn ở đường chuyển đổi), nên con số đưa ra là
**giới hạn TRÊN** của khả năng hiện lại, đúng yêu cầu brief.

### Hai con đường quay ngược

- **`sequential_wait_*`** (TUẦN TỰ): k phần tử quay ngược lần lượt, mỗi bước trung tính (gene đã hỏng
  sẵn nên 1 phần tử quay ngược một mình không cho lợi thế chọn lọc nào) và phải **cố định** trong quần
  thể trước khi bước sau xảy ra. Thời gian chờ kỳ vọng mỗi bước ≈ `1/p_rev` thế hệ, **KHÔNG phụ thuộc
  N**. Tổng k bước ĐỘC LẬP → tổng kỳ vọng = tổng các kỳ vọng = `k/p_rev` (đúng với mọi phân phối, không
  cần giả thiết gì thêm ngoài "độc lập + kỳ vọng hữu hạn"). **Việc mỗi bước ≈ 1/p_rev the he VÀ độc lập
  N là GIẢ ĐỊNH — xem mục "Giả định" bên dưới, CHƯA được coi là CONFIRMED.**
- **`concurrent_wait_*`** (ĐỒNG THỜI): cả k phần tử cùng quay ngược trong 1 giao tử, 1 thế hệ. Xác suất
  mỗi giao tử = `p_rev**k`. Số ca kỳ vọng mỗi thế hệ trong CẢ quần thể = `2*N*p_rev**k` (2N bản sao
  nhiễm sắc thể trong quần thể lưỡng bội cỡ N) → thời gian chờ ≈ `1/(2*N*p_rev**k)`, **PHỤ THUỘC N**.

Đường **tuần tự cho thời gian NGẮN hơn đường đồng thời khi k≥2** (vì `p_rev<1` nên `p_rev**k` giảm
nhanh hơn hẳn `k` lần) — đây là kịch bản **thuận lợi nhất**, nên `sequential_wait_*` là **CON SỐ CHÍNH**
trả lời câu hỏi 1 của chủ dự án (kiểm trong `test_sequential_shorter_than_concurrent_when_k_at_least_two`;
lệnh `locus`/`report` đều ghi chú rõ điều này).

**Lưu ý khi k=1 (không phải lỗi):** hai công thức đo **hai việc khác nhau**. `concurrent` đo thời gian
để một đột biến mới **XUẤT HIỆN** ở đâu đó trong quần thể (`1/(2*N*p_rev)`, nhanh khi N lớn vì có nhiều
bản sao cùng "thử"), còn `sequential` (dùng công thức thay thế/cố định trung tính) đo thời gian để nó
**CỐ ĐỊNH VĨNH VIỄN** trong quần thể — hai mốc thời gian khác hẳn nhau về ý nghĩa (đây chính là điều
CẢNH BÁO BẮT BUỘC 1 nhấn mạnh: "xuất hiện" khác "lan ra quần thể"). Ví dụ `example/locus_A.txt`
(k=1, N=1e6): `concurrent_wait_generations`≈217 thế hệ (xuất hiện, rất nhanh) nhưng
`sequential_wait_generations`≈4,35×10⁸ thế hệ (cố định, rất chậm) — **KHÔNG mâu thuẫn** với phát biểu
"tuần tự ngắn hơn khi k≥2" ở trên, vì phát biểu đó chỉ áp dụng cho k≥2.

### So sánh bắt buộc — câu hỏi 2 của chủ dự án

Xác suất xuất hiện **một đột biến mới bất kỳ** ở đâu đó trong locus dài `L` base, 1 thế hệ, 1 bản sao
NST:

```
p_new = 1 - (1-mu)**L  ~=  L*mu
```

Lệnh `compare` in ti so `p_new / p_rev` (`p_rev` ở đây dùng công thức ĐỒNG THỜI với
`k = --mutated-count`, tức `(mu*p_specific)**mutated_count`) — trả lời THẲNG câu hỏi 2: đột biến mới dễ
xảy ra hơn quay ngược đúng 1 (hoặc k) phần tử cụ thể **bao nhiêu lần**.

### 1. `locus` — trạng thái từng phần tử + xác suất/thời gian chờ

In 2 bảng TSV cách nhau 1 dòng trống: `element/state` (trạng thái từng phần tử trong `--elements`, `goc`
hoặc `dot_bien` theo `--mutated`), rồi `metric/value/unit` (mu, kappa, đường dùng, `p_specific`, `p_rev`,
k, N, và xác suất/thời gian chờ **CẢ HAI** đường, tính **cả thế hệ VÀ năm** qua `--generation-years`).
`k` = số phần tử trong `--mutated` (không bắt buộc là tiền tố A,B,C... của `--elements`, xem "Giới hạn").

### 2. `compare` — trả lời câu hỏi 2

In 1 bảng TSV `metric/value/unit`: `p_new`, `p_rev`, `ratio_p_new_over_p_rev`.

### 3. `sweep` — lưới tham số → bền vững trước khoảng giá trị

Vì hằng số THẬT cho chào mào (`mu`, `kappa`) còn đang chờ brief R14
(`research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md`, **CHƯA xong** lúc viết script này — xem
`research/briefs/R14-lai-to-va-hang-so-toc-do-dot-bien.md`), lệnh này quét lưới `--mu-list × --k-list ×
--n-list` → `sweep.tsv` (cột: xem `SWEEP_HEADER` trong code), để kết quả không phụ thuộc 1 con số đinh.

### 4. `report` — bảng gọn + 3 cảnh báo bắt buộc

Đọc `sweep.tsv`, xuất `report.md` gồm 3 cảnh báo bắt buộc (nguyên văn dưới đây), 1 đoạn diễn giải
**tự động từ chính dữ liệu** (bao nhiêu dòng k≥2 có tuần tự ngắn hơn đồng thời — không hardcode), và
bảng đầy đủ theo tổ hợp (k, mu, N).

## Ổn định số học (yêu cầu riêng của brief X8)

- **`p_new`**: tính bằng `-math.expm1(L*math.log1p(-mu))` thay vì `1-(1-mu)**L` trực tiếp — tránh **hủy
  số** (catastrophic cancellation): khi `mu` rất nhỏ, `1-mu` có thể làm tròn về ĐÚNG `1.0` trong dấu
  phảy động, khiến `(1-mu)**L` cũng trả về `1.0` và `p_new` trả về ĐÚNG `0.0` dù đáp số thật là 1 số
  dương nhỏ (~`L*mu`). Kiểm bằng `mu=1e-18` (`test_p_new_does_not_catastrophically_cancel_for_tiny_mu`):
  `1.0-(1.0-1e-18)**1` trả về đúng `0.0`, nhưng `p_new_any_mutation` trả về đúng `1e-18`.
- **Xác suất đường đồng thời (`p_rev**k`)**: tính qua không gian LOG (`k*math.log(p_rev)` rồi mới
  `exp()`) thay vì lũy thừa trực tiếp — `p_rev**k` có thể tràn về ĐÚNG `0.0` khi k lớn/p_rev nhỏ
  (float64 không biểu diễn được số dương nhỏ hơn ~4,9×10⁻³²⁴), khiến bước tính `1/rate` tiếp theo ném
  `ZeroDivisionError`, hoặc `math.log(p_rev**k)` ném `ValueError` — dù đáp số thật (log-xác-suất, hoặc
  thời gian chờ) vẫn tính được đúng qua log. Kiểm bằng k=10 (yêu cầu rõ trong brief):
  `test_log_prob_avoids_domain_error_when_direct_power_underflows_at_k_ten` (p_rev tổng hợp = 1e-33,
  `p_rev**10` tràn về đúng `0.0`, `math.log(0.0)` ném lỗi thật, nhưng
  `log_concurrent_prob_per_gamete` không lỗi) và
  `test_concurrent_wait_finite_and_cross_checked_with_decimal_at_k_ten_realistic_params` (k=10 với
  `mu`/`kappa` mặc định THỰC TẾ của CLI — ở đây `p_rev**10`≈4,1×10⁻⁸⁷, KHÔNG tràn — đối chiếu kết quả
  với `decimal.Decimal` độ chính xác tùy ý, sai số tương đối < 10⁻⁹).
- Khi chính **rate** (`2*N*p_rev**k`) nhỏ tới mức bản thân nó cũng không còn biểu diễn nổi trong float64
  (dưới ~4,9×10⁻³²⁴), thời gian chờ toán học thật (nghịch đảo của rate đó) **vượt quá** giá trị lớn nhất
  float64 biểu diễn được (~1,8×10³⁰⁸) — lúc đó hàm trả về `math.inf` một cách **có chủ đích** (không
  `ZeroDivisionError`/`OverflowError`), vì đó là đáp số đúng, không phải giới hạn của cách tính. Kiểm
  bằng `test_concurrent_wait_returns_inf_not_crash_for_astronomically_tiny_rate`.

## Giả định

**Thời gian chờ mỗi bước của đường TUẦN TỰ ≈ 1/p_rev thế hệ, KHÔNG phụ thuộc N** — đây là kết quả CHUẨN
của lý thuyết trung tính (Kimura): số đột biến trung tính mới xuất hiện mỗi thế hệ trong quần thể =
`2N*p_rev`, xác suất cố định của MỖI đột biến trung tính = `1/(2N)`, nhân lại **N triệt tiêu đúng**, cho
tốc độ cố định = `p_rev`/thế hệ bất kể N. Đây là **GIẢ ĐỊNH — CHƯA TỰ MỞ NGUỒN XÁC MINH TRONG PHIÊN VIẾT
SCRIPT NÀY** (brief X8 cấm tải dữ liệu mạng ở bước này; kết quả là kiến thức giáo trình phổ biến nhưng
theo đúng luật bằng chứng của dự án — xem `tasks/lessons.md` mục "Trí nhớ của Main cũng KHÔNG phải bằng
chứng" — KHÔNG được gắn nhãn CONFIRMED/STRONG_INFERENCE khi chưa tự mở 1 nguồn thật trong chính phiên
làm việc này). Brief X8 yêu cầu rõ: "nếu bạn không xác minh được thì gắn nhãn giả định trong README,
đừng khẳng định" — đây chính là nhãn đó. `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md` (đang
chờ subagent khác) có thể mang về nguồn xác nhận cho kết quả này.

**Các giả định khác:**
- Mọi phần tử `dot_bien` dùng CHUNG 1 cặp `mu`/`kappa`/lựa chọn đường (đồng nhất/hoán đổi được) — đơn
  giản hóa; thực tế mỗi base cụ thể có thể nằm trên đường chuyển đổi hay đảo chuyển khác nhau tùy base
  gốc là gì.
- `--mutated` không bắt buộc là tiền tố A,B,C... của `--elements` như mô tả "để dễ hình dung" trong
  brief — CLI nhận trực tiếp tập hợp tên phần tử nào đang `dot_bien`; toán học chỉ phụ thuộc **số lượng**
  k (do giả định đồng nhất ở trên), không phụ thuộc thứ tự/danh tính, nên không cần ép prefix.
- `sequential_wait_generations`/`concurrent_wait_generations` dùng CHUNG 1 chữ ký gọi hàm (cả 2 đều
  nhận `n_diploid`, dù `sequential` không dùng tới) — chỉ để tiện vòng lặp `sweep` và để viết được phép
  thử phân biệt trực tiếp (xem `tests/`), không mang ý nghĩa sinh học.
- Công thức `p_new = 1-(1-mu)**L` giả định các vị trí trong locus đột biến ĐỘC LẬP nhau (không có
  thiên lệch cục bộ dọc locus) — lấy nguyên văn từ brief X8.

## Giới hạn (in lại nguyên văn trong `report.md`, mục CẢNH BÁO BẮT BUỘC)

1. **"Xuất hiện" khác "nhìn thấy được" khác "lan ra quần thể".** Công cụ chỉ tính xác suất MỘT bản sao
   nhiễm sắc thể quay lại đúng base gốc trong MỘT thế hệ. Nếu tính trạng đó LẶN, 1 bản quay ngược ở THỂ
   DỊ HỢP sẽ KHÔNG biểu hiện ra kiểu hình — cần cả 2 bản cùng mang allele gốc (đồng hợp) mới thấy được ở
   1 cá thể, và cần đột biến đó LAN RA (tăng tần số) qua nhiều thế hệ nữa mới có ý nghĩa hơn 1 cá thể
   đơn lẻ. Ba câu hỏi này là BA CÂU HỎI KHÁC NHAU — công cụ CHỈ trả lời câu đầu tiên (xem mục "Lưu ý khi
   k=1" ở trên, minh họa cụ thể sự khác biệt "xuất hiện" vs "cố định").
2. **Mô hình BỎ QUA**: chọn lọc (selection — giả định mọi phần tử trung tính), đột biến **MẤT ĐOẠN**
   (deletion — không quay ngược được bằng bất kỳ con đường nào ở đây), và **CHUYỂN ĐỔI GENE** (gene
   conversion — có thể thay đổi 1 base nhanh hơn đột biến điểm đơn thuần, không được mô hình hóa).
3. Nếu phần tử bị **XÓA HẲN** (deletion thật sự, không còn base nào ở vị trí đó) thì xác suất quay
   ngược qua CHÍNH con đường đột biến điểm này bằng **ĐÚNG SỐ 0**, không phải một số nhỏ — vì không còn
   gì để "đột biến trở lại đúng base gốc". Công cụ CHỈ mô hình đột biến ĐIỂM (thay 1 base bằng 1 base
   khác ở 1 vị trí VẪN CÒN tồn tại).

Ngoài ra: hằng số mặc định `mu=4.6e-9`, `kappa=2.0` (và `--n-diploid 1e6`, `--generation-years 2`) là
**VÍ DỤ MINH HỌA lấy NGUYÊN VĂN từ brief X8**, KHÔNG PHẢI số đã xác nhận riêng cho chào mào — dùng lệnh
`sweep` để không phụ thuộc 1 con số đinh, và chờ `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md`.

## Test

`tests/test_reversion_prob.py` — unittest thuần, **KHÔNG có `random.Random` nào cả** (toàn bộ mô hình
là công thức đóng đại số, không mô phỏng ngẫu nhiên — khác `dilution_sim.py`), nên không cần seed, kết
quả tất định 100%. Chạy:

```
python -B -m unittest discover -s pipeline/p06_atavism/tests -v
```

**41/41 test PASS** (đã xác nhận chạy 3 lần liên tiếp, giống hệt nhau cả 3 lần — tất định, không cần
kiểm flaky như test có random).

Bao phủ (tóm tắt; xem tên hàm test để biết chi tiết): `p_specific` chuyển đổi/đảo chuyển cộng đúng theo
kappa (`transition + 2*transversion = 1`), kappa=1 cho đúng 1/3 mọi đường; **phép thử BẮT BUỘC phân
biệt 2 con đường** (`TestSequentialVsConcurrentDependenceOnN`): gọi CÙNG 1 cặp hàm với N chênh nhau
10⁸ lần — `sequential` cho kết quả GIỐNG HỆT, `concurrent` giảm đúng TỈ LỆ NGHỊCH với N (kiểm tỉ lệ
chính xác, không chỉ "có khác nhau"); tăng `mu` giảm thời gian chờ đúng tỉ lệ nghịch (`sequential` ~
1/mu, `concurrent` ~ 1/mu**k — kiểm cả số mũ k); `compare` với locus 1 base, 1 phần tử cho tỉ số đúng
bằng `1/p_specific`; k=0 (`--mutated ""`, `--mutated-count 0`, `--k-list` chứa 0) đều `raise SystemExit`
rõ ràng ở cả cấp hàm lẫn cấp CLI, không trả về vô cực im lặng; ổn định số học ở k=10 (xem mục "Ổn định
số học" ở trên) gồm cả kiểm tra "tiền đề" (phép tính ngay thẳng THẬT SỰ lỗi/tràn) lẫn đối chiếu độc lập
bằng `decimal.Decimal`; `p_new` không hủy số về 0 ở `mu` cực nhỏ (tiền đề: phép tính ngay thẳng THẬT SỰ
trả về 0); lệnh `locus` gán đúng trạng thái `goc`/`dot_bien`, validate tên phần tử trùng/không tồn tại;
`sweep` sinh đúng `|mu-list|×|k-list|×|n-list|` dòng, cột `wait_gen_sequential` giống hệt qua các N khác
nhau trong khi `wait_gen_concurrent` khác nhau; `report` in đủ 3 cảnh báo bắt buộc + dòng "TẠM
THỜI"/"R14" + diễn giải đúng số dòng k≥2 có tuần tự ngắn hơn (kiểm với dữ liệu dựng sẵn biết trước đáp
số, không chỉ đếm chuỗi xuất hiện); `--help` mọi lệnh (kể cả top-level) thoát mã 0.

## `example/`

Chạy thật cả 4 lệnh, tham số theo đúng ví dụ CLI của brief X8:

```
python -B pipeline/p06_atavism/reversion_prob.py locus --elements A,B,C --mutated A \
    --mu 4.6e-9 --kappa 2.0 --n-diploid 1e6 --generation-years 2 --out example/locus_A.txt

python -B pipeline/p06_atavism/reversion_prob.py compare --locus-length 1000 \
    --mu 4.6e-9 --kappa 2.0 --mutated-count 1 --out example/compare_locus1000.tsv

python -B pipeline/p06_atavism/reversion_prob.py sweep --mu-list 1e-9,4.6e-9,1e-8 \
    --k-list 1,2,3,5 --n-list 1e4,1e6,1e8 --kappa 2.0 --generation-years 2 --out example/sweep.tsv

python -B pipeline/p06_atavism/reversion_prob.py report --sweep example/sweep.tsv --out example/report.md
```

Kết quả cụ thể (mu=4,6×10⁻⁹, kappa=2,0 → đường chuyển đổi, `p_specific`=0,5 → `p_rev`=2,3×10⁻⁹):

- `sequential_wait_years` k=1 ≈ **8,696×10⁸ năm**; k=3 ≈ **2,609×10⁹ năm** (2 thế hệ/năm, `example/locus_A.txt`
  cho k=1; k=3 tính từ cùng `p_rev` — xem `example/sweep.tsv` dòng `mu=4.6e-09, k=3`).
- `compare --locus-length 1000 --mutated-count 1` → `ratio_p_new_over_p_rev` = **2000** (giá trị đúng
  trước khi làm tròn 6 chữ số có nghĩa: 1999,9954 — xấp xỉ `L/p_specific = 1000/0,5`).
- `sweep.tsv`: 36 dòng dữ liệu = 3 `mu` × 4 `k` × 3 `N` (đúng khoảng brief yêu cầu: mu∈[1e-9,1e-8],
  k∈{1,2,3,5}, N∈{1e4,1e6,1e8}).

**Mọi con số trên đều dùng hằng số TẠM THỜI** (xem "Giới hạn") — chỉ minh họa cơ chế tính, KHÔNG PHẢI
kết luận định lượng cho chào mào thật cho tới khi có `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md`.
