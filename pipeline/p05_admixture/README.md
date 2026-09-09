# p05_admixture — brief X7: mô phỏng "pha nước đường" (introgression dilution)

Bối cảnh đầy đủ: `research/briefs/X7-mo-phong-pha-loang-introgression.md`. Tham khảo phong cách CLI và
cách ghi TSV từ `pipeline/p04_anomaly/composition_scan.py` (không import module đó — mỗi lối trong
`pipeline/` tự chứa, đúng quy ước hiện có của dự án).

## Câu hỏi khoa học

Chủ dự án mô tả một mô hình: một bộ gene "gốc" trộn với một bộ gene "bản địa" rồi loãng dần qua các
đời lai ngược, giống pha nước đường, giống khoai lang Nhật qua F1→F2→F3→F4. Câu hỏi trả lời được
bằng mô phỏng: **sau bao nhiêu đời thì một sự kiện lai trở nên KHÔNG THỂ phát hiện, và cái gì quyết
định ngưỡng đó?**

Đây là mô hình chuẩn của introgression. Hai đại lượng suy giảm theo hai cách khác nhau:

- **Tỉ lệ tổ tiên** giảm theo cấp số nhân dưới lai ngược liên tiếp (mỗi đời còn lại ~một nửa), nhưng
  nếu quần thể lai giao phối nội bộ (không lai ngược) thì tỉ lệ trung bình toàn quần thể **không đổi**
  về mặt kỳ vọng, chỉ dao động do trôi dạt.
- **Chiều dài đoạn tổ tiên** giảm do tái tổ hợp cắt nhỏ dần, ở **cả hai** chế độ — đây là "đồng hồ"
  phân tử: đoạn dài = lai gần đây, đoạn ngắn = lai xa xưa.

## Chạy — `dilution_sim.py`, ba lệnh con

```
python -B pipeline/p05_admixture/dilution_sim.py simulate --generations 20 --mode backcross \
    --init-fraction 0.5 --chrom-length-cm 100 --n-chrom 30 --replicates 100 --seed 42 \
    --out tracts.tsv --stats sim_stats.tsv

python -B pipeline/p05_admixture/dilution_sim.py detect --stats sim_stats.tsv \
    --min-tract-cm 5 --min-fraction 0.01 --out detect.tsv

python -B pipeline/p05_admixture/dilution_sim.py report --stats sim_stats.tsv --detect detect.tsv \
    --out report.md
```

Chỉ Python 3.11 stdlib (`argparse`, `math`, `random`, `statistics`). Không numpy/scipy, không docker,
không tải dữ liệu mạng.

### 1. `simulate` — Markov trên mỗi nhiễm sắc thể

Mỗi replicate × mỗi nhiễm sắc thể là **1 cá thể lưỡng bội thật sự**, giữ **hai** bản sao (homolog A,
homolog B), mỗi bản sao 1 danh sách đoạn có nguồn gốc "gốc". Đời 0 khởi tạo **xác định** (không ngẫu
nhiên) kiểu F1 khi `--init-fraction 0.5` (homolog A = toàn bộ "gốc", homolog B = toàn bộ "bản địa");
tổng quát cho `--init-fraction` bất kỳ bằng cách lấp đầy homolog A trước rồi mới sang homolog B (xem
`_init_individual`) — nhờ vậy `--generations 0` trả về **đúng bằng** `--init-fraction`, không chỉ xấp
xỉ thống kê.

Mỗi đời, **mỗi cá thể tự giảm phân** (meiosis): sinh số điểm tái tổ hợp theo Poisson (kỳ vọng
`chrom_length_cm/100` giao chéo/nhiễm sắc thể — thuật toán Knuth, tích lũy tích các số ngẫu nhiên đều
đến khi vượt `e^-λ`, `O(λ)` nên chỉ phù hợp λ nhỏ như thang cM thông thường; xem docstring
`_sample_poisson`), với mỗi khoảng giữa 2 điểm chéo kế tiếp **chọn ngẫu nhiên (đồng xu công bằng) một
trong hai bản sao của chính cá thể đó** — đây mới đúng nghĩa đen "chọn ngẫu nhiên một trong hai bản
sao" của brief. Giao tử đó ghép với 1 giao tử thứ hai để tạo thế hệ kế tiếp:

- **`backcross`**: giao tử thứ hai luôn của 1 cá thể bản địa THUẦN (2 homolog rỗng) → giao tử chắc
  chắn rỗng, không cần giảm phân thật.
- **`isolated`**: giao tử thứ hai là sản phẩm giảm phân **độc lập** của **một cá thể khác** (`mate`)
  rút ngẫu nhiên đều từ chính quần thể đang mô phỏng ở thế hệ hiện tại (loại trừ chính nó khi
  `--replicates > 1`; khi `--replicates 1` phải tự phối — trường hợp suy biến, xem Giới hạn).

Xuất `tracts.tsv`: `replicate generation chrom start_cm end_cm len_cm` (đoạn của **cả hai** homolog
đều được ghi, không có cột phân biệt homolog — đúng schema brief; `replicate`/`chrom` đánh số từ 1).
Xuất `sim_stats.tsv`, mỗi đời 1 dòng: `generation mean_ancestry_fraction sd_ancestry_fraction
n_tracts mean_tract_cm median_tract_cm max_tract_cm`. `mean_ancestry_fraction`/`sd_ancestry_fraction`
tính trên **ti lệ trung bình mỗi replicate** (mỗi replicate = trung bình trên toàn bộ `n_chrom` × 2
homolog của nó), rồi lấy trung bình/độ lệch chuẩn (population, `statistics.pstdev`) của các giá trị
đó qua các replicate — đo đúng "trôi dạt giữa các dòng lai", không phải phương sai trong-1-cá-thể.
`n_tracts`/`mean_tract_cm`/`median_tract_cm`/`max_tract_cm` gộp tất cả đoạn của mọi replicate/nhiễm
sắc thể/homolog ở đời đó (0.0 mặc định nếu không còn đoạn nào).

**Hằng đẳng thức kiểm chứng được** (test `test_tracts_within_bounds_and_length_matches_fraction`):
tổng `len_cm` mọi dòng ở 1 đời = `mean_ancestry_fraction × replicates × n_chrom × chrom_length_cm × 2`
(nhân 2 vì mỗi nhiễm sắc thể có 2 homolog) — trong sai số làm tròn 6 chữ số có nghĩa của `_fmt_cell`.

### 2. `detect` — đánh giá còn "phát hiện được" hay không

Đọc `sim_stats.tsv`. Mỗi đời: `fraction_ok=1` nếu `mean_ancestry_fraction >= --min-fraction`;
`tract_ok=1` nếu `max_tract_cm >= --min-tract-cm` (cả hai dùng ngưỡng **kèm theo**, `>=`, không phải
`>` — quy ước của script, giống `outliers` của `composition_scan.py`, vì biến liên tục nên xác suất
trúng đúng ngưỡng gần như 0). `detectable=1` **chỉ khi cả hai** đều =1. Xuất `detect.tsv`:
`generation fraction_ok tract_ok detectable`.

**QUAN TRỌNG — ngưỡng phát hiện:** `--min-fraction` và `--min-tract-cm` là **THAM SỐ NGƯỜI DÙNG ĐẶT**,
**KHÔNG PHẢI hằng số tự nhiên**. Giá trị thực tế (bao nhiêu cM là "còn phát hiện được" với phương
pháp phân tích nào) phải lấy từ tài liệu phương pháp thật. `research/raw/R11-introgression-ghost-
lineage.md` (brief R11, **đã xong**, không phải tôi làm) đã tìm — nhưng **không tìm ra được 1 ngưỡng
phổ quát đã kiểm chứng trực tiếp**: chỉ có mốc tham khảo **UNVERIFIED** (tổng hợp Perplexity, chưa tự
mở từng bài trích dẫn gốc) — "tỉ lệ pha trộn ~1-5% có thể phát hiện ở độ sâu thời gian vừa phải, công
suất giảm mạnh khi <1-2%"; nguyên lý định tính (đánh đổi 2 trục tỉ lệ × thời gian) là STRONG_INFERENCE
tốt hơn. R11 cũng nói rõ: **không tìm được con số riêng cho kịch bản lai ngược liên tục** (đúng kịch
bản chế độ `backcross` của script này — chỉ có số cho 1 đợt pha trộn đơn/pulse). Cho tới khi có số đã
kiểm chứng trực tiếp, mọi kết luận từ `detect` (và từ `report`, vốn đọc lại `detect.tsv`) **chỉ mang
tính minh họa**, không dùng để phát biểu định lượng cho trường hợp chào mào thật. Ví dụ trong
`example/` dùng `--min-fraction 0.01 --min-tract-cm 5` chỉ để MINH HỌA cơ chế `detect`, lấy đại khái
trong khoảng mốc tham khảo UNVERIFIED nói trên — không phải số đã được xác nhận.

### 3. `report` — bảng theo đời + diễn giải + cảnh báo bắt buộc

Đọc `sim_stats.tsv` + `detect.tsv` (nối theo cột `generation`), xuất `report.md` gồm: **cảnh báo bắt
buộc** (không có cách nào tắt), phần **diễn giải tự động từ chính dữ liệu** (không hardcode: tỉ lệ
đầu/cuối, tỉ số cuối/đầu, phân loại xu hướng "giảm mạnh kiểu backcross" hay "phẳng kiểu isolated" hay
"chưa rõ", đời `detectable=1` xa nhất tìm được), đoạn giải thích khái niệm backcross vs isolated khác
nhau ra sao (văn bản cố định, không phụ thuộc dữ liệu — vì `report` chỉ nhận 1 cặp `--stats`/`--detect`
tại một thời điểm, không biết trước đó là chế độ nào), và bảng đầy đủ theo đời.

## Mô hình và giả định

**Vì sao lưỡng bội (2 homolog/cá thể) chứ không phải 1 danh sách đoạn/cá thể đơn giản hơn:** bản nháp
đầu tiên của script này dùng mô hình đơn giản hơn — mỗi replicate/nhiễm sắc thể chỉ giữ **1** danh
sách đoạn (khớp sát nghĩa đen "tập CÁC ĐOẠN" số ít của brief). Mô hình đó **đủ đúng cho `backcross`**
(đã kiểm chứng toán học: kỳ vọng tỉ lệ tại 1 vị trí tuân theo đúng phép truy hồi nhân 0,5 mỗi đời) —
nhưng khi tự chạy thử phát hiện **`isolated` bị ĐÓNG BĂNG vĩnh viễn** nếu quần thể khởi tạo đồng nhất
(mọi cá thể giống hệt nhau ở đời 0, đúng như cách khởi tạo xác định cần cho test `--generations 0`):
giảm phân giữa 2 "bản sao" GIỐNG HỆT NHAU (self và mate, khi mate chỉ là 1 bản sao y hệt của quần thể
đồng nhất) luôn trả về đúng chính nó bất kể điểm chéo/đồng xu ra sao — quần thể không bao giờ tiến
triển, vi phạm trực tiếp yêu cầu "chiều dài đoạn trung bình giảm đơn điệu ở CẢ HAI chế độ". Đã sửa
bằng cách chuyển sang mô hình lưỡng bội thật (2 homolog/cá thể): homolog A và homolog B của MỘT cá
thể khác nhau ngay từ đời 0 (kiểu F1), nên giảm phân (tái tổ hợp giữa A và B của chính nó) luôn sinh
ra biến dị thật giữa các replicate ngay từ đời 1, dù quần thể khởi tạo đồng nhất — cả 2 chế độ đều
tiến triển đúng. Regression test cho đúng lỗi này: `tests/test_dilution_sim.py::TestMeiosis::
test_identical_homologs_reproduce_themselves_exactly`.

**Vì sao vẫn đúng toán học đã nêu ở trên** (suy ra từ cấu trúc mô hình, tự kiểm chứng bằng mô phỏng
trong test — brief yêu cầu "tự kiểm chứng công thức đúng bằng chính mô phỏng thay vì tin dòng có sẵn"):
xét 1 vị trí `x` cố định trên nhiễm sắc thể, gọi `p_g(x)` = xác suất `x` còn "gốc" ở đời `g` (trung
bình trên 1 bản sao ngẫu nhiên của 1 cá thể ngẫu nhiên trong quần thể). Đời `g+1`, giao tử tự thân
(self-gamete) giữ trạng thái `x` với xác suất đúng bằng `p_g(x)` (vì giảm phân chỉ chọn LẠI giữa 2
homolog có sẵn, không đổi giá trị). `backcross`: giao tử thứ hai luôn "bản địa" (0) → `p_{g+1}(x) =
0,5·p_g(x) + 0,5·0 = 0,5·p_g(x)` — suy giảm cấp số nhân đúng `init_fraction × 0,5^g`. `isolated`:
giao tử thứ hai là sản phẩm giảm phân của 1 cá thể KHÁC nhưng **cùng phân phối** (đối xứng giữa các
replicate) → `p_{g+1}(x) = 0,5·p_g(x) + 0,5·p_g(x) = p_g(x)` — hằng số theo kỳ vọng (martingale),
đúng "không đổi, chỉ dao động do trôi dạt". Cả hai đã **kiểm chứng thực nghiệm** trong
`tests/test_dilution_sim.py` (200 replicate × 20 nhiễm sắc thể, seed cố định).

**Các giả định khác:**
- Poisson bằng thuật toán Knuth (tích lũy tích số ngẫu nhiên đều) — thuật toán phổ biến trong nhiều
  giáo trình mô phỏng; **KHÔNG tự tra cứu được nguồn/DOI cụ thể trong phiên viết script này** vì brief
  cấm tải dữ liệu mạng, nên không gán một trích dẫn ngoài cụ thể nào cho CHÍNH thuật toán rút mẫu —
  chỉ dựa vào việc đây là thuật toán chuẩn/phổ biến, và kiểm chứng gián tiếp qua test
  (`test_empirical_mean_close_to_lambda`: trung bình mẫu ~5000 lần rút khớp λ trong sai số chuẩn).
- Công thức tỉ lệ tổ tiên `init_fraction × 0,5^đời` dưới backcross: **tự suy ra độc lập** từ cấu trúc
  mô hình (xem chứng minh ở trên) + xác nhận thực nghiệm trong test. Sau khi suy ra, đối chiếu lại với
  nghiên cứu đã có sẵn trong dự án — **research/raw/R11-introgression-ghost-lineage.md** (brief R11,
  đã hoàn thành, KHÔNG sửa tay) — cho thấy CÙNG công thức `(1/2)^n` đã được ghi nhận
  **[CONFIRMED]**: "Mendel cơ bản, hình thức hóa định lượng trong Gravel, S. 2012. 'Population
  genetics models of local ancestry.' *Genetics* 191(2):607-19, DOI
  [10.1534/genetics.112.139808](https://doi.org/10.1534/genetics.112.139808) (PMID 22491189, R11 ghi
  đã đọc trực tiếp)" — khớp đúng kết quả tự suy ra ở trên, không mâu thuẫn.
- Công thức chiều dài đoạn kỳ vọng `~ 1/(đời × tỉ lệ tái tổ hợp)` Morgan trong mô tả câu hỏi khoa học
  ở trên (dạng kỳ vọng của phân phối Mũ với tốc độ = đời × tỉ lệ tái tổ hợp): **lấy nguyên văn từ brief
  X7**. Tự bản thân phiên viết script này **KHÔNG tra cứu được nguồn** (cấm tải mạng) — nhưng
  **research/raw/R11-introgression-ghost-lineage.md** (đã hoàn thành trước, không phải tôi làm) ghi
  nhãn **[CONFIRMED]** cho đúng dạng công thức này ("chiều dài đoạn tổ tiên ~ phân phối Mũ, tốc độ
  phân rã ∝ số thế hệ kể từ pha trộn × tần số tái tổ hợp"), nguồn: Gravel 2012 (trên) + Corbett-Detig,
  R. & Nielsen, R. 2017. *PLoS Genetics* 13(1):e1006529, DOI
  [10.1371/journal.pgen.1006529](https://doi.org/10.1371/journal.pgen.1006529) (PMID 28045893, R11
  ghi đã đọc trực tiếp). Script này (X7) tự nó chỉ kiểm chứng được **xu hướng giảm đơn điệu** theo đời
  bằng mô phỏng (`test_mean_tract_length_decreases_monotonically_in_both_modes`), KHÔNG tự kiểm chứng
  hằng số tỉ lệ chính xác của công thức — độ tin cậy CONFIRMED ở trên là **kế thừa từ R11**, không phải
  do X7 tự xác minh lại.
- `--replicates 1` ở chế độ `isolated`: không còn cá thể nào khác để ghép, tự phối với chính nó
  (`_pick_mate_index` trả về chính `r`) — trường hợp suy biến, chỉ để không crash, không có ý nghĩa
  sinh học thật.
- Mỗi lần cần 1 giao tử (của self hay của mate) đều gọi giảm phân **độc lập mới** — giống 1 cá thể
  thật sinh nhiều giao tử độc lập qua nhiều lần giảm phân, kể cả khi cùng 1 cá thể được nhiều cá thể
  khác cùng rút làm `mate` trong cùng 1 đời.
- `chrom_length_cm` là **một giá trị duy nhất dùng chung cho mọi nhiễm sắc thể** (không có bản đồ
  riêng từng nhiễm sắc thể) — xem Giới hạn.

## Giới hạn

Mô hình mô phỏng này **BỎ QUA**:

1. **Chọn lọc** (selection) — không có đoạn tổ tiên nào được ưu tiên giữ lại hay đào thải; mọi vị trí
   trên nhiễm sắc thể trung tính như nhau.
2. **Thay đổi kích thước quần thể** (population size change) — `--replicates` cố định suốt mô phỏng,
   không có thắt cổ chai (bottleneck) hay tăng trưởng quần thể thật như quần thể chim ngoài tự nhiên.
3. **Biến thiên tốc độ tái tổ hợp dọc nhiễm sắc thể** (recombination rate heterogeneity) — dùng đúng
   một tỉ lệ trung bình `chrom_length_cm/100` giao chéo/nhiễm sắc thể/đời, không có "điểm nóng" hay
   vùng ít tái tổ hợp gần tâm động (centromere) như nhiễm sắc thể thật.
4. **Bản đồ di truyền đồng nhất** (uniform genetic map) — 1 cM tương ứng xác suất tái tổ hợp như nhau
   ở mọi vị trí, không khớp bản đồ di truyền thực của bất kỳ loài chim nào (kể cả chào mào).

Ngoài 4 điều bắt buộc nêu trên (in lại nguyên văn trong `report.md`), script còn:

- Không mô hình hóa **liên kết giữa các nhiễm sắc thể** (mỗi nhiễm sắc thể mô phỏng hoàn toàn độc
  lập) và không có nhiễm sắc thể giới tính riêng (Z/W) — mọi nhiễm sắc thể xử lý như thường (autosome).
  không có drift/di cư/đột biến ngoài khung tái tổ hợp+ancestry đã mô tả.
- `--replicates` không tương đương kích thước quần thể thực (không có chồng lấp thế hệ, không có số
  con khác nhau giữa các cặp, không có giới tính ảnh hưởng tới ai ghép với ai).
- Kết quả chỉ cho **BẬC ĐỘ LỚN** (order of magnitude) về số đời cần để một sự kiện lai trở nên khó
  phát hiện — **KHÔNG PHẢI con số dùng để công bố hay kết luận định lượng** cho trường hợp chào mào
  thật.
- Ngưỡng phát hiện của `detect` là tham số người dùng đặt (xem mục `detect` ở trên) — R11 (đã xong)
  không tìm ra được số phổ quát đã kiểm chứng trực tiếp, chỉ có mốc tham khảo UNVERIFIED.

## Test

`tests/test_dilution_sim.py` — unittest thuần, seed cố định, phần lớn dữ liệu nhúng ngay trong file.
Chạy:

```
python -B -m unittest discover -s pipeline/p05_admixture/tests -v
```

37/37 test PASS (đã xác nhận chạy 3 lần liên tiếp, không flaky — mọi test dùng seed cố định).

Bao phủ (tóm tắt; xem tên hàm test để biết chi tiết): `_clip_tracts`/`_merge_tracts` đơn vị;
`_sample_poisson` không âm luôn, trung bình mẫu khớp λ trong sai số chuẩn; `_init_individual` đúng
tỉ lệ trung bình cho cả 2 nhánh (`<=0,5` và `>0,5`) và 2 biên (0,0/1,0); `_meiosis` xác định theo seed,
**hồi quy trực tiếp cho lỗi "đóng băng" đã tìm thấy** (2 bản sao giống hệt nhau → giữ nguyên chính
nó), và trung bình thống kê ~một nửa chiều dài khi 1 bản sao đầy/1 bản sao rỗng; `_pick_mate_index`
không bao giờ tự chọn chính nó khi quần thể >1. Ở mức lệnh (`simulate`/`detect`/`report` qua `main()`):
**backcross giảm đơn điệu và khớp cấp số nhân `init_fraction × 0,5^đời`** trong sai số thống kê (6 lần
sai số chuẩn của trung bình, ≥200 replicate); **isolated KHÔNG giảm hệ thống** (nằm trong dải rộng
quanh `init_fraction`) **VÀ khác biệt rõ rệt (≥20 lần) so với dự đoán backcross ở cùng đời** — đây là
phép thử BẮT BUỘC phân biệt 2 chế độ theo brief, và sẽ THẬT SỰ FAIL nếu logic 2 chế độ bị hoán đổi
hay isolated vô tình suy giảm giống backcross; chiều dài đoạn trung bình giảm đơn điệu ở CẢ HAI chế
độ; không đoạn nào vượt `[0, chrom_length_cm]`; tổng `len_cm` khớp `mean_ancestry_fraction × tổng độ
dài (× 2 vì lưỡng bội)` trong sai số làm tròn 6 chữ số có nghĩa; cùng `--seed` cho kết quả giống hệt
(so cả 2 file) qua 2 lần chạy, seed khác cho kết quả khác; `--generations 0` trả đúng `--init-fraction`
(kiểm cả 2 nhánh `_init_individual`) với `sd_ancestry_fraction=0` chính xác; đối số không hợp lệ
(`--generations -1`, `--init-fraction 1.5`, `--replicates 0`, `--n-chrom 0`, `--chrom-length-cm 0`)
đều raise `SystemExit`; `detect` áp đúng ngưỡng kèm theo (`>=`) ở biên; `report` tự đọc đúng xu hướng
dữ liệu (phân biệt được kịch bản dựng sẵn kiểu "giảm mạnh" vs "phẳng"), luôn in cảnh báo bắt buộc +
giải thích 2 chế độ; `--help` mọi lệnh (kể cả top-level) thoát mã 0.

## `example/`

Kết quả chạy thật `simulate` (+ `detect` + `report` để minh họa đủ 3 lệnh) cho cả 2 chế độ, 20 đời:

```
python -B pipeline/p05_admixture/dilution_sim.py simulate --generations 20 --mode backcross \
    --replicates 40 --n-chrom 12 --out example/backcross_tracts.tsv --stats example/backcross_stats.tsv
```

(và tương tự cho `isolated`, rồi `detect --min-tract-cm 5 --min-fraction 0.01` + `report` trên mỗi
`*_stats.tsv`). Dùng `--replicates 40 --n-chrom 12` (nhỏ hơn mặc định 100/30 của CLI) **chỉ để giữ
tổng dung lượng `example/` dưới 5 MB** — mô phỏng `isolated` tích lũy càng nhiều đoạn (tái tổ hợp
không ngừng) chứ không mất ancestry như `backcross`, nên `tracts.tsv` của nó lớn hơn nhiều (~1 MB ở
quy mô này; ở mặc định 100×30 sẽ ra ~7 MB, vượt ngân sách) — công thức/kết luận không đổi, chỉ ít
lineage hơn nên nhiễu thống kê cao hơn một chút. Muốn chạy quy mô đầy đủ như brief minh họa (100
replicate × 30 nhiễm sắc thể), bỏ `--replicates`/`--n-chrom` để dùng mặc định.

Với tham số mặc định (`--init-fraction 0.5`), mô phỏng thực tế (200 replicate × 20 nhiễm sắc thể,
seed 42 — xem `tests/test_dilution_sim.py`) cho thấy tỉ lệ tổ tiên trung bình tụt xuống dưới 1% lần
đầu ở **đời 6** trong chế độ `backcross` (đời 5: 1,42%; đời 6: 0,78%) — khớp dự đoán lý thuyết
`0,5^7 = 0,78% < 1% < 0,5^6 = 1,56%`. Đây là số của MÔ HÌNH MINH HỌA (xem Giới hạn), không phải số
sinh học thật cho chào mào.
