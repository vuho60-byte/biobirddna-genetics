# Thuật toán và công thức — BioBirdDNA Genetics

Tài liệu này mô tả lõi toán học để model khác có thể tái triển khai mà không làm thay đổi semantics.

---

## 1. Quy ước

- Chim trống: `ZZ`.
- Chim mái: `ZW`.
- `N`: allele ký hiệu bình thường.
- `m`: allele ký hiệu biến đổi.
- `M`: mô hình di truyền.
- `G♂`, `G♀`: genotype cha/mẹ.
- `Y`: phenotype đời con.
- `D`: dữ liệu quan sát.

`N` và `m` không phải tên gene thật.

---

## 2. Phân bố giao tử

Với genotype autosomal hoặc ZZ:

```text
N/N → P(N)=1
N/m → P(N)=0.5, P(m)=0.5
m/m → P(m)=1
```

Với mái ZW:

```text
N/W → P(ZN)=0.5, P(W)=0.5
m/W → P(Zm)=0.5, P(W)=0.5
```

Prototype giả định segregation đều. Nếu có meiotic drive hoặc transmission distortion, phải thêm plugin khác và không sửa ngầm engine hiện tại.

---

## 3. Cross genotype chắc chắn

### 3.1. NST thường


a) Chọn giao tử từ cha và mẹ.  
 b) Kết hợp thành genotype đời con.  
 c) Nhân thêm xác suất giới.

\[
P(G_c, S_c) = P(a_{cha}) P(a_{mẹ}) P(S_c)
\]

Mặc định:

\[
P(S_c=mái)=P(S_c=trống)=0.5
\]

### 3.2. Liên kết Z

- Mẹ truyền `W` → con mái; Z của con mái đến từ cha.
- Mẹ truyền `Z` → con trống; con trống nhận một Z từ mỗi parent.

\[
P(G_c,S_c)=P(Z_{cha})P(Z/W_{mẹ})
\]

Không nhân sex ratio riêng vì sex đã phát sinh từ giao tử mẹ.

---

## 4. Genotype → phenotype

### 4.1. Lặn thường

```text
N/N → normal
N/m → normal visible; carrier state
m/m → target
```

### 4.2. Trội thường

```text
N/N → normal
N/m → target
m/m → target
```

### 4.3. Trội không hoàn toàn

```text
N/N → normal
N/m → intermediate
m/m → target
```

### 4.4. Lặn liên kết Z

Trống:

```text
N/N → normal
N/m → normal visible; carrier state
m/m → target
```

Mái:

```text
N/W → normal
m/W → target
```

### 4.5. Trội liên kết Z

Trống:

```text
N/N → normal
N/m hoặc m/m → target
```

Mái:

```text
N/W → normal
m/W → target
```

---

## 5. Penetrance

Nếu genotype thuộc state affected nhưng penetrance là `π`:

\[
P(Y=target|G)=\pi
\]

\[
P(Y=normal|G)=1-\pi
\]

Prototype dùng một penetrance chung. Production nên cho phép penetrance theo:

- genotype;
- giới;
- tuổi;
- mùa lông;
- môi trường;
- dòng gia đình.

Không dùng penetrance để che lấp chấm phenotype kém; hai nguồn bất định phải tách.

---

## 6. Parent genotype chưa chắc chắn

Khi cha/mẹ có posterior:

\[
P(G_{cha}|D_{cha}), P(G_{mẹ}|D_{mẹ})
\]

phân bố đời con:

\[
P(Y)=\sum_{G_{cha},G_{mẹ},G_c}
P(Y|G_c,M)
P(G_c|G_{cha},G_{mẹ},M)
P(G_{cha}|D_{cha})
P(G_{mẹ}|D_{mẹ})
\]

Code phải normalize posterior đầu vào và lỗi nếu tổng bằng 0.

---

## 7. Likelihood đời con

Với các category `k=1..K`, số quan sát `n_k`, xác suất model dự đoán `p_k`:

\[
P(D|M,G_{cha},G_{mẹ}) =
\frac{N!}{\prod_k n_k!}\prod_k p_k^{n_k}
\]

Log-likelihood:

\[
\log L = \log N! - \sum_k \log n_k! + \sum_k n_k \log p_k
\]

Code dùng `lgamma(n+1)` để tránh overflow.

---

## 8. Observation error

Để tránh một lỗi chấm duy nhất làm probability bằng 0 tuyệt đối, prototype dùng:

\[
\tilde{p}_k = (1-\epsilon)p_k + \frac{\epsilon}{K}
\]

Sau đó normalize lại trên category set.

`ε` không được đặt lớn để cứu model sai. Nó phải được version/config và sau này estimate từ reviewer disagreement.

---

## 9. Marginalize parent configuration

Nếu có nhiều genotype cha/mẹ tương thích phenotype:

\[
P(D|M)=\sum_c P(D|M,c)P(c|M)
\]

Trong log-space:

\[
\log P(D|M)=\operatorname{logsumexp}
(\log P(D|M,c)+\log P(c|M))
\]

Prototype dùng prior đều giữa configuration khi không có prior riêng.

---

## 10. Model posterior

\[
P(M|D)=\frac{P(D|M)P(M)}{\sum_jP(D|M_j)P(M_j)}
\]

Report phải nêu rõ:

- model set được so sánh;
- prior;
- observation error;
- parent phenotype assumptions;
- dữ liệu bị thiếu.

Posterior chỉ tương đối trong model set. Nếu model thật không nằm trong tập, posterior cao vẫn có thể gây hiểu nhầm.

---

## 11. Non-genetic baseline

Prototype `M0_NONGENETIC` dùng target probability nền `p0` độc lập parent/genotype:

\[
P(Y=target)=p_0
\]

Đây chưa phải model sinh học hoàn chỉnh. Nó là guard để model comparison không bắt buộc mọi data phải mang tính Mendel.

Production nên thay bằng các candidate:

- age-dependent;
- injury-localized;
- progressive greying;
- environment-associated;
- random classification error.

---

## 12. Relationship matrix

Numerator relationship matrix `A`:

Với founder không biết parent:

\[
A_{ii}=1
\]

Với cá thể `i` có sire `s`, dam `d`:

\[
A_{ij}=\frac{1}{2}(A_{sj}+A_{dj}), \quad j<i
\]

\[
A_{ii}=1+\frac{1}{2}A_{sd}
\]

Nếu một parent không biết, contribution parent đó là 0 trong prototype.

### 12.1. Inbreeding cá thể

\[
F_i=A_{ii}-1
\]

### 12.2. Inbreeding đời con dự kiến

\[
F_{con}=\frac{1}{2}A_{cha,mẹ}
\]

Founders “không biết” được coi không liên quan là giả định tiện dụng, có thể đánh giá thấp cận huyết. Report production phải gắn flag `founder_relatedness_unknown`.

---

## 13. Expected information gain

Với model hypothesis `M` có weight `w_m` và outcome distribution `P(Y|M,pair)`:

Mixture:

\[
P(Y)=\sum_m w_mP(Y|M,pair)
\]

Expected information gain cho một offspring:

\[
EIG = H(P(Y)) - \sum_m w_mH(P(Y|M,pair))
\]

Trong đó:

\[
H(P)=-\sum_yP(y)\log_2P(y)
\]

EIG cao khi các model dự đoán outcome khác nhau. Nó không đồng nghĩa cặp đó tốt về sức khỏe hoặc đạo đức.

---

## 14. Genotype uncertainty

Với posterior genotype `q_g`:

\[
H_G=-\sum_g q_g\log_2q_g
\]

Normalize:

\[
U_G=\frac{H_G}{\log_2 |G|}
\]

`U_G` nằm `[0,1]` khi có trên một genotype state.

---

## 15. Pair score

Prototype:

\[
Score =
 w_TP(target)
 + w_IEIG
 - w_FF_{con}
 - w_HRisk
 - w_UU_G
\]

Thứ tự xử lý:

1. hard block cùng cá thể (`same_individual_pair`) và relation/welfare;
2. tính score;
3. sort eligible trước;
4. giữ blocked pair trong report.

`same_individual_pair` là block vô điều kiện, không phụ thuộc
`max_offspring_inbreeding` hoặc cờ cho phép bỏ block quan hệ. Không được “phạt
nhẹ” cha–con thay vì block.

---

## 16. Rule engine sắc tố

Feature đầu vào là số đã đo hoặc reviewer coding:

- white area ratio;
- brown/grey/dark shift;
- red-domain change;
- pattern retained;
- iris red confirmed;
- flash risk;
- progressive across molts;
- first plumage;
- injury localization;
- Q-level.

Rule cộng score rồi normalize:

\[
P_{heuristic}(label)=\frac{s_{label}}{\sum_js_j}
\]

Đây không phải posterior sinh học. Tên trường `scores` được dùng thay vì `probability_of_gene`.

---

## 17. Monte Carlo tương lai

Analytic cross đủ cho one-locus. Monte Carlo chỉ cần khi:

- multi-locus lớn;
- uncertain survival;
- family-level random effects;
- complex age penetrance;
- distribution of brood outcomes.

Mỗi run phải lưu seed. Không dùng Monte Carlo nếu analytic exact có thể tính.

---

## 18. Confidence interval cho một lứa tương lai

Nếu target probability là `p`, số offspring `n`:

\[
K\sim Binomial(n,p)
\]

Có thể xuất distribution `P(K=k)`, không hứa “n × p con”.

Khi nhiều phenotype category:

\[
\mathbf{K}\sim Multinomial(n,\mathbf{p})
\]

---

## 19. Survival model tương lai

Cần tách:

\[
P(phenotype\ observed)=
P(genotype\ at\ conception)
\times P(survival|genotype)
\times P(classified|phenotype)
\]

Nếu chỉ dùng chim sống sau thay lông, không thể mặc định tỷ lệ đó bằng segregation lúc thụ tinh.

---

## 20. Acceptance tests toán học

- mọi distribution sum gần 1;
- AR `N/m × m/m` → target 0.5;
- ZR `m/m × N/W` → female target 0.5, male normal 0.5;
- full siblings A=0.5, offspring F=0.25;
- unrelated founders F=0;
- posterior models sum 1;
- evidence guard chặn gene claim ở G4.

---

## 21. Synthetic recovery, calibration và absolute fit

Với model sinh dữ liệu `M*`, brood size `n` và seed đã khóa, mỗi replicate lấy:

\[
\mathbf{K}_r \sim Multinomial(n,\mathbf{p}_{M*})
\]

Sau đó chạy chính inference engine trên `K_r`. Recovery top-1 chỉ được tính khi
true model là unique posterior maximum. Nếu nhiều model bằng nhau trong tolerance,
report giữ toàn bộ `best_models` và không dùng thứ tự candidate array để phá tie.

Wilson 95% interval cho recovery rate:

\[
\frac{\hat p + z^2/(2R) \pm z
\sqrt{\hat p(1-\hat p)/R + z^2/(4R^2)}}{1+z^2/R}
\]

với `z = 1.95996398454`.

Multiclass Brier score:

\[
BS_r=\sum_j (P(M_j|D_r)-I[M_j=M*])^2
\]

Calibration chỉ hợp lệ trên ensemble design có nhiều realized true models. Trong
mỗi bin posterior của model `j`, so `mean_posterior` với tỷ lệ
`I[M_j=M*]`; top-1 accuracy chỉ là discrimination metric, không phải calibration
target.

Absolute-fit deviance:

\[
D=2\sum_{k:n_k>0}n_k\log\frac{n_k}{Np_k}
\]

Với composite candidate model, chọn cấu hình bố mẹ tốt nhất trên dữ liệu quan sát,
sinh bootstrap từ cấu hình fitted đó, rồi **refit lại cấu hình** ở từng bootstrap
replicate trước khi tính `D_b`. Add-one p-value:

\[
p_{boot}=\frac{1+\sum_b I[D_b\ge D_{obs}]}{B+1}
\]

Do đó phải có `B + 1 > 1/alpha`; nếu không, warning ở alpha đã chọn là bất khả thi.
Chỉ cảnh báo model set misspecified khi mọi candidate, gồm `M0_NONGENETIC`, đều
có `p_boot < alpha`.

Seed dùng SHA-256 của base seed, experiment ID, brood size và replicate index rồi
đọc 64 bit đầu theo big-endian và mask về 53 bit. Cả Python và Node dùng LCG
`state=(6364136223846793005*state+1442695040888963407) mod 2^53`; multinomial
bucket được sort theo byte UTF-8. Bootstrap thêm suffix `goodness_of_fit` vào cùng
seed material. Vì vậy cùng request và engine version tạo cùng stream/output ở hai
runtime. Work-budget validation chạy trước mọi sampling.
Known built-in model IDs là immutable semantic identifiers: nếu generator mang
một ID built-in thì snapshot của nó phải khớp model canonical dù ID đó không nằm
trong candidate set. Generator cố ý nằm ngoài model set phải có ID tùy biến.
