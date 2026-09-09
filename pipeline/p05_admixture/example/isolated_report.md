# Bao cao mo phong pha loang gene 'goc' qua cac doi lai (X7 dilution_sim)

Nguon: `example/isolated_stats.tsv` + `example/isolated_detect.tsv` -- 21 dong du lieu.

## CANH BAO BAT BUOC -- doc truoc khi dung so trong bao cao nay cho bat ky ket luan nao

Mo hinh mo phong nay (`dilution_sim.py`) **BO QUA**:

1. **Chon loc** (selection) -- khong co doan to tien nao duoc uu tien giu
   lai hay dao thai; moi vi tri tren nhiem sac the trung tinh nhu nhau.
2. **Thay doi kich thuoc quan the** (population size change) -- `--replicates`
   co dinh suot mo phong, khong co that co (bottleneck) hay tang truong quan
   the thuc nhu quan the chim ngoai tu nhien.
3. **Bien thien toc do tai to hop doc nhiem sac the** (recombination rate
   heterogeneity) -- dung DUNG MOT ti le trung binh `chrom_length_cm/100`
   giao cheo/nhiem sac the/doi, khong co "diem nong" hay vung it tai to hop
   gan tam dong (centromere) nhu nhiem sac the that.
4. **Ban do di truyen dong nhat** (uniform genetic map) -- 1 cM tuong ung
   xac suat tai to hop nhu nhau o moi vi tri, khong khop ban do di truyen
   thuc cua bat ky loai chim nao (kem ca chao mao).

Vi vay ket qua trong bao cao nay chi cho **BAC DO LON** (order of magnitude)
ve so doi can de mot su kien lai tro nen kho phat hien -- **KHONG PHAI con
so dung de cong bo hay ket luan dinh luong** cho truong hop that. Nguong
`--min-fraction` / `--min-tract-cm` cua lenh `detect` la **THAM SO NGUOI
DUNG DAT**, khong phai hang so sinh hoc/tu nhien -- gia tri thuc te (bao
nhieu cM la "con phat hien duoc" voi phuong phap nao) phai lay tu tai lieu
phuong phap that. `research/raw/R11-introgression-ghost-lineage.md` (brief
R11, DA XONG) khong tim ra 1 nguong pho quat da kiem chung truc tiep -- chi
co moc tham khao UNVERIFIED (tong hop Perplexity, chua tu mo tung bai goc);
cho toi khi co so DA KIEM CHUNG, cot `detectable` chi mang tinh **minh hoa**.

## Dien giai (tu chinh du lieu trong file nay)

Du lieu chay tu doi 0 den doi 20. Ti le to tien trung binh di tu 0.5 (doi 0) xuong 0.479608 (doi 20) -- ti so cuoi/dau = 0.959216.
Xu huong nay **khong giam he thong** (con lai gan muc ban dau, chi dao dong) -- phu hop kieu 'giao phoi noi bo' (isolated).
Voi nguong da chon o lenh `detect` (xem CANH BAO BAT BUOC ben tren), su kien lai **con phat hien duoc toi da den doi 20** trong du lieu nay (doi `detectable=1` xa nhat); tu doi sau do trong du lieu nay khong con dong nao `detectable=1`.

## Backcross vs isolated -- khac nhau ra sao

Hai che do chi khac nhau o **ban phoi (mate)** duoc ghep voi giao tu tai to
hop moi doi -- moi thu khac (tai to hop, chon self/mate tung khoang) giong
het nhau:

- **backcross** (lai nguoc): ban phoi LUON la 1 ca the ban dia THUAN (0 doan
  goc). Vi vay moi doi, xac suat 1 vi tri con giu goc chi con **mot nua** so
  voi doi truoc (mat theo cap so nhan, xap xi `init_fraction x 0.5^doi`) --
  giong pha loang lien tuc bang nuoc sach, hoac khoai lang Nhat lai nguoc
  F1->F2->F3->F4. Ti le to tien **giam don dieu** qua cac doi.
- **isolated** (giao phoi noi bo): ban phoi la **1 replicate KHAC rut ngau
  nhien** tu chinh quan the dang mo phong (ban than no cung mang mot phan
  gene "goc", ti le trung binh xap xi ti le CA quan the). Ve ky vong toan
  hoc, ti le to tien trung binh cua CA quan the **khong doi** qua cac doi --
  chi dao dong do troi dat ngau nhien (kich thuoc mau huu han). Day la diem
  khac biet then chot so voi backcross, va la co so cua phep thu phan biet 2
  che do trong `tests/test_dilution_sim.py`.

**CA HAI che do** deu bi tai to hop **cat nho doan to tien moi doi** (du ti
le tong the co giam hay khong), nen **chieu dai doan trung binh giam don
dieu qua cac doi o ca hai che do** -- day la "dong ho" phan tu cua su kien
lai: doan dai = lai gan day, doan ngan = lai xa xua.

## Bang theo doi

| generation | mean_ancestry_fraction | sd_ancestry_fraction | n_tracts | mean_tract_cm | median_tract_cm | max_tract_cm | fraction_ok | tract_ok | detectable |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.5 | 0 | 480 | 100 | 100 | 100 | 1 | 1 | 1 |
| 1 | 0.497805 | 0.0811448 | 738 | 64.7551 | 72.0231 | 100 | 1 | 1 | 1 |
| 2 | 0.505232 | 0.0850034 | 887 | 54.6812 | 51.6501 | 100 | 1 | 1 | 1 |
| 3 | 0.484422 | 0.0772629 | 992 | 46.8796 | 38.7041 | 100 | 1 | 1 | 1 |
| 4 | 0.474245 | 0.0715449 | 1094 | 41.6156 | 31.9499 | 100 | 1 | 1 | 1 |
| 5 | 0.478984 | 0.0702649 | 1241 | 37.0528 | 28.2869 | 100 | 1 | 1 | 1 |
| 6 | 0.476603 | 0.060606 | 1325 | 34.5312 | 26.2192 | 100 | 1 | 1 | 1 |
| 7 | 0.478177 | 0.0625745 | 1413 | 32.4876 | 26.1419 | 100 | 1 | 1 | 1 |
| 8 | 0.481621 | 0.0518555 | 1548 | 29.8679 | 23.7725 | 100 | 1 | 1 | 1 |
| 9 | 0.48074 | 0.0583104 | 1651 | 27.9534 | 22.1563 | 100 | 1 | 1 | 1 |
| 10 | 0.48073 | 0.051587 | 1771 | 26.0588 | 20.7561 | 100 | 1 | 1 | 1 |
| 11 | 0.480179 | 0.0495361 | 1878 | 24.5459 | 20.0697 | 100 | 1 | 1 | 1 |
| 12 | 0.483105 | 0.0470138 | 2033 | 22.8126 | 17.4471 | 100 | 1 | 1 | 1 |
| 13 | 0.478114 | 0.0451927 | 2160 | 21.2495 | 15.0558 | 100 | 1 | 1 | 1 |
| 14 | 0.480938 | 0.0303705 | 2243 | 20.5841 | 14.349 | 100 | 1 | 1 | 1 |
| 15 | 0.486621 | 0.0387661 | 2310 | 20.2232 | 13.9066 | 100 | 1 | 1 | 1 |
| 16 | 0.487304 | 0.0440771 | 2431 | 19.2436 | 12.9022 | 100 | 1 | 1 | 1 |
| 17 | 0.491056 | 0.0389674 | 2517 | 18.7292 | 12.7135 | 100 | 1 | 1 | 1 |
| 18 | 0.489091 | 0.0361926 | 2637 | 17.8054 | 11.7123 | 100 | 1 | 1 | 1 |
| 19 | 0.481143 | 0.0364989 | 2697 | 17.1263 | 11.439 | 100 | 1 | 1 | 1 |
| 20 | 0.479608 | 0.0322661 | 2753 | 16.7244 | 11.146 | 100 | 1 | 1 | 1 |

