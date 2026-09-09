# Bao cao mo phong pha loang gene 'goc' qua cac doi lai (X7 dilution_sim)

Nguon: `example/backcross_stats.tsv` + `example/backcross_detect.tsv` -- 21 dong du lieu.

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

Du lieu chay tu doi 0 den doi 20. Ti le to tien trung binh di tu 0.5 (doi 0) xuong 0 (doi 20) -- ti so cuoi/dau = 0.
Xu huong nay **giam manh, don dieu** -- phu hop kieu 'lai nguoc' (backcross): moi doi ti le mat xap xi mot nua.
Voi nguong da chon o lenh `detect` (xem CANH BAO BAT BUOC ben tren), su kien lai **con phat hien duoc toi da den doi 5** trong du lieu nay (doi `detectable=1` xa nhat); tu doi sau do trong du lieu nay khong con dong nao `detectable=1`.

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
| 1 | 0.239484 | 0.0582792 | 354 | 64.9449 | 70.0455 | 100 | 1 | 1 | 1 |
| 2 | 0.120245 | 0.0600307 | 232 | 49.7567 | 47.1058 | 100 | 1 | 1 | 1 |
| 3 | 0.0581727 | 0.0378534 | 160 | 34.9036 | 28.4745 | 100 | 1 | 1 | 1 |
| 4 | 0.024352 | 0.0219964 | 84 | 27.8309 | 22.5431 | 100 | 1 | 1 | 1 |
| 5 | 0.011874 | 0.0128268 | 56 | 20.3554 | 15.2765 | 100 | 1 | 1 | 1 |
| 6 | 0.00693704 | 0.00972331 | 35 | 19.0273 | 12.5879 | 100 | 0 | 1 | 0 |
| 7 | 0.00339919 | 0.00811848 | 16 | 20.3952 | 13.3682 | 100 | 0 | 1 | 0 |
| 8 | 0.00166557 | 0.00714484 | 7 | 22.8421 | 9.51839 | 96.4413 | 0 | 1 | 0 |
| 9 | 0.00056182 | 0.00313402 | 5 | 10.7869 | 3.49195 | 23.5566 | 0 | 1 | 0 |
| 10 | 0.000537921 | 0.00313457 | 4 | 12.9101 | 12.8258 | 23.5566 | 0 | 1 | 0 |
| 11 | 3.63744e-05 | 0.000227158 | 1 | 3.49195 | 3.49195 | 3.49195 | 0 | 0 | 0 |
| 12 | 3.63744e-05 | 0.000227158 | 1 | 3.49195 | 3.49195 | 3.49195 | 0 | 0 | 0 |
| 13 | 3.63744e-05 | 0.000227158 | 1 | 3.49195 | 3.49195 | 3.49195 | 0 | 0 | 0 |
| 14 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 15 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 16 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 17 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 18 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 19 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

