# Bao cao xac suat gene co hien lai qua dot bien diem (X8 reversion_prob)

Nguon: `pipeline/p06_atavism/example/sweep.tsv` -- 36 dong du lieu.

## CANH BAO BAT BUOC -- doc truoc khi dung so trong bao cao nay cho bat ky ket luan nao

1. **"Xuat hien" khac "nhin thay duoc" khac "lan ra quan the".** Cong cu nay
   chi tinh xac suat MOT ban sao nhiem sac the quay lai DUNG base goc trong
   MOT the he. Neu tinh trang do la LAN (recessive), 1 ban quay nguoc o THE
   DI HOP (ket voi 1 ban con mang allele dot bien) se KHONG bieu hien ra
   kieu hinh -- can ca 2 ban cung mang allele goc (dong hop) moi thay duoc
   o 1 ca the, va can dot bien do LAN RA (tang tan so) trong quan the qua
   nhieu the he nua moi co y nghia hon 1 ca the don le. Ba cau hoi "xuat hien
   1 lan" / "nhin thay duoc o 1 ca the" / "lan ra quan the" la BA CAU HOI
   KHAC NHAU -- cong cu nay CHI tra loi cau dau tien.
2. **Mo hinh BO QUA**: chon loc (selection -- gia dinh moi phan tu trung
   tinh, khong loi the/bat loi khi quay nguoc), dot bien MAT DOAN (deletion
   -- neu phan tu bi xoa han thi KHONG con base nao de "quay nguoc" ve, xem
   canh bao 3 duoi day), va CHUYEN DOI GENE (gene conversion -- co che sao
   chep giua cac ban sao/paralog co the thay doi 1 base nhanh hon dot bien
   diem don thuan, khong duoc mo hinh hoa o day).
3. **Neu phan tu bi XOA HAN** (deletion that su, khong con base nao o vi tri
   do) thi xac suat quay nguoc qua CHINH CON DUONG dot bien diem (point
   mutation) nay bang **DUNG SO 0**, khong phai 1 so nho -- vi khong con gi
   de "dot bien tro lai dung base goc". Cong cu nay CHI mo hinh dot bien
   DIEM (thay 1 base bang 1 base khac o 1 vi tri VAN CON TON TAI); no khong
   ap dung cho truong hop mat doan/xoa gene.

## Hang so dung trong bang duoi day la TAM THOI

`mu`/`kappa` trong bang la THAM SO NGUOI DUNG TRUYEN VAO khi chay lenh `sweep`, KHONG PHAI hang so sinh hoc da xac nhan cho chao mao. Gia tri mac dinh cua CLI (`mu=4.6e-09`, `kappa=2`) lay NGUYEN VAN tu vi du minh hoa trong brief X8, cho TOI KHI `research/raw/R14-lai-to-va-hang-so-toc-do-dot-bien.md` (dang cho, xem `research/briefs/R14-lai-to-va-hang-so-toc-do-dot-bien.md`) mang ve so da co nguon.

## Duong tuan tu co luon ngan hon duong dong thoi khi k>=2? (tu chinh du lieu nay)

26/27 dong co k>=2 trong file nay co thoi gian cho TUAN TU ngan hon DONG THOI -- khop du doan cua brief X8 rang duong tuan tu la kich ban thuan loi nhat (con so CHINH) khi k>=2.

## Bang day du theo to hop (k, mu, n_diploid)

| mu | kappa | transversion | k | n_diploid | generation_years | p_specific | p_rev | wait_gen_sequential | wait_years_sequential | wait_gen_concurrent | wait_years_concurrent |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1e-09 | 2 | 0 | 1 | 10000 | 2 | 0.5 | 5e-10 | 2e+09 | 4e+09 | 100000 | 200000 |
| 1e-09 | 2 | 0 | 1 | 1e+06 | 2 | 0.5 | 5e-10 | 2e+09 | 4e+09 | 1000 | 2000 |
| 1e-09 | 2 | 0 | 1 | 1e+08 | 2 | 0.5 | 5e-10 | 2e+09 | 4e+09 | 10 | 20 |
| 4.6e-09 | 2 | 0 | 1 | 10000 | 2 | 0.5 | 2.3e-09 | 4.34783e+08 | 8.69565e+08 | 21739.1 | 43478.3 |
| 4.6e-09 | 2 | 0 | 1 | 1e+06 | 2 | 0.5 | 2.3e-09 | 4.34783e+08 | 8.69565e+08 | 217.391 | 434.783 |
| 4.6e-09 | 2 | 0 | 1 | 1e+08 | 2 | 0.5 | 2.3e-09 | 4.34783e+08 | 8.69565e+08 | 2.17391 | 4.34783 |
| 1e-08 | 2 | 0 | 1 | 10000 | 2 | 0.5 | 5e-09 | 2e+08 | 4e+08 | 10000 | 20000 |
| 1e-08 | 2 | 0 | 1 | 1e+06 | 2 | 0.5 | 5e-09 | 2e+08 | 4e+08 | 100 | 200 |
| 1e-08 | 2 | 0 | 1 | 1e+08 | 2 | 0.5 | 5e-09 | 2e+08 | 4e+08 | 1 | 2 |
| 1e-09 | 2 | 0 | 2 | 10000 | 2 | 0.5 | 5e-10 | 4e+09 | 8e+09 | 2e+14 | 4e+14 |
| 1e-09 | 2 | 0 | 2 | 1e+06 | 2 | 0.5 | 5e-10 | 4e+09 | 8e+09 | 2e+12 | 4e+12 |
| 1e-09 | 2 | 0 | 2 | 1e+08 | 2 | 0.5 | 5e-10 | 4e+09 | 8e+09 | 2e+10 | 4e+10 |
| 4.6e-09 | 2 | 0 | 2 | 10000 | 2 | 0.5 | 2.3e-09 | 8.69565e+08 | 1.73913e+09 | 9.4518e+12 | 1.89036e+13 |
| 4.6e-09 | 2 | 0 | 2 | 1e+06 | 2 | 0.5 | 2.3e-09 | 8.69565e+08 | 1.73913e+09 | 9.4518e+10 | 1.89036e+11 |
| 4.6e-09 | 2 | 0 | 2 | 1e+08 | 2 | 0.5 | 2.3e-09 | 8.69565e+08 | 1.73913e+09 | 9.4518e+08 | 1.89036e+09 |
| 1e-08 | 2 | 0 | 2 | 10000 | 2 | 0.5 | 5e-09 | 4e+08 | 8e+08 | 2e+12 | 4e+12 |
| 1e-08 | 2 | 0 | 2 | 1e+06 | 2 | 0.5 | 5e-09 | 4e+08 | 8e+08 | 2e+10 | 4e+10 |
| 1e-08 | 2 | 0 | 2 | 1e+08 | 2 | 0.5 | 5e-09 | 4e+08 | 8e+08 | 2e+08 | 4e+08 |
| 1e-09 | 2 | 0 | 3 | 10000 | 2 | 0.5 | 5e-10 | 6e+09 | 1.2e+10 | 4e+23 | 8e+23 |
| 1e-09 | 2 | 0 | 3 | 1e+06 | 2 | 0.5 | 5e-10 | 6e+09 | 1.2e+10 | 4e+21 | 8e+21 |
| 1e-09 | 2 | 0 | 3 | 1e+08 | 2 | 0.5 | 5e-10 | 6e+09 | 1.2e+10 | 4e+19 | 8e+19 |
| 4.6e-09 | 2 | 0 | 3 | 10000 | 2 | 0.5 | 2.3e-09 | 1.30435e+09 | 2.6087e+09 | 4.10948e+21 | 8.21895e+21 |
| 4.6e-09 | 2 | 0 | 3 | 1e+06 | 2 | 0.5 | 2.3e-09 | 1.30435e+09 | 2.6087e+09 | 4.10948e+19 | 8.21895e+19 |
| 4.6e-09 | 2 | 0 | 3 | 1e+08 | 2 | 0.5 | 2.3e-09 | 1.30435e+09 | 2.6087e+09 | 4.10948e+17 | 8.21895e+17 |
| 1e-08 | 2 | 0 | 3 | 10000 | 2 | 0.5 | 5e-09 | 6e+08 | 1.2e+09 | 4e+20 | 8e+20 |
| 1e-08 | 2 | 0 | 3 | 1e+06 | 2 | 0.5 | 5e-09 | 6e+08 | 1.2e+09 | 4e+18 | 8e+18 |
| 1e-08 | 2 | 0 | 3 | 1e+08 | 2 | 0.5 | 5e-09 | 6e+08 | 1.2e+09 | 4e+16 | 8e+16 |
| 1e-09 | 2 | 0 | 5 | 10000 | 2 | 0.5 | 5e-10 | 1e+10 | 2e+10 | 1.6e+42 | 3.2e+42 |
| 1e-09 | 2 | 0 | 5 | 1e+06 | 2 | 0.5 | 5e-10 | 1e+10 | 2e+10 | 1.6e+40 | 3.2e+40 |
| 1e-09 | 2 | 0 | 5 | 1e+08 | 2 | 0.5 | 5e-10 | 1e+10 | 2e+10 | 1.6e+38 | 3.2e+38 |
| 4.6e-09 | 2 | 0 | 5 | 10000 | 2 | 0.5 | 2.3e-09 | 2.17391e+09 | 4.34783e+09 | 7.76839e+38 | 1.55368e+39 |
| 4.6e-09 | 2 | 0 | 5 | 1e+06 | 2 | 0.5 | 2.3e-09 | 2.17391e+09 | 4.34783e+09 | 7.76839e+36 | 1.55368e+37 |
| 4.6e-09 | 2 | 0 | 5 | 1e+08 | 2 | 0.5 | 2.3e-09 | 2.17391e+09 | 4.34783e+09 | 7.76839e+34 | 1.55368e+35 |
| 1e-08 | 2 | 0 | 5 | 10000 | 2 | 0.5 | 5e-09 | 1e+09 | 2e+09 | 1.6e+37 | 3.2e+37 |
| 1e-08 | 2 | 0 | 5 | 1e+06 | 2 | 0.5 | 5e-09 | 1e+09 | 2e+09 | 1.6e+35 | 3.2e+35 |
| 1e-08 | 2 | 0 | 5 | 1e+08 | 2 | 0.5 | 5e-09 | 1e+09 | 2e+09 | 1.6e+33 | 3.2e+33 |

