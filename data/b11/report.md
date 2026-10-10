# B11 — phép thử vị trí câm: Hox và bộ gene công cụ, người so với chào mào

Ngày tạo: 2026-10-10 · công cụ: `pipeline/p08_silent_sites/silent_sites.py report` · đầu vào: `data/b11/silent_sites.tsv`

Chữ ký mọi kết quả: **STRONG_INFERENCE tối đa** — đọc mục cuối trước khi trích dẫn.

## 1. Gene Hox

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| HOXA1 | HOXA1|human|NM_005522.5|NP_005513.2 | lcl|VWYP01010181.1_cds_NXR74635.1_4218 | 311 | 74.3% | 115 | 41.7% | 31.1% | 1.1248 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXA10 | HOXA10|human|NM_018951.4|NP_061824.3 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXA11 | HOXA11|human|NM_005523.6|NP_005514.1 | lcl|VWYP01010181.1_cds_NXR74638.1_4210 | 281 | 88.6% | 106 | 65.1% | 33.2% | 0.46969 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA13 | HOXA13|human|NM_000522.5|NP_000513.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXA2 | HOXA2|human|NM_006735.4|NP_006726.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXA3 | HOXA3|human|NM_001384335.1|NP_001371264.1 | lcl|VWYP01010181.1_cds_NXR74631.1_4216 | 406 | 80.8% | 162 | 35.8% | 23.7% | 1.45329 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXA4 | HOXA4|human|NM_002141.5|NP_002132.3 | lcl|VWYP01012267.1_cds_NXR75481.1_5139 | 56 | 94.6% | 16 | 81.2% | 43.8% | 0.215762 | không xét (few_4fold_sites) |
| HOXA5 | HOXA5|human|NM_019102.4|NP_061975.2 | lcl|VWYP01010181.1_cds_NXR74632.1_4215 | 265 | 84.9% | 94 | 56.4% | 30.8% | 0.653417 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA6 | HOXA6|human|NM_024014.4|NP_076919.1 | lcl|VWYP01010181.1_cds_NXR74640.1_4214 | 231 | 88.7% | 84 | 50.0% | 30.2% | 0.823959 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA7 | HOXA7|human|NM_006896.4|NP_008827.2 | lcl|VWYP01010181.1_cds_NXR74639.1_4213 | 190 | 85.8% | 69 | 58.0% | 30.5% | 0.616394 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA9 | HOXA9|human|NM_152739.4|NP_689952.1 | lcl|VWYP01010181.1_cds_NXR74634.1_4212 | 169 | 88.2% | 57 | 47.4% | 25.5% | 0.907378 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB1 | HOXB1|human|NM_002144.4|NP_002135.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXB13 | HOXB13|human|NM_006361.6|NP_006352.2 | lcl|VWYP01021258.1_cds_NXR79240.1_9312 | 77 | 83.1% | 25 | 52.0% | 44.8% | 0.766238 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB2 | HOXB2|human|NM_002145.4|NP_002136.1 | lcl|VWYP01011832.1_cds_NXR75345.1_4989 | 191 | 74.9% | 65 | 52.3% | 36.1% | 0.75774 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB3 | HOXB3|human|XM_017024560.2|XP_016880049.1 | lcl|VWYP01016679.1_cds_NXR77252.1_7095 | 249 | 77.9% | 91 | 51.6% | 39.2% | 0.776071 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB4 | HOXB4|human|NM_024015.5|NP_076920.1 | lcl|VWYP01011250.1_cds_NXR75089.1_4708 | 52 | 86.5% | 18 | 27.8% | 26.5% | 2.47188 | không xét (few_4fold_sites) |
| HOXB5 | HOXB5|human|NM_002147.4|NP_002138.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXB6 | HOXB6|human|XM_011524727.4|XP_011523029.1 | lcl|VWYP01029682.1_cds_NXR82451.1_12935 | 192 | 84.4% | 68 | 57.4% | 40.0% | 0.630587 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB7 | HOXB7|human|NM_004502.4|NP_004493.3 | lcl|VWYP01030243.1_cds_NXR82765.1_13284 | 85 | 89.4% | 26 | 46.2% | 38.3% | 0.94925 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB8 | HOXB8|human|NM_024016.4|NP_076921.1 | lcl|VWYP01005537.1_cds_NXR72952.1_2296 | 241 | 92.5% | 92 | 60.9% | 37.5% | 0.553199 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB9 | HOXB9|human|NM_024017.5|NP_076922.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC10 | HOXC10|human|NM_017409.4|NP_059105.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC11 | HOXC11|human|NM_014212.4|NP_055027.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC12 | HOXC12|human|NM_173860.3|NP_776272.1 | lcl|VWYP01011900.1_cds_NXR75363.1_5010 | 59 | 100.0% | 21 | 28.6% | 39.5% | 2.28339 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXC13 | HOXC13|human|NM_017410.3|NP_059106.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC4 | HOXC4|human|NM_014620.6|NP_055435.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC5 | HOXC5|human|NM_018953.4|NP_061826.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC6 | HOXC6|human|NM_004503.4|NP_004494.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC8 | HOXC8|human|NM_022658.4|NP_073149.1 | lcl|VWYP01025071.1_cds_NXR80596.1_10846 | 84 | 94.0% | 25 | 32.0% | 24.5% | 1.77868 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXC9 | HOXC9|human|NM_006897.3|NP_008828.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXD1 | HOXD1|human|NM_024501.3|NP_078777.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXD10 | HOXD10|human|NM_002148.4|NP_002139.2 | lcl|VWYP01016378.1_cds_NXR77150.1_6987 | 338 | 92.6% | 121 | 61.2% | 27.0% | 0.547213 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD11 | HOXD11|human|NM_021192.3|NP_067015.2 | lcl|VWYP01016378.1_cds_NXR77148.1_6988 | 231 | 80.5% | 71 | 39.4% | 28.5% | 1.23579 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD12 | HOXD12|human|NM_021193.4|NP_067016.3 | lcl|VWYP01016378.1_cds_NXR77149.1_6989 | 180 | 79.4% | 53 | 30.2% | 24.8% | 2.00326 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD13 | HOXD13|human|NM_000523.4|NP_000514.2 | lcl|VWYP01016378.1_cds_NXR77152.1_6990 | 111 | 96.4% | 39 | 46.2% | 26.3% | 0.94925 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD3 | HOXD3|human|NM_006898.5|NP_008829.3 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXD4 | HOXD4|human|NM_014621.3|NP_055436.2 | lcl|VWYP01016378.1_cds_NXR77154.1_6984 | 203 | 83.7% | 71 | 47.9% | 26.9% | 0.890179 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD8 | HOXD8|human|NM_019558.4|NP_062458.1 | lcl|VWYP01016378.1_cds_NXR77151.1_6985 | 148 | 91.9% | 43 | 74.4% | 25.6% | 0.312871 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD9 | HOXD9|human|NM_014213.4|NP_055028.3 | lcl|VWYP01016378.1_cds_NXR77153.1_6986 | 272 | 75.7% | 102 | 44.1% | 32.5% | 1.02516 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |

Trung vị lớp này: aa giống 85.3%, vị trí 4-fold giống 48.9%, ngẫu nhiên 30.7% (22 gene đủ dữ liệu).

## 2. Bộ gene công cụ phát triển

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| BMP4 | BMP4|human|NM_001347912.1|NP_001334841.1 | lcl|VWYP01008177.1_cds_NXR73813.1_3284 | 385 | 85.5% | 157 | 45.2% | 30.8% | 0.983003 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| FGF8 | FGF8|human|NM_033163.5|NP_149353.1 | lcl|VWYP01021592.1_cds_NXR79274.1_9355 | 180 | 88.3% | 70 | 62.9% | 34.1% | 0.512751 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| GLI3 | GLI3|human|NM_000168.6|NP_000159.3 | lcl|VWYP01016894.1_cds_NXR77324.1_7176 | 202 | 96.0% | 75 | 48.0% | 25.4% | 0.886495 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HAND2 | HAND2|human|NM_021973.3|NP_068808.1 | lcl|VWYP01017394.1_cds_NXR77531.1_7425 | 213 | 96.2% | 99 | 60.6% | 40.1% | 0.558729 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PAX6 | PAX6|human|NM_001368910.2|NP_001355839.1 | lcl|VWYP01012011.1_cds_NXR75403.1_5053 | 380 | 98.4% | 174 | 53.4% | 26.5% | 0.72705 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PITX1 | PITX1|human|XM_054352830.1|XP_054208805.1 | lcl|VWYP01025178.1_cds_NXR80680.1_10940 | 301 | 85.7% | 109 | 50.5% | 39.0% | 0.810322 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| RUNX2 | RUNX2|human|NM_001024630.4|NP_001019801.3 | lcl|VWYP01007754.1_cds_NXR73638.1_3085 | 401 | 93.8% | 204 | 59.3% | 26.2% | 0.586457 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| SHH | SHH|human|NM_000193.4|NP_000184.1 | lcl|VWYP01034182.1_cds_NXR84230.1_14898 | 279 | 90.7% | 128 | 56.2% | 33.8% | 0.656602 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| SOX9 | SOX9|human|NM_000346.4|NP_000337.1 | lcl|VWYP01028539.1_cds_NXR81957.1_12379 | 436 | 88.1% | 171 | 61.4% | 45.9% | 0.542117 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| TBX4 | TBX4|human|XM_011525490.3|XP_011523792.1 | lcl|VWYP01002178.1_cds_NXR71645.1_839 | 507 | 87.0% | 195 | 44.1% | 31.3% | 1.02575 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| TBX5 | TBX5|human|XM_017019912.2|XP_016875401.1 | lcl|VWYP01028264.1_cds_NXR81877.1_12292 | 459 | 88.7% | 178 | 43.3% | 32.5% | 1.05965 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |

Trung vị lớp này: aa giống 88.7%, vị trí 4-fold giống 53.4%, ngẫu nhiên 32.5% (11 gene đủ dữ liệu).

## 3. Gene đối chứng (housekeeping)

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| ACTB | ACTB|human|NM_001101.5|NP_001092.1 | lcl|VWYP01013830.1_cds_NXR76134.1_5848 | 375 | 100.0% | 178 | 59.0% | 32.7% | 0.593594 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| ALDOA | ALDOA|human|NM_001243177.4|NP_001230106.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| EEF1A1 | EEF1A1|human|NM_001402.6|NP_001393.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| ENO1 | ENO1|human|NM_001353346.3|NP_001340275.1 | lcl|VWYP01020671.1_cds_NXR78965.1_9008 | 191 | 91.6% | 79 | 44.3% | 28.8% | 1.01789 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| GAPDH | GAPDH|human|NM_001289745.3|NP_001276674.1 | lcl|VWYP01008447.1_cds_NXR73957.1_3453 | 325 | 91.7% | 150 | 53.3% | 31.3% | 0.730087 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HSPA8 | HSPA8|human|NM_006597.6|NP_006588.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| LDHA | LDHA|human|NM_001165414.2|NP_001158886.1 | lcl|VWYP01032620.1_cds_NXR83570.1_14174 | 324 | 85.8% | 123 | 48.8% | 25.8% | 0.861467 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PGK1 | PGK1|human|NM_000291.4|NP_000282.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| RPLP0 | RPLP0|human|NM_001002.4|NP_000993.1 | lcl|VWYP01007773.1_cds_NXR73702.1_3147 | 317 | 95.6% | 153 | 49.7% | 26.1% | 0.833828 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| TPI1 | TPI1|human|NM_001159287.1|NP_001152759.1 | lcl|VWYP01031335.1_cds_NXR83054.1_13613 | 211 | 91.5% | 92 | 41.3% | 25.0% | 1.14454 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |

Trung vị lớp này: aa giống 91.7%, vị trí 4-fold giống 49.2%, ngẫu nhiên 27.5% (6 gene đủ dữ liệu).

## 4. So Hox với đối chứng

Trung vị giống 4-fold: Hox 48.9% so với đối chứng 49.2%. Nếu hai số này cùng nằm gần mức ngẫu nhiên, Hox không khác gì phần còn lại của bộ gene ở vị trí câm.

## Cách đọc — BẮT BUỘC

Phép thử này đo **một** đại lượng: độ giống ở vị trí câm (base thứ 3 của codon 4-fold, nơi mọi thay đổi đều
không đổi axit amin) giữa người và chào mào, so với **mức ngẫu nhiên** tính từ thành phần base của chính hai
trình tự. Hai kịch bản cho hai kỳ vọng được viết ra **trước khi** nhìn số:

| Kịch bản | Protein | Vị trí câm (4-fold) | JC |
|---|---|---|---|
| Thừa kế từ tổ tiên chung ~320 Ma, chọn lọc giữ protein (H1/H2) | giống cao ở vùng lõi | tiến về mức ngẫu nhiên | lớn hoặc bão hòa |
| Chép/ghép gần đây giữa hai bộ gene | bất kỳ | gần 1,0 | gần 0 |

Giới hạn:
1. **Một cặp loài, một bộ gene mức scaffold.** CDS chào mào là mô hình gene dự đoán; exon sai sẽ hạ `aa_ident`
   và kéo `ident_4fold` về ngẫu nhiên một cách giả tạo. Vì vậy thống kê chính tính trên **cột lõi** (cửa sổ
   21 aa có ≥ 50 % aa giống) và chỉ gene có ≥ 20 vị trí 4-fold lõi mới được xét. Lặp lại với gà/sẻ vằn
   trước khi ghi nhãn cao hơn STRONG_INFERENCE.
2. **Thành phần GC khác nhau** giữa hai loài làm mức ngẫu nhiên lệch khỏi 0,25; cột `chance_4fold_core` đã
   tính từ chính dữ liệu, so với cột đó, không so với 0,25.
3. **Chọn lọc lên codon** (codon usage bias, vùng điều hòa trong CDS) có thể giữ vị trí câm trên mức ngẫu nhiên
   một chút; điều đó **không** tạo ra giá trị gần 1,0 trên cả gene.
4. Gene đối chứng (housekeeping) phải cho cùng bức tranh; nếu Hox khác hẳn đối chứng theo hướng "giống hệt
   vị trí câm", đó mới là tín hiệu đáng xem.
5. Kết quả **không** phân biệt được "thừa kế" với một kịch bản thiết kế cố tình làm ngẫu nhiên vị trí câm; nhưng
   kịch bản đó không còn là "ghép/chép" và phải tự nêu dự đoán riêng (đúng tinh thần ô H3 trong `docs/00`).

Cấm viết "đã chứng minh". Không nâng nhãn.
