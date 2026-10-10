# B11 — phép thử vị trí câm: Hox và bộ gene công cụ, người so với chào mào

Ngày tạo: 2026-10-10 · công cụ: `pipeline/p08_silent_sites/silent_sites.py report` · đầu vào: `data/b11/chicken/silent_sites.tsv`

Chữ ký mọi kết quả: **STRONG_INFERENCE tối đa** — đọc mục cuối trước khi trích dẫn.

## 1. Gene Hox

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| HOXA1 | HOXA1|human|NM_005522.5|NP_005513.2 | lcl|NC_052533.1_cds_XP_003640745.1_10930 | 300 | 73.0% | 106 | 45.3% | 34.0% | 0.980778 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXA10 | HOXA10|human|NM_018951.4|NP_061824.3 | lcl|NC_052533.1_cds_NP_001292574.1_10960 | 322 | 82.6% | 143 | 51.7% | 31.7% | 0.773264 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA11 | HOXA11|human|NM_005523.6|NP_005514.1 | lcl|NC_052533.1_cds_NP_989950.1_10961 | 286 | 87.4% | 110 | 63.6% | 33.6% | 0.497471 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA13 | HOXA13|human|NM_000522.5|NP_000513.2 | lcl|NC_052533.1_cds_NP_989470.2_10964 | 305 | 86.6% | 126 | 49.2% | 40.7% | 0.848155 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXA2 | HOXA2|human|NM_006735.4|NP_006726.1 | lcl|NC_052533.1_cds_NP_990481.1_10931 | 373 | 87.4% | 145 | 68.3% | 28.0% | 0.41242 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA3 | HOXA3|human|NM_001384335.1|NP_001371264.1 | lcl|NC_052533.1_cds_NP_001336689.2_10943 | 410 | 82.4% | 172 | 40.7% | 25.5% | 1.17298 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA4 | HOXA4|human|NM_002141.5|NP_002132.3 | lcl|NC_052533.1_cds_NP_001025517.1_10946 | 189 | 81.5% | 67 | 61.2% | 45.3% | 0.546445 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA5 | HOXA5|human|NM_019102.4|NP_061975.2 | lcl|NC_052533.1_cds_NP_001305348.1_10949 | 265 | 85.7% | 92 | 58.7% | 32.1% | 0.600089 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA6 | HOXA6|human|NM_024014.4|NP_076919.1 | lcl|NC_052533.1_cds_NP_001026158.1_10950 | 228 | 87.3% | 81 | 59.3% | 33.2% | 0.587648 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA7 | HOXA7|human|NM_006896.4|NP_008827.2 | lcl|NC_052533.1_cds_XP_046773755.1_10954 | 204 | 82.8% | 73 | 53.4% | 32.5% | 0.727673 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXA9 | HOXA9|human|NM_152739.4|NP_689952.1 | lcl|NC_052533.1_cds_NP_001292153.1_10958 | 260 | 80.4% | 95 | 56.8% | 31.5% | 0.642524 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB1 | HOXB1|human|NM_002144.4|NP_002135.2 | lcl|NC_052558.1_cds_NP_001074328.2_61187 | 225 | 72.4% | 75 | 50.7% | 32.7% | 0.804221 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB13 | HOXB13|human|NM_006361.6|NP_006352.2 | lcl|NC_052558.1_cds_XP_001234927.1_61135 | 192 | 72.4% | 66 | 42.4% | 32.6% | 1.09472 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB2 | HOXB2|human|NM_002145.4|NP_002136.1 | lcl|NC_052558.1_cds_XP_003642840.1_61186 | 244 | 78.3% | 87 | 41.4% | 27.4% | 1.1411 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB3 | HOXB3|human|XM_017024560.2|XP_016880049.1 | lcl|NC_052558.1_cds_XP_015154916.1_61172 | 333 | 78.1% | 116 | 44.8% | 34.3% | 0.99781 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB4 | HOXB4|human|NM_024015.5|NP_076920.1 | lcl|NC_052558.1_cds_XP_015154923.1_61146 | 183 | 82.0% | 58 | 44.8% | 33.9% | 0.99781 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB5 | HOXB5|human|NM_002147.4|NP_002138.1 | lcl|NC_052558.1_cds_NP_001020526.1_61144 | 253 | 87.7% | 92 | 46.7% | 28.9% | 0.928781 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB6 | HOXB6|human|XM_011524727.4|XP_011523029.1 | lcl|NC_052558.1_cds_XP_024999768.1_61140 | 197 | 84.8% | 67 | 52.2% | 40.9% | 0.759634 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXB7 | HOXB7|human|NM_004502.4|NP_004493.3 | lcl|NC_052558.1_cds_XP_040509241.1_61139 | 181 | 82.3% | 56 | 53.6% | 35.0% | 0.723811 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB8 | HOXB8|human|NM_024016.4|NP_076921.1 | lcl|NC_052558.1_cds_NP_990242.1_61138 | 241 | 91.7% | 85 | 51.8% | 33.3% | 0.772803 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXB9 | HOXB9|human|NM_024017.5|NP_076922.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC10 | HOXC10|human|NM_017409.4|NP_059105.2 | lcl|NC_052565.1_cds_NP_001292189.1_63993 | 323 | 76.2% | 82 | 45.1% | 25.6% | 0.986758 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXC11 | HOXC11|human|NM_014212.4|NP_055027.1 | lcl|NC_052565.1_cds_XP_428250.2_63994 | 278 | 78.4% | 87 | 39.1% | 30.5% | 1.25453 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXC12 | HOXC12|human|NM_173860.3|NP_776272.1 | lcl|NC_052565.1_cds_XP_426957.3_63995 | 259 | 83.4% | 103 | 32.0% | 34.6% | 1.77453 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXC13 | HOXC13|human|NM_017410.3|NP_059106.2 | lcl|NC_052565.1_cds_XP_001235166.4_63996 | 241 | 84.6% | 99 | 49.5% | 38.7% | 0.839266 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXC4 | HOXC4|human|NM_014620.6|NP_055435.2 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC5 | HOXC5|human|NM_018953.4|NP_061826.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC6 | HOXC6|human|NM_004503.4|NP_004494.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| HOXC8 | HOXC8|human|NM_022658.4|NP_073149.1 | lcl|NC_052565.1_cds_NP_990224.1_63991 | 242 | 93.0% | 79 | 49.4% | 26.4% | 0.843191 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXC9 | HOXC9|human|NM_006897.3|NP_008828.1 | lcl|NC_052565.1_cds_NP_001264211.1_63992 | 260 | 83.8% | 97 | 42.3% | 37.0% | 1.10147 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD1 | HOXD1|human|NM_024501.3|NP_078777.1 | lcl|NC_052538.1_cds_NP_001383931.1_33846 | 69 | 71.0% | 19 | 31.6% | 34.9% | 1.82521 | không xét (few_4fold_sites) |
| HOXD10 | HOXD10|human|NM_002148.4|NP_002139.2 | lcl|NC_052538.1_cds_NP_001383552.1_33869 | 338 | 92.3% | 120 | 65.0% | 28.3% | 0.471456 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD11 | HOXD11|human|NM_021192.3|NP_067015.2 | lcl|NC_052538.1_cds_NP_989951.1_33872 | 236 | 80.9% | 75 | 46.7% | 29.3% | 0.931285 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD12 | HOXD12|human|NM_021193.4|NP_067016.3 | lcl|NC_052538.1_cds_NP_990580.1_33873 | 182 | 79.7% | 56 | 39.3% | 28.7% | 1.24367 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD13 | HOXD13|human|NM_000523.4|NP_000514.2 | lcl|NC_052538.1_cds_NP_990765.1_33874 | 290 | 87.9% | 123 | 35.0% | 33.3% | 1.51423 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD3 | HOXD3|human|NM_006898.5|NP_008829.3 | lcl|NC_052538.1_cds_XP_015145074.1_33856 | 407 | 79.9% | 168 | 47.0% | 34.6% | 0.919023 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| HOXD4 | HOXD4|human|NM_014621.3|NP_055436.2 | lcl|NC_052538.1_cds_NP_001012293.1_33859 | 198 | 85.9% | 72 | 52.8% | 26.9% | 0.744939 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD8 | HOXD8|human|NM_019558.4|NP_062458.1 | lcl|NC_052538.1_cds_XP_015144874.1_33862 | 266 | 85.0% | 100 | 64.0% | 31.3% | 0.490445 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HOXD9 | HOXD9|human|NM_014213.4|NP_055028.3 | lcl|NC_052538.1_cds_NP_001383928.1_33863 | 280 | 76.4% | 108 | 48.1% | 36.9% | 0.88168 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |

Trung vị lớp này: aa giống 82.7%, vị trí 4-fold giống 49.3%, ngẫu nhiên 32.6% (34 gene đủ dữ liệu).

## 2. Bộ gene công cụ phát triển

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| BMP4 | BMP4|human|NM_001347912.1|NP_001334841.1 | lcl|NC_052536.1_cds_NP_990568.4_30021 | 385 | 85.2% | 153 | 49.7% | 31.8% | 0.833828 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| FGF8 | FGF8|human|NM_033163.5|NP_149353.1 | lcl|NC_052537.1_cds_NP_001012785.2_31680 | 185 | 89.7% | 75 | 64.0% | 39.2% | 0.490445 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| GLI3 | GLI3|human|NM_000168.6|NP_000159.3 | lcl|NC_052533.1_cds_NP_001258832.2_11902 | 1380 | 84.5% | 535 | 46.4% | 25.0% | 0.942147 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HAND2 | HAND2|human|NM_021973.3|NP_068808.1 | lcl|NC_052535.1_cds_NP_990297.2_22993 | 215 | 96.3% | 103 | 66.0% | 42.7% | 0.452582 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PAX6 | PAX6|human|NM_001368910.2|NP_001355839.1 | lcl|NC_052536.1_cds_XP_025006366.1_26170 | 423 | 98.8% | 197 | 56.3% | 26.8% | 0.654321 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PITX1 | PITX1|human|XM_054352830.1|XP_054208805.1 | lcl|NC_052544.1_cds_NP_001161156.2_46045 | 302 | 85.8% | 110 | 50.0% | 39.4% | 0.823959 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| RUNX2 | RUNX2|human|NM_001024630.4|NP_001019801.3 | lcl|NC_052534.1_cds_XP_015140566.2_20554 | 323 | 95.7% | 156 | 59.0% | 28.4% | 0.593912 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| SHH | SHH|human|NM_000193.4|NP_000184.1 | lcl|NC_052533.1_cds_NP_990152.1_9677 | 413 | 85.5% | 189 | 52.4% | 36.6% | 0.75573 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| SOX9 | SOX9|human|NM_000346.4|NP_000337.1 | lcl|NC_052549.1_cds_NP_989612.1_52258 | 483 | 90.1% | 208 | 55.3% | 41.8% | 0.680041 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| TBX4 | TBX4|human|XM_011525490.3|XP_011523792.1 | lcl|NC_052550.1_cds_NP_001025708.1_53927 | 507 | 89.2% | 200 | 46.5% | 28.4% | 0.937076 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| TBX5 | TBX5|human|XM_017019912.2|XP_016875401.1 | lcl|NC_052546.1_cds_XP_015150317.1_49525 | 494 | 91.3% | 193 | 41.5% | 30.4% | 1.13784 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |

Trung vị lớp này: aa giống 89.7%, vị trí 4-fold giống 52.4%, ngẫu nhiên 31.8% (11 gene đủ dữ liệu).

## 3. Gene đối chứng (housekeeping)

| Gene | CDS người | CDS chào mào | Codon lõi | aa giống (lõi) | Vị trí 4-fold | Giống 4-fold | Ngẫu nhiên | JC | Đọc |
|---|---|---|---|---|---|---|---|---|---|
| ACTB | ACTB|human|NM_001101.5|NP_001092.1 | lcl|NC_052545.1_cds_NP_990849.1_46697 | 375 | 100.0% | 177 | 58.2% | 32.8% | 0.611382 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| ALDOA | ALDOA|human|NM_001243177.4|NP_001230106.1 | lcl|NW_024095963.1_cds_XP_040512781.1_64599 | 192 | 90.1% | 90 | 44.4% | 38.7% | 1.01245 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |
| EEF1A1 | EEF1A1|human|NM_001402.6|NP_001393.1 | lcl|NC_052534.1_cds_NP_001308445.1_19383 | 462 | 99.4% | 223 | 55.2% | 31.4% | 0.683305 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| ENO1 | ENO1|human|NM_001353346.3|NP_001340275.1 | lcl|NC_052552.1_cds_XP_024998255.2_56004 | 423 | 92.2% | 184 | 48.4% | 27.6% | 0.87454 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| GAPDH | GAPDH|human|NM_001289745.3|NP_001276674.1 | lcl|NC_052532.1_cds_NP_989636.1_3943 | 333 | 92.2% | 156 | 52.6% | 29.6% | 0.75073 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| HSPA8 | HSPA8|human|NM_006597.6|NP_006588.1 | lcl|NC_052555.1_cds_NP_990334.2_58333 | 646 | 99.4% | 295 | 47.1% | 28.2% | 0.9158 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| LDHA | LDHA|human|NM_001165414.2|NP_001158886.1 | lcl|NC_052536.1_cds_NP_990615.2_26693 | 326 | 85.6% | 122 | 53.3% | 26.1% | 0.731535 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| PGK1 | PGK1|human|NM_000291.4|NP_000282.1 | — | — | — | — | — | — | — | không có CDS chào mào mang tên gene này (chú giải thiếu, không kết luận) |
| RPLP0 | RPLP0|human|NM_001002.4|NP_000993.1 | lcl|NC_052546.1_cds_NP_990318.1_49232 | 316 | 95.6% | 156 | 52.6% | 28.9% | 0.75073 | trên ngẫu nhiên nhưng xa 1,0 → thừa kế + ràng buộc codon, không phải chép gần đây |
| TPI1 | TPI1|human|NM_001159287.1|NP_001152759.1 | lcl|NC_052532.1_cds_NP_990782.1_4014 | 248 | 89.5% | 102 | 39.2% | 25.4% | 1.24736 | vị trí câm ở mức ngẫu nhiên → khớp thừa kế lâu đời dưới chọn lọc |

Trung vị lớp này: aa giống 92.2%, vị trí 4-fold giống 52.6%, ngẫu nhiên 28.9% (9 gene đủ dữ liệu).

## 4. So Hox với đối chứng

Trung vị giống 4-fold: Hox 49.3% so với đối chứng 52.6%. Nếu hai số này cùng nằm gần mức ngẫu nhiên, Hox không khác gì phần còn lại của bộ gene ở vị trí câm.

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
