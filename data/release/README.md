# data/release — kết quả chính, nén sẵn để tải

| File | Nội dung | Số vùng |
|---|---|---|
| `bulbul_core.bed.gz` | Vùng lõi bảo tồn trên bộ gene chào mào, ngưỡng phủ 0,8 | 248.461 |
| `bulbul_core.cov50.bed.gz` | Cùng loại, ngưỡng phủ 0,5 (rộng hơn, nhiễu hơn) | 411.275 |

Cột: `scaffold  start  end  id_vùng_trên_gà  identity(0–1000)  strand  qlen  coverage  mapq`
Tọa độ theo bộ gene `GCA_013400435.1` (*Pycnonotus jocosus*), hệ 0-based nửa mở như chuẩn BED.

Mã băm của bản chưa nén nằm trong `data/a2/SHA256SUMS.txt`. Tái lập bằng `pipeline/p02_core_map/launch_full.sh`.

**Không kèm trong kho** (quá lớn, tái tạo được): `elements.bed` 130 MB, `elements.fa` 303 MB, `elements.paf` 80 MB, `core_annot*.tsv` 22–36 MB, `unmapped*.txt` 77–85 MB, `anomaly/profile.tsv` 41 MB.
