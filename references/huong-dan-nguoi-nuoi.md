# Nhập hồ sơ và đọc báo cáo cặp ghép

Luồng này đọc file JSON trên máy và xuất Markdown tiếng Việt. Không gửi hồ sơ
lên mạng. Chỉ dùng mã chim/mã lứa; không cần tên, điện thoại hoặc địa chỉ người nuôi.

1. Sao chép `examples/input/breeder_report.json` ra thư mục làm việc riêng.
2. Điền `records` theo `assets/schemas/breeding-record.schema.json`: mã cặp,
   mã bố/mẹ, kiểu hình quan sát, cách xác định giới tính và độ tin do nguồn khai báo;
   từng lứa, số trứng theo giai đoạn, từng con kể cả chết hoặc mất theo dõi.
3. Điền `provenance.source` bằng nguồn sổ ghi chép và `synthetic` bằng `false`
   cho hồ sơ thực tế. Đây là nguồn tự khai, chưa được phần mềm xác minh.
4. Nếu có phả hệ, điền `pedigree` theo schema hiện có của `rank-mates`.
   Không tự điền bố mẹ chưa biết. `welfare_blocked_pairing_ids` ghi mã cặp có
   cờ dừng phúc lợi. Danh sách trống không xác nhận sức khỏe an toàn.
5. Chỉ khi đã có giả định nghiên cứu rõ ràng, điền `simulations` với mã cặp
   và `cross_request` theo schema cross hiện có. Allele N/m là ký hiệu giả định;
   không suy ra N/m từ tên “Bông”, “Nâu” hay ảnh. Bỏ trường này nếu chưa biết.

```text
python -B "<skill-root>/scripts/breeder_report.py" --input <absolute-request.json> --output <absolute-report.md>
```

Function Python: `biobirddna_genetics.breeder_report.build_breeder_report(request)`
trả hồ sơ chuẩn hóa, mô phỏng, kiểm tra cặp, provenance/hash và `markdown`.
Đây là adapter Python riêng, không thêm operation vào sáu core operation hoặc MCP.
Chưa hỗ trợ CSV phẳng: JSON giữ nguyên timeline và các giai đoạn trứng/con.

## Những điểm cần đọc đúng

- `status: missing` là chưa ghi nhận; `status: observed, count: 0` là đã thấy
  không có con/trứng ở giai đoạn đó. Không thay missing bằng 0.
- Mỗi probability `female:...`/`male:...` là phần trên **tổng số con**. Ví dụ
  25% mái màu mục tiêu khác 25% trong riêng nhóm mái. Lứa nhỏ có thể lệch tỷ lệ.
- Mô phỏng giữ đúng các giả định nhập vào; số đời con trong sổ không tự fit mô hình.
  Giới tính và kiểu hình nguồn chưa được xác minh; simulation không xác nhận sự
  phù hợp của giả định với cặp thật. Dùng core `infer` riêng cho nghiên cứu đủ dữ liệu.
- F là ước tính theo phả hệ khai báo. Founder chưa biết được coi không liên quan
  chỉ trong phép tính; F = 0 không chứng minh không cận huyết.
- Trực hệ, anh chị em ruột, cờ phúc lợi hoặc F > 6,25% sẽ hiện DỪNG. Không có
  nút bỏ chặn trong adapter. Cặp còn lại vẫn chưa được phê duyệt; cần phúc tra hồ sơ,
  sức khỏe và người phụ trách/thú y duyệt. Không tự xếp hạng “cặp tốt nhất”.
- Báo cáo luôn `evidence_grade: not_assessed`; không cấp G-level từ ảnh, số con
  hoặc xác suất. Toàn bộ `result.limits` được giữ trong phụ lục, cùng hồ sơ chuẩn hóa.

## Sửa lỗi nhập liệu

CLI trả exit code 2 và thông báo JSON khi hồ sơ lỗi; không tạo báo cáo mới cho input
lỗi. Đường dẫn như `pedigree/0/sex` hoặc mã lỗi core kèm ID lứa/trường chỉ vị trí cần
sửa. Các lỗi core giữ nguyên để đối chiếu sổ nguồn. Input và output phải khác file.
Không sửa dữ liệu thật chỉ để làm test hoặc mô hình khớp.

Golden `examples/output/breeder_report.md` dùng dữ liệu tổng hợp: một cặp anh chị em
ruột bị chặn, một cặp thiếu tổ tiên chờ duyệt, một con chết phôi, một con mất dấu và
giai đoạn số 0/missing được giữ riêng. Đây không phải kết quả nghiên cứu trên chào mào.
