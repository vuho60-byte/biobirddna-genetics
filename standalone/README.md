# Bản code độc lập một file

File `biobirddna_genetics_standalone.py` là bản tham chiếu thứ hai, không phụ thuộc package và không dùng thư viện ngoài Python.

Chạy:

```bash
python standalone/biobirddna_genetics_standalone.py cross --input examples/input/cross_autosomal_recessive.json
python standalone/biobirddna_genetics_standalone.py infer --input examples/input/infer_z_linked_case.json
python standalone/biobirddna_genetics_standalone.py pigment --input examples/input/pigment_progressive_greying.json
```

Bản một file thuận tiện để model khác đọc toàn bộ logic trong một lượt. Bản modular trong `src/` mới là nền để phát triển production.

## Bản Node.js một file

`biobirddna_genetics_node.mjs` chạy trên Node.js hiện đại, không dependency, hữu ích để tích hợp n8n hoặc worker JavaScript. Phạm vi hiện tại là cross simulation.

```bash
node standalone/biobirddna_genetics_node.mjs cross examples/input/cross_autosomal_recessive.json
```
