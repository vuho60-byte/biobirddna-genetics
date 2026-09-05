"""Read a local JSON breeder report request; write Vietnamese Markdown."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from biobirddna_genetics.breeder_report import build_breeder_report
from biobirddna_genetics.errors import BioBirdDNAError


def main() -> int:
    parser = argparse.ArgumentParser(description='Xuất báo cáo hồ sơ cặp ghép tiếng Việt; dữ liệu JSON cục bộ.')
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.input.resolve() == args.output.resolve():
            raise ValueError('input_output_must_differ')
        request = json.loads(args.input.read_text(encoding='utf-8-sig'))
        result = build_breeder_report(request)
        args.output.write_text(result['markdown'], encoding='utf-8', newline='\n')
    except (BioBirdDNAError, ValueError, TypeError, OSError) as exc:
        print(json.dumps({'ok': False, 'message_vi': 'Không tạo báo cáo: hồ sơ không hợp lệ hoặc không đọc/ghi được file.',
                          'error': str(exc)}, ensure_ascii=True), file=sys.stderr)
        return 2
    print(json.dumps({'ok': True, 'output': str(args.output), 'input_sha256': result['input_sha256']}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
