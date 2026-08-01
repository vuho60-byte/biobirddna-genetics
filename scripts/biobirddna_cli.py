from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True

SKILL_ROOT = Path(__file__).resolve().parents[1]
SRC = SKILL_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from biobirddna_genetics.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
