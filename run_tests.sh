#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python -B -m unittest discover -s "$ROOT/tests" -v
(cd "$ROOT/npm" && npm test)
