#!/usr/bin/env sh
set -eu
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <destination-root>" >&2
  exit 2
fi
SKILL_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
DESTINATION_ROOT=$1
TARGET="$DESTINATION_ROOT/biobirddna-genetics"
LOCK="$TARGET.install-lock"
LOCK_HELD=0
STAGE=
mkdir -p "$DESTINATION_ROOT"
if ! mkdir "$LOCK"; then
  echo "Another BioBirdDNA Genetics installation is active: $LOCK" >&2
  exit 3
fi
LOCK_HELD=1
cleanup() {
  if [ -n "$STAGE" ] && [ -d "$STAGE" ]; then
    rm -rf -- "$STAGE"
  fi
  if [ "$LOCK_HELD" -eq 1 ] && [ -d "$LOCK" ]; then
    rmdir -- "$LOCK"
  fi
}
trap cleanup EXIT HUP INT TERM

if [ -e "$TARGET" ]; then
  echo "Target already exists: $TARGET" >&2
  exit 3
fi
STAGE=$(mktemp -d "$DESTINATION_ROOT/.biobirddna-genetics.install.XXXXXX")
cp -R "$SKILL_ROOT"/. "$STAGE"/
find "$STAGE" -type d \( -name __pycache__ -o -name .pytest_cache -o -name .git \) -prune -exec rm -rf -- {} +
find "$STAGE" -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete

if command -v python3 >/dev/null 2>&1; then
  PYTHON=python3
elif command -v python >/dev/null 2>&1; then
  PYTHON=python
else
  echo "Python 3.11+ is required for package verification." >&2
  exit 4
fi

"$PYTHON" -B "$STAGE/scripts/verify_package.py" --root "$STAGE"
"$PYTHON" -B "$STAGE/scripts/atomic_rename.py" "$STAGE" "$TARGET"
STAGE=
"$PYTHON" -B "$TARGET/scripts/verify_package.py" --root "$TARGET"
rmdir -- "$LOCK"
LOCK_HELD=0
trap - EXIT HUP INT TERM
echo "Installed verified BioBirdDNA Genetics skill to $TARGET"
