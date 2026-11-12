#!/usr/bin/env bash
# Render every scene of one aula (or all aulas) with manimgl.
#
#   ./render_all.sh            # all aulas, 1080p
#   ./render_all.sh 5          # only aula5.py
#   ./render_all.sh 3 -l       # aula3.py, low quality (fast preview)
#
# Extra flags after the aula number are passed to manimgl (-l, -m, --hd, --uhd).
# Videos go to videos/, subtitles (.srt) go to subtitles/.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"  # manimgl does not put the cwd on sys.path

aulas=(3 4 5 6)
if [[ $# -ge 1 && "$1" =~ ^[0-9]+$ ]]; then
  aulas=("$1")
  shift
fi
flags=("$@")
[[ ${#flags[@]} -eq 0 ]] && flags=(--hd)

for a in "${aulas[@]}"; do
  f="aula${a}.py"
  for scene in $(grep -oE "^class A${a}_[A-Za-z0-9_]+" "$f" | awk '{print $2}'); do
    echo ">>> $f $scene"
    uv run manimgl "$f" "$scene" -w "${flags[@]}"
  done
done
