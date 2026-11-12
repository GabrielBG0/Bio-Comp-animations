#!/usr/bin/env python3
"""Number the rendered scenes in lecture order and optionally join them.

    python number_videos.py          # hardlink videos/<Scene>.mp4 -> videos/ordered/NNN_<Scene>.mp4
    python number_videos.py --join   # also write videos/all.mp4 (ffmpeg concat, no re-encode)

The number is the scene's position across aula3..aula6 (order of the `class A<n>_...`
definitions), so it stays stable even if some scenes are not rendered yet.
"""
import os, re, subprocess, sys

AULAS = [3, 4, 5, 6]
os.chdir(os.path.dirname(os.path.abspath(__file__)))
scenes = []
for a in AULAS:
    with open(f"aula{a}.py") as f:
        scenes += re.findall(rf"^class (A{a}_\w+)", f.read(), re.M)

out = "videos/ordered"
os.makedirs(out, exist_ok=True)
for old in os.listdir(out):
    os.remove(os.path.join(out, old))
made = []
for i, s in enumerate(scenes, 1):
    src = f"videos/{s}.mp4"
    if not os.path.exists(src):
        print(f"missing  {i:03d} {s}")
        continue
    dst = f"{out}/{i:03d}_{s}.mp4"
    os.link(src, dst)
    made.append(os.path.basename(dst))
    print(f"ok       {i:03d} {s}")

if "--join" in sys.argv and made:
    with open(f"{out}/.list.txt", "w") as f:
        f.writelines(f"file '{m}'\n" for m in made)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0",
                    "-i", f"{out}/.list.txt", "-c", "copy", "videos/all.mp4"], check=True)
    print("wrote videos/all.mp4")
