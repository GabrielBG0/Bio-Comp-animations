# Bio-Comp-animations

Manim animations for MO640 / MC668 Bioinformatica (Prof. Zanoni Dias), aulas 3 to 6.
Every algorithm explained in the masterclasses is animated step by step, with the
masterclass text shown as subtitles.

## Setup

Uses [3b1b/manim](https://github.com/3b1b/manim) (`manimgl`) via [uv](https://docs.astral.sh/uv/).

System deps (Debian/Ubuntu): `sudo apt install libpango1.0-dev libcairo2-dev pkg-config python3-dev ffmpeg texlive texlive-latex-extra texlive-fonts-extra dvisvgm`

```
uv sync
uv run manimgl scenes.py Hello      # interactive preview
uv run manimgl scenes.py Hello -w   # write video file
```

The `aula*.py` files import `common`/`rearr` from the repo root, which manimgl does not put on `sys.path`, so use `PYTHONPATH=. uv run manimgl aula6.py A6_Hurdles -w` (or `./render_all.sh 6`).


## Rendering the aulas

```
uv run manimgl aula3.py A3_SeedList          # interactive window
uv run manimgl aula3.py A3_SeedList -w       # write videos/A3_SeedList.mp4 (1080p is the default of manimgl; add -l for a fast 480p preview)
./render_all.sh 5                            # every scene of aula5.py at 1080p
./render_all.sh 6 -l                         # every scene of aula6.py, low quality
./render_all.sh                              # everything
```

Videos go to `videos/` (git ignored). Each scene also writes `subtitles/<Scene>.srt`.

## Subtitles

manimgl has no subtitle feature, so `common.py` implements one:

* `NarratedScene.say(text, *animations)` shows `text` as a subtitle, plays the animations, and then waits long enough to read it (about 3 words per second). Long texts are split into several consecutive subtitle chunks.
* The subtitles are burned into the video (a dark box at the bottom of the frame).
* The same text and timing are exported to `subtitles/<Scene>.srt`, in case you prefer to attach them to the video in a player instead.

## Scenes

| File | Scene | What it animates |
| --- | --- | --- |
| aula3.py | A3_Overview | BLAST pipeline, the four steps |
| | A3_SeedList | seeds, neighbors and the score threshold T |
| | A3_FindHits | hit table / automaton scan of the database |
| | A3_ExtendHSP | gapless extension with drop-off |
| | A3_TwoHit | the two-hit method (diagonals and window A) |
| | A3_GappedBLAST | X-drop dynamic programming, pruned cells |
| | A3_Significance | bit score, E-value and P-value |
| aula4.py | A4_Shotgun | fragments, coverage and expected contigs |
| | A4_OverlapGraph | overlap graph and the shortest common superstring |
| | A4_GreedySCS | greedy algorithm on the overlap graph |
| | A4_GreedyVsTopological | greedy failing where a topological sort works |
| | A4_DeBruijn | the de Bruijn graph, k = 3 |
| | A4_KmerAssembly | assembly as an Eulerian walk of k-mers |
| | A4_Consensus | Phred quality and the CAP3 consensus |
| aula5.py | A5_AdditiveBasics | three point and four point conditions |
| | A5_AdditiveInsertion | additive tree by inserting one object at a time |
| | A5_UPGMA | UPGMA clustering and the dendrogram |
| | A5_NeighborJoining | NJ step by step, and NJ against UPGMA |
| | A5_Sandwich | ultrametric sandwich: Kruskal, cut weights, tree |
| | A5_PerfectPhylogeny | the matrix L and the perfect phylogeny tree |
| | A5_Fitch | small parsimony (an extra example, not on the slides) |
| aula6.py | A6_Events | reversal, transposition, signs, mitochondrial data |
| | A6_ReversalBasics | d(pi, sigma) = d(sigma^-1 pi) |
| | A6_SelectionSort | selection sort with reversals, worst case |
| | A6_Pancakes | pancake sorting, bp - 1 bound, strips |
| | A6_Breakpoints | breakpoints, strips, greedy 2-approximation |
| | A6_SignedReversals | oriented pairs and Bergeron's theorem |
| | A6_Hurdles | framed intervals, hurdles, cut and merge |
| | A6_CycleGraph | cycle graph and d = n + 1 - c + h + f |
| | A6_DCJ | adjacency graph and the DCJ distance |
| | A6_Transpositions | transpositions, breakpoints, naive algorithm |
| | A6_RevTransp | reversals and transpositions together |
| | A6_Summary | the table of the whole lecture |

## How the animations are built

The traces are computed, not typed. Scenes call real code for the algorithm they
show (and assert that the result equals the number printed in the masterclass), so
if an example on screen is wrong the render fails instead of lying.

* `common.py`: palette, `NarratedScene` (subtitles), small widgets.
* `rearr.py`: plain Python for aula 6 (breakpoints, strips, oriented pairs, cycle graph, hurdles, DCJ, transpositions). Run `python rearr.py` to run its self tests.

## Numbering and joining

`python number_videos.py` hardlinks the rendered videos to `videos/ordered/NNN_<Scene>.mp4` (lecture order across aula3 to aula6; unrendered scenes leave a gap). Add `--join` to also write `videos/all.mp4` with ffmpeg (no re-encode).
