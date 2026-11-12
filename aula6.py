"""MO640 / MC668, aula 6: genome rearrangements.

Run one scene:   uv run manimgl aula6.py A6_Breakpoints -w
All the maths (breakpoints, oriented pairs, cycle graph, hurdles, DCJ, ...) comes
from rearr.py, so every permutation you see on screen is computed, not typed.
"""
from manimlib import *

from common import *
import rearr as R


# =====================================================================
# helpers
# =====================================================================
def mkrow(vals, signed=False, y=0.0, x=0.0, size=0.56, font=None, ends=True, colors=None, grey_ends=True):
    """A row of cells. `vals` is the full extended list; with ends=False the
    two extremes (0 and n+1) are not drawn but keep their index slots."""
    show = vals if ends else vals[1:-1]
    N = len(vals)
    if font is None:
        font = 20 if max(len(R.fmt(v, signed)) for v in show) <= 2 else 15
    cells = []
    for idx, v in enumerate(show):
        gi = idx if ends else idx + 1
        fill = None
        if grey_ends and (gi == 0 or gi == N - 1):
            fill = GREY
        elif signed and v < 0:
            fill = C_BAD
        elif signed:
            fill = BLUE
        if colors and gi in colors:
            fill = colors[gi]
        cells.append(cell(R.fmt(v, signed), size=size, font=font, fill=fill))
    g = VGroup(*cells).arrange(RIGHT, buff=0.0)
    g.move_to([x, y, 0])
    g.cells, g.off, g.vals, g.signed, g.size, g.font, g.ends = cells, (0 if ends else 1), list(vals), signed, size, font, ends
    return g


def seg(row, lo, hi):
    return VGroup(*row.cells[lo - row.off: hi - row.off + 1])


def box(mob, color=YELLOW, buff=0.05, width=3):
    return SurroundingRectangle(mob, buff=buff).set_stroke(color, width)


def side(text, x, y, size=22, color=WHITE, anchor=LEFT):
    t = T(text, size, color)
    t.move_to([x, y, 0], aligned_edge=anchor) if anchor is not ORIGIN else t.move_to([x, y, 0])
    return t


class RearrScene(NarratedScene):
    """NarratedScene plus the permutation-row toolkit."""

    def _label(self, label, new, label_x, size, color):
        if not label:
            return VGroup()
        lab = T(label, size, color)
        if label_x is None:
            lab.next_to(new, RIGHT, buff=0.3)
        else:
            lab.move_to([label_x, new.get_y(), 0], aligned_edge=LEFT)
        return lab

    def reverse_step(self, prev, new_vals, lo, hi, label, color=YELLOW, dy=0.62, rt=0.9,
                     label_x=None, label_size=21, y=None, x=None, signed=None):
        signed = prev.signed if signed is None else signed
        rect = box(seg(prev, lo, hi), color)
        new = mkrow(new_vals, signed, y=(prev.get_y() - dy) if y is None else y,
                    x=prev.get_x() if x is None else x, size=prev.size, font=prev.font, ends=prev.ends)
        lab = self._label(label, new, label_x, label_size, color)
        self.play(ShowCreation(rect), run_time=0.35)
        self.play(TransformFromCopy(prev, new), FadeIn(lab), run_time=rt)
        self.play(FadeOut(rect), run_time=0.15)
        return new, lab

    def transp_step(self, prev, new_vals, i, j, k, label, color=YELLOW, dy=0.62, rt=0.9,
                    label_x=None, label_size=21):
        r1 = box(seg(prev, i, j - 1), color)
        r2 = box(seg(prev, j, k - 1), ORANGE)
        new = mkrow(new_vals, prev.signed, y=prev.get_y() - dy, x=prev.get_x(), size=prev.size,
                    font=prev.font, ends=prev.ends)
        lab = self._label(label, new, label_x, label_size, color)
        self.play(ShowCreation(r1), ShowCreation(r2), run_time=0.4)
        self.play(TransformFromCopy(prev, new), FadeIn(lab), run_time=rt)
        self.play(FadeOut(r1), FadeOut(r2), run_time=0.15)
        return new, lab

    def wipe(self, title, keep=()):
        """Fade out everything above the subtitle area except the title."""
        ms = [m for m in self.mobjects if m is not title and m not in keep and m.get_y() > -2.6]
        if ms:
            self.play(FadeOut(Group(*ms)))

    def bp_marks(self, row, color=C_BAD):
        marks = VGroup()
        ext = row.vals
        for x in R.breakpoints(ext):
            a, b = row.cells[x - row.off], row.cells[x + 1 - row.off]
            xm = (a.get_right()[0] + b.get_left()[0]) / 2
            marks.add(Line([xm, row.get_y() - 0.4, 0], [xm, row.get_y() + 0.4, 0]).set_stroke(color, 5))
        return marks


# =====================================================================
class A6_Events(RearrScene):
    """Slides 4 to 15: why gene order matters, the events, and the signs."""

    def construct(self):
        title = self.corner_title("Genome rearrangements: the events")
        self.add(title)
        self.say("Until now we compared sequences letter by letter. For whole genomes that measures the wrong thing. The cabbage and turnip mitochondrial genomes have almost identical genes, but the genes appear in a different order. An alignment would say nearly equal and miss the story.")
        self.say("So we ignore the letters and look at conserved blocks. A genome becomes a list of block numbers, and the question is how many large scale events turn one list into the other.")

        base = mkrow([0, 1, 2, 3, 4, 5, 6], True, y=1.6, x=-1.5, ends=False)
        lbl = T("start", 22, GREY_B).next_to(base, LEFT, buff=0.4)
        self.play(FadeIn(base), FadeIn(lbl))
        self.say("The conservative events change order and orientation but not the set of genes. A reversal cuts a segment, turns it around and puts it back, so the order inverts and every sign flips.")
        v = [0, 1, -4, -3, -2, 5, 6, 7]
        assert R.revx([0, 1, 2, 3, 4, 5, 6][:0] + [0, 1, 2, 3, 4, 5, 6], 2, 4, True)[2:5] == [-4, -3, -2]
        r1, l1 = self.reverse_step(base, [0, 1, -4, -3, -2, 5, 6][:7], 2, 4, "reversal", dy=0.9, label_x=2.4)
        self.hold()
        self.say("A transposition swaps two adjacent segments, with no sign change. A reverse transposition does the same, but the moved segment also turns around.",
                 run_time=0.1)
        e = [0, 1, 2, 3, 4, 5, 6]
        t1 = R.transpx(e, 2, 3, 5)
        assert t1 == [0, 1, 3, 4, 2, 5, 6]
        r2, l2 = self.transp_step(base, t1, 2, 3, 5, "transposition", dy=1.8, label_x=2.4)
        rt_vals = [0, 1, 3, 4, -2, 5, 6]
        r3, l3 = self.reverse_step(base, rt_vals, 2, 4, "reverse transposition", color=TEAL, dy=2.7, label_x=2.4)
        self.hold()
        self.say("Fission splits one chromosome into two, fusion joins two into one, and a translocation exchanges the ends of two chromosomes. The non conservative events change the content: insertion, deletion and duplication.")
        self.say("Why do signs appear? A gene sits on one of the two DNA strands. When a segment is reversed, every gene inside it is read in the opposite direction, and we write that as a minus sign. This small detail will decide whether the problem is easy or NP-hard.")
        self.play(FadeOut(VGroup(base, lbl, r1, l1, r2, l2, r3, l3)))

        # the real data: animal mitochondria
        fly = [0, 1, 2, 3, 4, 5, -8, -7, -6, -9, 10, -13, -12, -11, 14]
        chicken = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9, 11, 12, 13, 14]
        human = list(range(15))
        sz = 0.46
        rh = mkrow(human, True, y=1.5, size=sz, font=14, grey_ends=False)
        rf = mkrow(fly, True, y=0.4, size=sz, font=14, grey_ends=False)
        rc = mkrow(chicken, True, y=-0.7, size=sz, font=14, grey_ends=False)
        names = VGroup(T("human", 22, GREY_B), T("fruit fly", 22, GREY_B), T("chicken", 22, GREY_B))
        for n_, r_ in zip(names, (rh, rf, rc)):
            n_.next_to(r_, UP, buff=0.12, aligned_edge=LEFT)
        self.say("Real data: the 15 non tRNA genes of animal mitochondria, numbered 0 to 14. The human order is the reference, and every other animal is a signed permutation. The fruit fly has three reversed pieces, so its reversal distance is 3.",
                 FadeIn(rh), FadeIn(rf), FadeIn(rc), FadeIn(names))
        hl = VGroup(box(seg(rf, 6, 8), C_HI), box(seg(rf, 9, 9), C_HI), box(seg(rf, 11, 13), C_HI))
        self.play(ShowCreation(hl))
        self.hold()
        self.play(FadeOut(hl))
        hl2 = box(seg(rc, 9, 10), C_HI)
        pair = [10, 9]
        rev_pair = [-x for x in pair[::-1]]
        assert rev_pair == [-9, -10]
        self.say("The chicken only swaps two neighbors, with no sign change. It looks like the smallest possible difference, yet it costs 3 reversals. A single reversal of that pair would give minus 9, minus 10, with the wrong signs. This is the first appearance of what will be called a hurdle.",
                 ShowCreation(hl2))
        self.hold()


# =====================================================================
class A6_ReversalBasics(RearrScene):
    """Slides 16 to 31: unsigned reversals, and why distance between two
    permutations equals sorting one permutation."""

    def construct(self):
        title = self.corner_title("Reversals: d(pi, sigma) = d(sigma^-1 pi)")
        self.add(title)
        self.say("A genome with n blocks and no orientation is a permutation of 1 to n. The reversal rho of i and j inverts the piece from position i to position j. Notice that i and j are positions, not values.")
        e = [0, 1, 4, 3, 2, 5, 6]
        r0 = mkrow(e, False, y=1.5, x=-1.0, ends=False)
        self.play(FadeIn(r0))
        r1, l1 = self.reverse_step(r0, R.revx(e, 2, 4), 2, 4, "rho(2,4)", dy=0.8, label_x=2.2)
        assert R.plain(R.revx(e, 2, 4)) == [1, 2, 3, 4, 5]
        self.say("The distance d of pi and sigma is the minimum number of reversals that turn pi into sigma. Sorting by reversals is the special case where sigma is the identity.")
        self.play(FadeOut(VGroup(r0, r1, l1)))

        pi = [5, 3, 1, 6, 7, 4, 2]
        sg = [6, 3, 2, 7, 5, 1, 4]
        ops = [(6, 7), (3, 6), (1, 5), (2, 4)]
        self.set_subtitle("Here is an example with four reversals that take pi to sigma. The two problems look different, but they are the same one.")
        row = mkrow(R.ext_of(pi), False, y=2.0, x=-3.6, size=0.5, ends=False)
        self.play(FadeIn(row))
        cur = R.ext_of(pi)
        rows_l = [row]
        labs = []
        for (i, j) in ops:
            cur = R.revx(cur, i, j)
            nr, lab = self.reverse_step(rows_l[-1], cur, i, j, f"rho({i},{j})", dy=0.6, label_x=-0.55, rt=0.7)
            rows_l.append(nr)
            labs.append(lab)
        assert R.plain(cur) == sg
        self.hold()

        pos = {v: k + 1 for k, v in enumerate(sg)}
        newp = [pos[v] for v in pi]
        assert newp == [5, 2, 6, 1, 4, 7, 3]
        mapping = VGroup(*[T(f"{v}->{pos[v]}", 22, YELLOW) for v in sg]).arrange(RIGHT, buff=0.35).move_to([2.6, 3.0, 0]).scale(0.85)
        mtxt = T("rename each block by its position in sigma", 20, GREY_B).next_to(mapping, DOWN, buff=0.1)
        self.say("The names of the blocks are arbitrary. Rename each block by the position it has in sigma. Sigma becomes the identity, and pi becomes pi prime, which equals sigma inverse times pi.",
                 FadeIn(mapping), FadeIn(mtxt))
        row2 = mkrow(R.ext_of(newp), False, y=2.0, x=3.1, size=0.5, ends=False)
        self.play(TransformFromCopy(rows_l[0], row2))
        cur = R.ext_of(newp)
        prev = row2
        for (i, j) in ops:
            cur = R.revx(cur, i, j)
            nr, _ = self.reverse_step(prev, cur, i, j, "", dy=0.6, rt=0.7, label_x=20)
            prev = nr
        assert R.plain(cur) == list(range(1, 8))
        self.say("Reversals act on positions, so the same four reversals sort pi prime into the identity. Hence the distance between pi and sigma equals the distance of sigma inverse times pi to the identity, and from now on we only study sorting. The slides add that here d is 4, with 8 different optimal sequences.")


# =====================================================================
class A6_SelectionSort(RearrScene):
    """Slides 32 to 46: selection sort with reversals, and how bad it can be."""

    def construct(self):
        title = self.corner_title("Selection sort with reversals")
        self.add(title)
        self.say("Any sorting algorithm based on swaps can be adapted to reversals. Selection sort never does more than n minus 1 swaps, so Watterson and colleagues used its strategy: put element 1 in position 1, then element 2 in position 2, and so on.")
        self.say("For i from 1 to n minus 1, find the position j of element i. If j is not i, apply rho of i and j. The reversal brings i to position i and never disturbs the positions before i.")

        pi = [6, 4, 3, 7, 2, 5, 8, 1]
        cur = R.ext_of(pi)
        row = mkrow(cur, False, y=2.2, x=-1.6, size=0.52, ends=False)
        self.play(FadeIn(row))
        prev = row
        count = 0
        for i in range(1, len(pi)):
            j = cur.index(i)
            if j == i:
                note = T(f"{i} is already in place", 20, C_GOOD).next_to(prev, RIGHT, buff=0.3)
                self.set_subtitle(f"Element {i} is already at position {i}, so nothing happens.")
                self.play(FadeIn(note), run_time=0.4)
                self.wait(0.6)
                self.play(FadeOut(note), run_time=0.2)
                continue
            cur = R.revx(cur, i, j)
            count += 1
            self.set_subtitle(f"Element {i} is at position {j}. Apply rho of {i} and {j} to bring it to position {i}.")
            prev, lab = self.reverse_step(prev, cur, i, j, f"rho({i},{j})  bring {i}", dy=0.62, label_x=1.6)
            for c in prev.cells[:i]:
                c.sq.set_fill(C_GOOD, 0.35)
        assert R.plain(cur) == list(range(1, 9)) and count == 5
        self.say("Five reversals, and for this input 5 is optimal. That was luck.", extra=0.5)
        self.wipe(title)

        n = 7
        worst = [n] + list(range(1, n))
        cur = R.ext_of(worst)
        left = mkrow(cur, False, y=2.2, x=-3.4, size=0.5, ends=False)
        right = mkrow(cur, False, y=2.2, x=3.2, size=0.5, ends=False)
        hdr = VGroup(T("selection sort", 22, C_HI).next_to(left, UP, buff=0.2), T("optimal", 22, C_GOOD).next_to(right, UP, buff=0.2))
        self.say("How bad can it be? Take n followed by 1, 2, up to n minus 1. Selection sort pushes n one step to the right at every iteration and spends n minus 1 reversals.",
                 FadeIn(left), FadeIn(right), FadeIn(hdr))
        prev = left
        cur = R.ext_of(worst)
        steps = 0
        for i in range(1, n):
            j = cur.index(i)
            if j == i:
                continue
            cur = R.revx(cur, i, j)
            steps += 1
            prev, _ = self.reverse_step(prev, cur, i, j, "", dy=0.55, rt=0.5, label_x=20)
        assert steps == n - 1
        cur2 = R.ext_of(worst)
        prev = right
        for (i, j) in ((2, n), (1, n)):
            cur2 = R.revx(cur2, i, j)
            prev, _ = self.reverse_step(prev, cur2, i, j, f"rho({i},{j})", dy=0.55, rt=0.6, label_x=5.55, label_size=18)
        assert R.plain(cur2) == list(range(1, n + 1))
        ratio = MathT(r"\text{ratio} = \frac{n-1}{2}", 34, YELLOW).move_to([0, -1.9, 0])
        self.say("But two reversals are enough: rho of 2 and n gives n down to 1, and rho of 1 and n gives the identity. The ratio is n minus 1 over 2, which grows with n, so this algorithm has no constant approximation factor.",
                 FadeIn(ratio))
        self.say("The lesson: fixing one element at a time is a local view. Good solutions create order in large pieces. We need a measure of how much order is missing, and that will be the breakpoint.")



# =====================================================================
class A6_Pancakes(RearrScene):
    """Slides 47 to 65: pancake sorting by prefix reversals."""

    def stack(self, vals, x=-3.2, ytop=2.3, h=0.42):
        g = VGroup()
        for k, v in enumerate(vals):
            r = Rectangle(width=0.55 + 0.34 * v, height=h)
            r.set_fill(interpolate_color(YELLOW_E, ORANGE, v / 7), 1).set_stroke(WHITE, 2)
            t = Text(str(v), font_size=20, color=BLACK).move_to(r)
            pc = VGroup(r, t)
            pc.move_to([x, ytop - k * (h + 0.04), 0])
            g.add(pc)
        return g

    def construct(self):
        title = self.corner_title("Pancake sorting")
        self.add(title)
        self.say("A waiter carries a stack of pancakes of different diameters. He wants the smallest on top and the largest at the bottom. The only move: insert a spatula somewhere and flip everything above it. The problem was posed in 1975 under the name Harry Dweighter, a pun on harried waiter.")
        self.say("In our language the stack, read from top to bottom, is a permutation, and a flip is a reversal that starts at position 1. So pancake sorting is sorting by prefix reversals: only rho of 1 and i is allowed.")

        stack0 = [1, 6, 3, 5, 2, 4]
        n = 6
        cur = R.ext_of(stack0)
        st = self.stack(stack0)
        log = VGroup()
        self.say("The greedy algorithm places the largest pancakes first. For i from n down to 2: if pancake i is not on top, flip it to the top, then flip the top down to position i. Pancakes below position i are never touched again.",
                 FadeIn(st))
        moves = []
        for i in range(n, 1, -1):
            j = cur.index(i)
            if j == i:
                continue
            if j != 1:
                moves.append((j, f"bring {i} to the top"))
            moves.append((i, f"send {i} to position {i}" if j != 1 or True else ""))
            # simulate to know: moves added below
            if j != 1:
                cur = R.flip(cur, j)
            cur = R.flip(cur, i)
        assert R.plain(cur) == list(range(1, 7))
        assert len(moves) == 7
        cur = R.ext_of(stack0)
        self.set_subtitle("Follow the spatula. Each flip is a reversal of a prefix.")
        for k, (i, txt) in enumerate(moves):
            cur = R.flip(cur, i)
            ytop = 2.3 - (i - 1) * 0.46 - 0.23
            spat = Line([-5.4, ytop, 0], [-1.0, ytop, 0]).set_stroke(YELLOW, 5)
            above = VGroup(*st[:i])
            self.play(ShowCreation(spat), *[pc[0].animate.set_stroke(YELLOW, 3) for pc in above], run_time=0.4)
            new = self.stack(R.plain(cur))
            entry = T(f"rho(1,{i})  {txt}", 20, WHITE).move_to([3.2, 2.3 - 0.4 * k, 0])
            self.play(Transform(st, new), FadeIn(entry), run_time=0.9)
            self.play(FadeOut(spat), run_time=0.15)
            for m in st:
                m[0].set_stroke(WHITE, 2)
        self.say("The greedy algorithm used 7 flips. It makes at most 2 per iteration, and the last iteration needs at most 1, so it never uses more than 2n minus 3 moves. The optimum for this stack is 6.",
                 extra=0.5)
        self.wipe(title)

        # breakpoints for pancakes and the optimal solution
        row0 = mkrow(R.ext_of(stack0), False, y=2.2, x=-1.4, ends=True)
        bpc = R.bp_pancake(R.ext_of(stack0))
        assert bpc == 7
        self.say("The same breakpoint tool works for flips, with one change: the pair of the extra 0 and the first pancake always counts, because every flip cuts there. A flip changes only one pair, so each move removes at most one breakpoint, and the distance is at least bp minus 1. Here bp is 7, so at least 6 moves are needed.",
                 FadeIn(row0))
        opt = [4, 5, 3, 6, 2, 5]
        cur = R.ext_of(stack0)
        prev = row0
        bpl = T("bp = 7", 20, C_HI).next_to(prev, RIGHT, buff=0.3)
        self.add(bpl)
        bps_seen = [7]
        for i in opt:
            cur = R.flip(cur, i)
            prev, lab = self.reverse_step(prev, cur, 1, i, f"rho(1,{i})", dy=0.62, label_x=2.7, rt=0.7)
            bp = R.bp_pancake(cur)
            bps_seen.append(bp)
            b = T(f"bp = {bp}", 20, C_HI).move_to([4.5, prev.get_y(), 0], aligned_edge=LEFT)
            self.add(b)
        assert R.plain(cur) == list(range(1, 7)) and bps_seen == [7, 6, 5, 4, 3, 2, 1], bps_seen
        self.say("This sequence removes exactly one breakpoint at every move, so 6 is optimal, and an optimal solution cannot waste a move. That is why it is unique.")
        self.hold()

        # a better idea: do not separate neighbors
        self.wipe(title)
        a = [0, 3, 5, 4, 2, 1, 6]
        g = mkrow(a, False, y=2.0, x=-3.4, ends=False, size=0.5)
        s = mkrow(a, False, y=2.0, x=3.0, ends=False, size=0.5)
        hdr = VGroup(T("greedy: 6 flips", 22, C_BAD).next_to(g, UP, buff=0.2), T("keep strips: 3 flips", 22, C_GOOD).next_to(s, UP, buff=0.2))
        self.say("Take the stack 3 5 4 2 1. Greedy starts by flipping the 5 to the top, which separates the neighbors 5 and 4, and ends with 6 flips. A better idea is to move whole strips and never cut at an adjacency.",
                 FadeIn(g), FadeIn(s), FadeIn(hdr))
        cur = list(a)
        prev = g
        gmoves = []
        for i in range(5, 1, -1):
            j = cur.index(i)
            if j == i:
                continue
            if j != 1:
                gmoves.append(j)
                cur = R.flip(cur, j)
            gmoves.append(i)
            cur = R.flip(cur, i)
        assert len(gmoves) == 6 and R.plain(cur) == [1, 2, 3, 4, 5]
        cur = list(a)
        for i in gmoves:
            cur = R.flip(cur, i)
            prev, _ = self.reverse_step(prev, cur, 1, i, "", dy=0.5, rt=0.4, label_x=20)
        cur = list(a)
        prev = s
        for i in (3, 2, 5):
            cur = R.flip(cur, i)
            prev, _ = self.reverse_step(prev, cur, 1, i, f"rho(1,{i})", dy=0.5, rt=0.6, label_x=4.9, label_size=18)
        assert R.plain(cur) == [1, 2, 3, 4, 5]
        self.say("Flip the strip 5 4 to the top with 3, join it into the decreasing run 5 4 3 2 1, then flip everything. Three moves. Cutting only at breakpoints gives a 3-approximation, while the plain greedy is a 4-approximation. For comparison, 5n plus 5 over 3 flips are enough for any stack, and computing the minimum is NP-hard.")


# =====================================================================
class A6_Breakpoints(RearrScene):
    """Slides 67 to 98: breakpoints, strips, and the greedy 2-approximation."""

    def construct(self):
        title = self.corner_title("Breakpoints and strips")
        self.add(title)
        e = [0, 4, 3, 2, 5, 6, 7, 1, 8, 9]
        row = mkrow(e, False, y=2.0, x=-0.8)
        self.say("We need a number that says how far a permutation is from sorted. Extend the permutation with a 0 in front and n plus 1 at the end. Now look at each pair of neighbors. If the two values are consecutive, the pair is an adjacency. Otherwise it is a breakpoint.",
                 FadeIn(row))
        marks = self.bp_marks(row)
        b = R.b_count(e)
        assert b == 4
        bl = MathT(rf"b(\pi)={b}", 32, C_BAD).next_to(row, RIGHT, buff=0.5)
        self.say("In 0 4 3 2 5 6 7 1 8 9 the breakpoints are 0 4, 2 5, 7 1 and 1 8, so b is 4. The identity is the only permutation with b equal to 0, so sorting is the same as removing all breakpoints.",
                 ShowCreation(marks), FadeIn(bl))

        seq = [(2, 4), (4, 5), (4, 7), (3, 4), (5, 6)]
        cur = list(e)
        prev = row
        self.set_subtitle("A reversal changes only two pairs, the ones at its two ends. So it removes at most 2 breakpoints, and we need at least b over 2 reversals. Watch what different reversals do to b.")
        bs = [b]
        for (i, j) in seq:
            cur = R.revx(cur, i, j)
            nb = R.b_count(cur)
            delta = nb - bs[-1]
            bs.append(nb)
            prev, lab = self.reverse_step(prev, cur, i, j, f"rho({i},{j})   b = {nb}   ({delta:+d})", dy=0.55, label_x=3.4, label_size=19, rt=0.7)
            mk = self.bp_marks(prev)
            self.add(mk)
        assert bs == [4, 5, 7, 6, 6, 4]
        self.say("The last step turns 2 7 3 6 into 2 3 7 6. It creates two adjacencies at once, which is the best a reversal can do.")
        self.wipe(title)

        # strips
        e2 = [0, 3, 1, 2, 6, 5, 4, 7]
        r2 = mkrow(e2, False, y=2.0, x=-1.0)
        self.say("A strip is a maximal piece in which all neighboring pairs are adjacencies. Breakpoints are the borders between strips. An increasing strip is already in the right order. A decreasing strip still has to be turned around, which is exactly what a reversal does, so decreasing strips are opportunities.",
                 FadeIn(r2))
        for (s, t, d) in R.strips(e2):
            col = C_GOOD if d == "inc" else C_HI
            self.play(*[c.sq.animate.set_fill(col, 0.5) for c in seg(r2, s, t)], run_time=0.5)
        leg = VGroup(T("increasing", 20, C_GOOD), T("decreasing", 20, C_HI)).arrange(RIGHT, buff=0.6).next_to(r2, DOWN, buff=0.5)
        self.add(leg)
        self.say("A single element strip counts as decreasing, except the strips that contain 0 or n plus 1, which are always increasing.")
        self.hold()
        # theorem k, k-1
        k = 4
        i, j = 4, 6
        self.say("Theorem: if k is in a decreasing strip and k minus 1 is in an increasing strip, some reversal removes a breakpoint. Reverse everything from just after k minus 1 up to and including k, so that k lands right next to k minus 1.",
                 run_time=0.1)
        self.hold()
        self.wipe(title)

        # the greedy algorithm
        start = [0, 3, 4, 1, 2, 7, 8, 5, 6, 9]
        steps = R.greedy_bp(start)
        assert [(s["i"], s["j"]) for s in steps] == [(1, 2), (3, 4), (1, 6), (3, 8), (1, 8)]
        row = mkrow(start, False, y=2.3, x=-1.4)
        prev = row
        self.say("Now the greedy algorithm. Let k be the smallest element in the decreasing strips and l the largest. Take the reversal that cuts after k and after k minus 1. If the result has no decreasing strip left, use instead the mirror reversal for l. If there is no decreasing strip at all, reverse the piece between the first two breakpoints, which removes nothing but prepares the next step.",
                 FadeIn(row))
        mk = self.bp_marks(row)
        self.add(mk)
        texts = {
            "none": "no decreasing strip",
            "k": "k = {k}",
            "l": "k would leave no decreasing strip, so l = {l}",
        }
        for s in steps:
            info = s["info"]
            why = texts[s["kind"]].format(**info) if s["kind"] != "none" else texts["none"]
            if s["kind"] == "l" and s["j"] == 8 and s["i"] == 1 and info.get("l") == 8:
                why = "k = 1 and l = 8 give the same reversal"
            self.set_subtitle(f"Next: {why}. The reversal is rho of {s['i']} and {s['j']}, and b changes by {s['db']}.")
            prev, lab = self.reverse_step(prev, s["after"], s["i"], s["j"],
                                          f"rho({s['i']},{s['j']})  b={R.b_count(s['after'])} ({s['db']:+d})",
                                          dy=0.55, label_x=3.3, label_size=19, rt=0.7)
            self.hold()
            m2 = self.bp_marks(prev)
            self.add(m2)
        self.say("Five reversals, exactly b of pi. The lemmas guarantee that the number of reversals A is at most b, because every free reversal was paid in advance by one that removed two breakpoints. Chain the inequalities: b over 2 is at most d, d is at most A, and A is at most b. So A over d is at most 2.",
                 extra=0.3)
        chain = MathT(r"\lceil b/2\rceil\ \le\ d\ \le\ A\ \le\ b\ \Rightarrow\ A/d\le 2", 34, YELLOW).move_to([0, -2.1, 0])
        self.say("This is the standard recipe for approximation proofs: a lower bound for the optimum, an upper bound for the algorithm, and the ratio of the two. The greedy algorithm is a 2-approximation.",
                 FadeIn(chain))
        self.hold()


# =====================================================================
class A6_SignedReversals(RearrScene):
    """Slides 103 to 144: signed reversals, oriented pairs, Bergeron's theorem."""

    def construct(self):
        title = self.corner_title("Signed reversals and oriented pairs")
        self.add(title)
        self.say("Now suppose we know the orientation of every block. A genome is a signed permutation, and a reversal inverts the order of the piece and flips every sign in it. So rho of 3 and 3 is allowed: it just flips one sign. Here is the surprise of the lecture: with more information the problem becomes easier. Unsigned is NP-hard, signed is polynomial.")

        cab = [0, 1, -5, 4, -3, 2, 6]
        row = mkrow(cab, True, y=2.1, x=-1.0, ends=False)
        self.say("Real data from Palmer and Herbon: the cabbage mitochondrial genome, with turnip as the reference. The plan is visible. First make the block minus 5 to minus 2 uniform, then one reversal fixes order and signs together.",
                 FadeIn(row))
        cur = list(cab)
        prev = row
        for (i, j) in ((3, 3), (5, 5), (2, 5)):
            cur = R.revx(cur, i, j, True)
            prev, _ = self.reverse_step(prev, cur, i, j, f"rho({i},{j})", dy=0.62, label_x=2.4)
        assert cur == list(range(7))
        self.say("Three reversals, and the slide says there are 6 optimal sequences.", extra=0.3)
        self.wipe(title)

        # oriented pairs
        pi = [0, 2, -1, 4, -6, -3, 7, -5, 8]
        row = mkrow(pi, True, y=2.2, x=-0.5)
        self.say("Extend with plus 0 and plus n plus 1, which count as positive. An oriented pair is a pair of elements whose absolute values are consecutive and whose signs are different. One reversal can make such a pair an adjacency: the reversal brings one element next to the other and flips its sign, so the signs now agree.",
                 FadeIn(row))
        pairs = R.oriented_pairs(pi)
        assert len(pairs) == 6
        scores = []
        for (p, q) in pairs:
            i, j = R.pair_reversal(pi, p, q)
            res = R.revx(pi, i, j, True)
            sc = R.score(res)
            scores.append(sc)
            pr = VGroup(box(seg(row, p, p), TEAL), box(seg(row, q, q), TEAL))
            rr = box(seg(row, i, j), YELLOW)
            out = mkrow(res, True, y=0.9, x=-0.5)
            txt = T(f"pair ({R.fmt(pi[p], True)},{R.fmt(pi[q], True)}) sum {pi[p] + pi[q]:+d}  rho({i},{j})  score {sc}", 21, WHITE).move_to([0, 0.0, 0])
            self.set_subtitle(f"The pair {R.fmt(pi[p], True)} and {R.fmt(pi[q], True)} gives the reversal rho of {i} and {j}. The score is the number of oriented pairs in the result, here {sc}.")
            self.play(ShowCreation(pr), ShowCreation(rr), run_time=0.4)
            self.play(FadeIn(out), FadeIn(txt), run_time=0.5)
            self.hold()
            self.play(FadeOut(VGroup(pr, rr, out, txt)), run_time=0.3)
        assert scores == [4, 4, 2, 4, 4, 4], scores
        self.say("The rule: if the sum of the pair is plus 1, reverse from the first element up to just before the second. If it is minus 1, reverse from just after the first up to the second. Bergeron proved that a reversal is oriented exactly when it removes breakpoints.")
        self.play(FadeOut(row))

        # Bergeron trace
        prefer = [(1, 2), (4, 7), (5, 7), (2, 2), (3, 5), (3, 6), (6, 6)]
        steps, fin = R.bergeron(pi, prefer=prefer)
        assert len(steps) == 7 and fin == list(range(9))
        row = mkrow(pi, True, y=2.4, x=-2.6, size=0.5)
        self.say("Bergeron's theorem: while there is an oriented pair, apply an oriented reversal of maximum score. Suppose the process stops after k reversals. Then none of those k reversals was wasted. Why maximum score? Because oriented pairs are our supply of good moves, and we want to keep them alive.",
                 FadeIn(row))
        prev = row
        for s in steps:
            p, q = s["pick"]["pair"]
            vals = s["before"]
            ties = sum(1 for c in s["cands"] if c["score"] == s["best"])
            i, j = s["pick"]["rev"]
            self.set_subtitle(f"The best score is {s['best']}, shared by {ties} pairs. We take the pair {R.fmt(vals[p], True)}, {R.fmt(vals[q], True)} and apply rho of {i} and {j}.")
            prev, lab = self.reverse_step(prev, s["pick"]["result"], i, j,
                                          f"({R.fmt(vals[p], True)},{R.fmt(vals[q], True)}) rho({i},{j}) score {s['pick']['score']}",
                                          dy=0.5, label_x=0.1, label_size=17, rt=0.7)
            self.hold()
        self.say("Seven reversals, and here the process ended exactly at the identity, so d of pi is 7. The slide adds that this instance has 8,144 optimal sequences.")
        self.wipe(title)

        # the trap
        t0 = [0, -3, -1, 2, 4]
        good, fg = R.bergeron(t0, prefer=[(2, 3), (1, 3), (1, 1)])
        assert len(good) == 3 and fg == list(range(5))
        c0 = R.bergeron(t0)[0][0]["cands"]
        assert sorted(c["score"] for c in c0) == [0, 2, 2, 2]
        assert R.signed_distance((1, 3, 2)) == 3
        left = mkrow(t0, True, y=2.0, x=-3.3, ends=False)
        right = mkrow(t0, True, y=2.0, x=3.0, ends=False)
        hdr = VGroup(T("maximum score", 22, C_GOOD).next_to(left, UP, buff=0.2), T("a score 0 move", 22, C_BAD).next_to(right, UP, buff=0.2))
        self.say("Why the words maximum score matter. Take minus 3, minus 1, plus 2. Four oriented pairs give four reversals, and one of them has score 0.",
                 FadeIn(left), FadeIn(right), FadeIn(hdr))
        prev = left
        for s in good:
            i, j = s["pick"]["rev"]
            prev, _ = self.reverse_step(prev, s["pick"]["result"], i, j, f"rho({i},{j}) score {s['pick']['score']}", dy=0.6, label_x=-1.6, label_size=17, rt=0.7)
        cur = R.revx(t0, 1, 2, True)
        assert cur == [0, 1, 3, 2, 4]
        prev = right
        prev, _ = self.reverse_step(prev, cur, 1, 2, "rho(1,2) score 0", dy=0.6, label_x=4.7, label_size=15, color=C_BAD, rt=0.7)
        self.set_subtitle("With the score 0 move we reach plus 1, plus 3, plus 2: all positive, no oriented pair, and it still needs 3 reversals.")
        for (i, j) in ((2, 2), (2, 3), (2, 2)):
            cur = R.revx(cur, i, j, True)
            prev, _ = self.reverse_step(prev, cur, i, j, f"rho({i},{j})", dy=0.6, label_x=4.7, label_size=15, rt=0.6)
        assert cur == list(range(5))
        self.say("To swap 3 and 2, one reversal must contain both and flip both signs, and fixing the signs afterwards needs one reversal for each. Total 4 against 3 with maximum score. This is the same situation as the chicken in the real data.")


# =====================================================================
class A6_Hurdles(RearrScene):
    """Slides 145 to 203: framed intervals, hurdles, cutting and merging."""

    def construct(self):
        title = self.corner_title("Framed intervals and hurdles")
        self.add(title)
        self.say("Bergeron's method stops when all elements are positive. If the permutation is the identity we are done. If not, we meet something like plus 1, plus 3, plus 2: no oriented pair, and no reversal that removes a breakpoint. We need to describe these obstacles and learn the cheapest way through them.")
        e = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8, 10, 12, 11, 13]
        N = len(e)
        row = mkrow(e, False, y=1.9, x=-0.6, size=0.46, font=17)
        self.say("Take a reduced permutation with only positive elements, extended with 0 and n plus 1, and think of it as circular. A framed interval is a piece from i to i plus k, with k at least 2, that contains all the integers in between. It starts with its smallest value, ends with its largest, and is scrambled inside. A hurdle is a framed interval that contains no other framed interval.",
                 FadeIn(row))
        for p in (0, 3):
            i = e[p]
            c = m = 0
            rect = info = None
            hdr = T(f"start at position {p}, value {i}", 22, YELLOW).move_to([4.7, 2.3, 0])
            self.add(hdr)
            if p == 0:
                self.set_subtitle("To find them, fix a start and rename each value v as v minus the start value. Walk right keeping c, the number of elements visited, and m, the largest renamed value. The piece is a framed interval when c is at least 2 and the current renamed value equals both m and c.")
            for step in range(1, N):
                x = (p + step) % N
                r = (e[x] - i) % N
                c += 1
                m = max(m, r)
                ok = c >= 2 and r == m and r == c
                new_rect = box(VGroup(*row.cells[p:x + 1]), C_GOOD if ok else YELLOW)
                new_info = T(f"c={c} m={m} r={r}" + (" framed!" if ok else ""), 22, C_GOOD if ok else WHITE).move_to([4.7, 1.6, 0])
                if rect is None:
                    rect, info = new_rect, new_info
                    self.play(ShowCreation(rect), FadeIn(info), run_time=0.35)
                else:
                    self.play(Transform(rect, new_rect), Transform(info, new_info), run_time=0.35)
                if ok:
                    break
            self.hold()
            self.play(FadeOut(rect), FadeOut(info), FadeOut(hdr))

        hs, cands = R.hurdles(e)
        assert sorted(hs) == [(0, 4), (3, 5), (7, 4), (10, 4)]
        br = VGroup()
        lvl = 0
        for (p, L) in sorted(cands, key=lambda t: t[1]):
            if p + L > N:
                continue
            a, b = row.cells[p], row.cells[p + L - 1]
            y = row.get_y() - 0.45 - 0.16 * lvl
            br.add(Line([a.get_left()[0] + 0.03, y, 0], [b.get_right()[0] - 0.03, y, 0]).set_stroke(TEAL, 3))
            lvl += 1
        self.say("Doing this for every start gives the shortest framed interval from each position. There are also bigger framed intervals, made of consecutive blocks, and circular ones that pass through 13 and 0.",
                 ShowCreation(br, lag_ratio=0.1), run_time=2.0)
        keep = VGroup(*[box(VGroup(*row.cells[p:p + L]), [C_GOOD, C_HI, C_QUERY, YELLOW][k], buff=0.04, width=4) for k, (p, L) in enumerate(sorted(hs))])
        self.say("A candidate is a hurdle exactly when it contains no other candidate. Here the hurdles are the four blocks 0 2 1 3, then 3 6 5 4 7, then 7 9 8 10, then 10 12 11 13. All the other framed intervals contain one of them.",
                 FadeOut(br), ShowCreation(keep))
        self.play(FadeOut(VGroup(row, keep)))

        # cutting and merging
        c1 = [0, 2, 4, 3, 1, 5]
        r1 = mkrow(c1, True, y=2.0, x=-2.0, ends=False)
        cut = R.revx(c1, 1, 3, True)
        assert cut == [0, -3, -4, -2, 1, 5]
        self.say("Why is a hurdle a cost? Inside it everything is positive, so no reversal removes a breakpoint. The first reversal applied to it is wasted: it only creates negative elements so that the oriented method works again. To cut a hurdle, reverse the elements strictly between i and i plus 1.",
                 FadeIn(r1))
        r1b, _ = self.reverse_step(r1, cut, 1, 3, "cut", dy=0.7, label_x=1.5, color=C_HI)
        m1 = [0, 2, 5, 4, 3, 6, 1, 7]
        r2 = mkrow(m1, True, y=-0.1, x=-2.0, ends=False)
        mer = R.revx(m1, 5, 5, True)
        assert mer == [0, 2, 5, 4, 3, -6, 1, 7]
        self.say("To merge two hurdles, reverse everything from the end of the first to the beginning of the second. Here hurdle 1 is 2 5 4 3 6 and hurdle 2 is the circular piece 6 1 7 0 2. They share the 6, so the merge reverses just that element, and one reversal disturbs two hurdles at once.",
                 FadeIn(r2))
        r2b, _ = self.reverse_step(r2, mer, 5, 5, "merge", dy=0.7, label_x=2.0, color=C_GOOD)
        self.say("A hurdle is simple if cutting it and applying oriented reversals of maximum score leaves fewer hurdles. Otherwise it is a super hurdle: you cut it and another hurdle appears in its place. The moral: never cut a super hurdle when you can merge, because a merge handles two hurdles with one wasted reversal.")
        self.wipe(title)

        # exercise 4: full algorithm on (+2 +1 +3 +5 +4)
        ex = [0, 2, 1, 3, 5, 4, 6]
        hs2, _ = R.hurdles(ex)
        assert sorted(hs2) == [(0, 4), (3, 4)]
        row = mkrow(ex, True, y=2.3, x=-2.6, size=0.5)
        self.say("Now the complete algorithm on plus 2, plus 1, plus 3, plus 5, plus 4. All elements are positive, so there is no oriented pair. It has two hurdles, 0 2 1 3 and 3 5 4 6. Two is even, so we merge them. Hurdle 1 ends at 3 and hurdle 2 begins at 3, so the merge reverses the single element at position 3.",
                 FadeIn(row))
        hb = VGroup(box(seg(row, 0, 3), C_GOOD), box(seg(row, 3, 6), C_HI))
        self.play(ShowCreation(hb))
        cur = R.revx(ex, 3, 3, True)
        prev, _ = self.reverse_step(row, cur, 3, 3, "merge rho(3,3)", dy=0.55, label_x=0.0, label_size=17, color=C_GOOD)
        self.play(FadeOut(hb))
        steps, fin = R.bergeron(cur, prefer=[(2, 3), (1, 3), (2, 4), (3, 5), (2, 5)])
        assert fin == list(range(7)) and len(steps) == 5
        self.set_subtitle("From here on, oriented reversals of maximum score.")
        for s in steps:
            i, j = s["pick"]["rev"]
            prev, _ = self.reverse_step(prev, s["pick"]["result"], i, j, f"rho({i},{j}) score {s['pick']['score']}", dy=0.55, label_x=0.0, label_size=17, rt=0.7)
        self.say("Six reversals in total, and an exhaustive search confirms that d is 6. At the third step the pair plus 1, minus 2 had score 0, and choosing it would have produced a new hurdle. The maximum score rule avoids it.")

        self.wipe(title)
        rules = [
            ("oriented pair exists", "oriented reversal of maximum score"),
            ("2 hurdles", "merge the two"),
            ("even number, more than 2", "merge two non-consecutive hurdles"),
            ("1 hurdle", "cut it"),
            ("odd number, one is simple", "cut a simple hurdle"),
            ("3 hurdles, all super", "merge two"),
            ("odd above 3, all super", "merge two non-consecutive"),
        ]
        tbl = VGroup()
        for k, (a, b) in enumerate(rules):
            ta = T(a, 21, YELLOW).move_to([-6.2, 2.4 - 0.55 * k, 0], aligned_edge=LEFT)
            tb = T(b, 21, WHITE).move_to([-1.0, 2.4 - 0.55 * k, 0], aligned_edge=LEFT)
            tbl.add(ta, tb)
        self.say("The whole algorithm is a table of cases. If there is an oriented pair, use it. Otherwise merge hurdles in pairs, choosing non-consecutive ones so that what lies between them does not turn into a new hurdle. If one hurdle is left over, cut a simple one. If all are super, even the best strategy loses one reversal, and that unlucky case is called a fortress.",
                 FadeIn(tbl), run_time=1.5)
        cost = MathT(r"\text{time } O(n^3):\ d(\pi)\le n+1\ \text{rounds, each } O(n^2)", 28, GREY_B).move_to([0, -1.8, 0])
        self.say("Each round costs O of n squared, to find the best oriented reversal or all the hurdles, and there are at most n plus 1 rounds, so the total is order n cubed.", FadeIn(cost))


# =====================================================================
class A6_CycleGraph(RearrScene):
    """Slides 204 to 214: the cycle graph and the exact distance formula."""

    DX = 0.88

    def draw_graph(self, e, y0=-0.1):
        blocks, black, gray = R.cg_edges(e)
        flat = [v for blk in blocks for v in blk]
        N = len(flat)
        x0 = -(N - 1) * self.DX / 2
        xs = {v: x0 + k * self.DX for k, v in enumerate(flat)}
        dots = VGroup()
        for v in flat:
            d = Dot([xs[v], y0, 0], radius=0.07).set_color(WHITE)
            lb = T(str(v), 17, GREY_B).move_to([xs[v], y0 - 0.3, 0])
            dots.add(d, lb)
        genes = VGroup()
        for blk, g in zip(blocks, e):
            xm = (xs[blk[0]] + xs[blk[-1]]) / 2
            txt = R.fmt(g, True) if g not in (0, len(e) - 1) else ("+0" if g == 0 else f"+{g}")
            genes.add(T(txt, 20, YELLOW).move_to([xm, y0 - 0.75, 0]))
        bmob, gmob = {}, {}
        bl = VGroup()
        for (a, b) in black:
            ln = Line([xs[a], y0, 0], [xs[b], y0, 0]).set_stroke(WHITE, 7)
            bmob[frozenset((a, b))] = ln
            bl.add(ln)
        gl = VGroup()
        for (a, b) in gray:
            arc = ArcBetweenPoints(np.array([xs[a], y0 + 0.05, 0]), np.array([xs[b], y0 + 0.05, 0]), angle=-1.6 if xs[a] < xs[b] else 1.6)
            arc.set_stroke(GREY_B, 3)
            gmob[frozenset((a, b))] = arc
            gl.add(arc)
        return dict(dots=dots, genes=genes, black=bl, gray=gl, bmob=bmob, gmob=gmob, all=VGroup(dots, genes, bl, gl))

    def color_cycles(self, g, e):
        cols = [C_HI, C_QUERY, C_GOOD, PURPLE_C, YELLOW, TEAL]
        cyc = R.cg_cycles(e)
        anims_all = []
        for k, cy in enumerate(cyc):
            col = cols[k % len(cols)]
            anims = []
            for (v, w) in cy:
                anims.append(g["bmob"][frozenset((v, w))].animate.set_stroke(col, 7))
                partner = w + 1 if w % 2 == 0 else w - 1
                anims.append(g["gmob"][frozenset((w, partner))].animate.set_stroke(col, 4))
            anims_all.append(anims)
        return cyc, anims_all

    def construct(self):
        title = self.corner_title("The cycle graph")
        self.add(title)
        e1 = [0, 1, -2, -5, -4, -6, 3, 7]
        self.say("Breakpoints look at one genome at a time. The cycle graph puts the neighborhoods of both genomes, pi and the identity, into one structure. Each gene x becomes two vertices, its extremities 2x minus 1 and 2x, in this order if x is positive and in the reverse order if x is negative.")
        g1 = self.draw_graph(e1)
        self.play(FadeIn(g1["dots"]), FadeIn(g1["genes"]))
        self.say("Black edges join the extremities that touch in pi: one for each pair of neighbors. Gray edges join the extremities that touch in the identity: 0 with 1, 2 with 3, 4 with 5, and so on.",
                 ShowCreation(g1["black"]), ShowCreation(g1["gray"]), run_time=1.5)
        cyc, anim_groups = self.color_cycles(g1, e1)
        assert len(cyc) == 4
        self.set_subtitle("Every vertex has exactly one black edge and one gray edge. So if you leave a vertex by its black edge, continue by the gray edge of the vertex you reach, and keep alternating, you must come back to the start. The graph splits uniquely into alternating cycles.")
        cnt = MathT(r"c(\pi)=0", 34, YELLOW).move_to([4.8, 2.4, 0])
        self.add(cnt)
        for k, an in enumerate(anim_groups):
            new = MathT(rf"c(\pi)={k + 1}", 34, YELLOW).move_to([4.8, 2.4, 0])
            self.play(*an, Transform(cnt, new), run_time=0.9)
            self.wait(0.3)
        self.hold()
        bound = MathT(r"d(\pi)\ \ge\ n+1-c(\pi)=7-4=3", 32, YELLOW).move_to([0, -2.0, 0])
        self.say("The identity has n plus 1 cycles, the maximum, because every black edge sits on top of a gray one. A reversal cuts two black edges and reconnects them, so it changes the number of cycles by at most one. Hence d is at least n plus 1 minus c, which is 3 here.",
                 FadeIn(bound))

        # flip the sign of gene 2
        e2 = [0, 1, 2, -5, -4, -6, 3, 7]
        g2 = self.draw_graph(e2)
        self.play(FadeOut(g1["all"]), FadeOut(cnt), FadeOut(bound))
        self.play(FadeIn(g2["dots"]), FadeIn(g2["genes"]), FadeIn(g2["black"]), FadeIn(g2["gray"]))
        cyc2, anim2 = self.color_cycles(g2, e2)
        assert len(cyc2) == 5
        self.say("Flip the sign of the minus 2. Only the vertices of gene 2 change order, and the cycle that had three black edges splits into one with a single black edge and one with two.",
                 *[a for an in anim2 for a in an], run_time=1.5)
        c2 = MathT(r"c=5,\quad d\ge 7-5=2", 34, YELLOW).move_to([0, -2.0, 0])
        seq = [0, 1, 2, -5, -4, -6, 3, 7]
        a1 = R.revx(seq, 5, 6, True)
        a2 = R.revx(a1, 3, 5, True)
        assert a2 == list(range(8))
        self.say("Now two reversals sort it: rho of 5 and 6, then rho of 3 and 5. So d is 2. This is the unsigned problem in miniature: for the unsigned permutation 1 2 5 4 6 3, the sign of the 2 decides between distance 3 and distance 2, and the best choice of signs is the NP-hard part.",
                 FadeIn(c2))
        self.play(FadeOut(g2["all"]), FadeOut(c2))

        # exact formula
        form = MathT(r"d(\pi)=n+1-c(\pi)+h(\pi)+f(\pi)", 44, YELLOW).move_to([0, 2.0, 0])
        self.say("Hannenhalli and Pevzner proved the exact formula for every signed permutation. The term n plus 1 minus c is the number of cycles still to create. The term h is the one wasted reversal per hurdle. The term f is 1 for a fortress, an odd number of hurdles that are all super, and 0 otherwise.",
                 FadeIn(form))
        rows = []
        for name, p in (("(+3,+2,+1)", [3, 2, 1]), ("(+2,+1,+3,+5,+4)", [2, 1, 3, 5, 4]), ("(+2,+1,+4,+3)", [2, 1, 4, 3])):
            ex = R.ext_of(p)
            c = R.c_count(ex)
            h = len(R.hurdles(ex)[0])
            n = len(p)
            rows.append((name, n + 1, c, h, 0, n + 1 - c + h))
        assert [r[5] for r in rows] == [3, 6, 5], rows
        hdr = T(f"{'permutation':<22}{'n+1':>4}{'c':>5}{'h':>5}{'f':>5}{'d':>5}", 20, GREY_B).move_to([0, 0.8, 0])
        lines = VGroup(hdr)
        for k, r in enumerate(rows):
            lines.add(T(f"{r[0]:<22}{r[1]:>4}{r[2]:>5}{r[3]:>5}{r[4]:>5}{r[5]:>5}", 20, WHITE).move_to([0, 0.2 - 0.5 * k, 0]))
        self.say("Test it. Plus 3, plus 2, plus 1 has two cycles and one hurdle, so d is 3. The two permutations of the hurdle exercise give 6 and 5, exactly what the algorithm produced.",
                 FadeIn(lines))
        self.say("How often do hurdles matter? Almost never. A random signed permutation has a hurdle with probability of order n to the minus 2, and a fortress with probability of order n to the minus 15. For real data, the distance is in practice n plus 1 minus c.")


# =====================================================================
class A6_DCJ(RearrScene):
    """Slides 215 to 220: several chromosomes and the DCJ model."""

    def construct(self):
        title = self.corner_title("Double cut and join")
        self.add(title)
        self.say("Real genomes have several chromosomes, so reversals are not the only events. Hannenhalli and Pevzner gave a polynomial algorithm for reversals, translocations, fusions and fissions, in the paper called Transforming men into mice. For human and mouse, the scenario found has 245 events: 149 reversals, 93 translocations and 3 fissions.")
        self.say("It also needs at least 190 reuses of breakpoints. If chromosomes broke at random places, hitting the same region twice would be rare, so some regions of the genome must be fragile and break again and again.")
        base = [0, 1, 2, 3, 4, 5]
        row = mkrow(base, True, y=1.7, x=-1.0, ends=False)
        gap = VGroup(*[Line([0, -0.4, 0], [0, 0.4, 0]).set_stroke(C_BAD, 5) for _ in range(2)])
        self.say("Yancopoulos, Attie and Friedberg proposed a single operation that covers everything: DCJ, double cut and join. Cut the genome at two points, which leaves four loose ends, and join the ends again in a different way.",
                 FadeIn(row))
        a, b = row.cells[0], row.cells[1]
        xa = (row.cells[0].get_right()[0] + row.cells[1].get_left()[0]) / 2
        xb = (row.cells[2].get_right()[0] + row.cells[3].get_left()[0]) / 2
        gap[0].move_to([xa, row.get_y(), 0])
        gap[1].move_to([xb, row.get_y(), 0])
        self.play(ShowCreation(gap))
        new = R.revx(base, 2, 3, True)
        assert new == [0, 1, -3, -2, 4, 5]
        out = mkrow(new, True, y=0.6, x=-1.0, ends=False)
        self.say("Cutting before gene 2 and after gene 3, and rejoining the ends the other way around, is a reversal. With the other genome arrangements, the same operation gives a translocation, a fusion or a fission. An exchange of blocks, such as a transposition, takes two DCJs.",
                 TransformFromCopy(row, out))
        self.play(FadeOut(VGroup(row, gap, out)))

        pi = [[1, 4], [-2, -3, -5]]
        sg = [[1, 2, 3, 4, 5]]
        va, vb, edges, comps = R.dcj_adjacency_graph(pi, sg)
        d, c, i = R.dcj_distance(pi, sg)
        assert (d, c, i) == (3, 1, 2)

        def vtxt(v):
            return " ".join(f"{g}{k}" for (g, k) in sorted(v, key=lambda t: (t[0], 0 if t[1] == 't' else 1)))
        self.say("The adjacency graph. Each gene has a tail t and a head h. A positive gene reads tail then head, a negative gene reads head then tail. Two extremities that touch inside a chromosome form an adjacency, and an extremity at the end of a chromosome is a telomere. Here pi has two chromosomes, plus 1 plus 4 and minus 2 minus 3 minus 5, and sigma has one.")
        top_y, bot_y = 1.5, -1.0
        px = [-5.7 + 1.9 * k for k in range(len(va))]
        sx = [-4.75 + 1.9 * k for k in range(len(vb))]

        def vbox(txt, x, y):
            r = Rectangle(width=1.5, height=0.5).set_stroke(WHITE, 2).move_to([x, y, 0])
            t = T(txt, 18, WHITE).move_to(r)
            return VGroup(r, t)
        pa = VGroup(*[vbox(vtxt(v), px[k], top_y) for k, v in enumerate(va)])
        pb = VGroup(*[vbox(vtxt(v), sx[k], bot_y) for k, v in enumerate(vb)])
        lp = T("pi", 24, C_HI).move_to([-6.7, top_y + 0.0, 0])
        ls = T("sigma", 24, C_QUERY).move_to([-6.7, bot_y, 0])
        self.play(FadeIn(pa), FadeIn(pb), FadeIn(lp), FadeIn(ls))

        def off(v, e):
            lst = sorted(v, key=lambda t: (t[0], 0 if t[1] == 't' else 1))
            k = lst.index(e)
            return (k - (len(lst) - 1) / 2) * 0.55
        elines = []
        for (x, y, ext) in edges:
            ln = Line([px[x] + off(va[x], ext), top_y - 0.25, 0], [sx[y] + off(vb[y], ext), bot_y + 0.25, 0]).set_stroke(GREY_B, 3)
            elines.append((x, y, ln))
        self.play(*[ShowCreation(l) for (_, _, l) in elines], run_time=1.5)
        cols = [C_GOOD, C_QUERY, C_HI, YELLOW, PURPLE_C]
        anims = []
        for k, comp in enumerate(comps):
            col = cols[k % len(cols)]
            members = set(comp["vertices"])
            for (x, y, ln) in elines:
                if ("a", x) in members:
                    anims.append(ln.animate.set_stroke(col, 5))
        legend = VGroup()
        for k, comp in enumerate(comps):
            legend.add(T(f"{comp['kind']}, {comp['edges']} edges", 20, cols[k % len(cols)]))
        legend.arrange(RIGHT, buff=0.5).move_to([0, -2.1, 0]).scale(0.9)
        self.say("Every vertex has degree 1 or 2, so the graph is a collection of cycles and paths. A path is odd if it has an odd number of edges. Here there is one cycle, two odd paths and one even path.",
                 *anims, FadeIn(legend), run_time=1.5)
        form = MathT(r"d_{DCJ}=n-\left(c+\frac{i}{2}\right)=5-\left(1+\frac{2}{2}\right)=3", 34, YELLOW).move_to([1.5, 2.75, 0])
        self.say("The reward for the generality is a very simple formula: the distance is n minus c plus i over 2. No hurdles, no fortresses, and it is computed in linear time. Since every reversal is a DCJ, the DCJ distance never exceeds the reversal distance.",
                 FadeIn(form))
        self.hold()


# =====================================================================
class A6_Transpositions(RearrScene):
    """Slides 221 to 232: sorting by transpositions."""

    def construct(self):
        title = self.corner_title("Sorting by transpositions")
        self.add(title)
        self.say("Now we go back to one linear chromosome and change the operation. A transposition rho of i, j, k exchanges two adjacent blocks: the block from position i to j minus 1, and the block from j to k minus 1. Nothing is inverted, so we work with unsigned permutations. The index i is where the first block starts, j where the second starts, and k the first position after the second.")
        e = [0, 6, 5, 4, 3, 2, 1, 7]
        row = mkrow(e, False, y=2.2, x=-2.6)
        self.play(FadeIn(row))
        ops = [(1, 2, 3), (1, 4, 6), (2, 5, 7), (1, 3, 5)]
        cur = list(e)
        prev = row
        cs = [R.c_count(cur)]
        bt = [R.bt_count(cur)]
        self.set_subtitle("The fully reversed permutation sorted with transpositions. Each label shows the number of breakpoints bt and the number of cycles c.")
        self.add(T(f"bt={bt[0]}  c={cs[0]}", 17, C_HI).move_to([2.9, row.get_y(), 0], aligned_edge=LEFT))
        for (i, j, k) in ops:
            cur = R.transpx(cur, i, j, k)
            cs.append(R.c_count(cur))
            bt.append(R.bt_count(cur))
            prev, lab = self.transp_step(prev, cur, i, j, k, f"rho({i},{j},{k})", dy=0.62, label_x=0.1, label_size=18, rt=0.7)
            self.add(T(f"bt={bt[-1]}  c={cs[-1]}", 17, C_HI).move_to([2.9, prev.get_y(), 0], aligned_edge=LEFT))
        assert cur == list(range(8)) and cs == [1, 3, 3, 5, 7], cs
        self.say("Four transpositions, and three are not enough. This permutation is one reversal away from the identity, but four transpositions away. The two operations measure very different things.")
        self.wipe(title)

        lb = MathT(r"d_t(\pi)\ \ge\ \left\lceil \frac{b_t(\pi)}{3}\right\rceil", 38, YELLOW).move_to([-2.5, 2.3, 0])
        lb2 = MathT(r"d_t(\pi)\ \ge\ \left\lceil \frac{n+1-c(\pi)}{2}\right\rceil", 38, YELLOW).move_to([-2.5, 1.0, 0])
        self.say("The breakpoint definition changes slightly: a pair is an adjacency only if the second value is exactly the first plus 1. So 4 3 is now a breakpoint, because a transposition cannot turn a block around. A transposition cuts three places, so it removes at most 3 breakpoints. A transposition also changes the number of cycles by minus 2, 0 or plus 2, never by 1.",
                 FadeIn(lb), FadeIn(lb2))
        self.play(FadeOut(VGroup(lb, lb2)))

        pi = [0, 5, 2, 1, 4, 3, 6]
        steps = R.naive_transp(pi)
        assert len(steps) == 4
        left = mkrow(pi, False, y=2.1, x=-3.4, size=0.5)
        right = mkrow(pi, False, y=2.1, x=3.0, size=0.5)
        hdr = VGroup(T("naive: 4", 22, C_HI).next_to(left, UP, buff=0.2), T("optimal: 2", 22, C_GOOD).next_to(right, UP, buff=0.2))
        self.say("The naive algorithm is again selection sort on strips. Let i be the first breakpoint, j the position of the element that belongs there, and k the end of the strip that starts at j. Apply rho of i, j, k. It cuts only at breakpoints, creates an adjacency at position i, and the last step removes three, so it uses at most bt minus 2 transpositions, a factor below 3.",
                 FadeIn(left), FadeIn(right), FadeIn(hdr))
        prev = left
        for s in steps:
            i, j, k = s["i"], s["j"], s["k"]
            prev, _ = self.transp_step(prev, s["after"], i, j, k, f"rho({i},{j},{k}) bt={R.bt_count(s['after'])}", dy=0.58, label_x=-1.55, label_size=16, rt=0.7)
        cur = list(pi)
        prev = right
        for (i, j, k) in ((1, 3, 5), (2, 4, 6)):
            cur = R.transpx(cur, i, j, k)
            prev, _ = self.transp_step(prev, cur, i, j, k, f"rho({i},{j},{k})", dy=0.58, label_x=4.8, label_size=15, rt=0.7)
        assert cur == list(range(7))
        self.say("Two transpositions are enough, the lower bound. To finish in two, each must remove three breakpoints, so each must create three adjacencies at once: think in blocks, we want 1 after 0 and 4 5 together. The bounds are not always reached, though. For 4 3 1 5 6 2 both bounds say 2, but the distance is 3.")
        self.say("The history: Bafna and Pevzner gave factors 2, 1.75 and 1.5 in 1995. Elias and Hartman reached 1.375 in 2005, and in 2011 the problem was proved NP-hard. With only the cycle lower bound, no factor below 1.375 can be guaranteed, so going below needs a new lower bound.")


# =====================================================================
class A6_RevTransp(RearrScene):
    """Slides 233 to 258: reversals and transpositions together."""

    def construct(self):
        title = self.corner_title("Reversals and transpositions together")
        self.add(title)
        self.say("In nature both events happen, so the last model allows both operations. Having more tools never hurts, so the distance is at most the reversal distance, and without signs also at most the transposition distance. The lower bound for unsigned permutations is b over 3, since a reversal removes at most 2 breakpoints and a transposition at most 3.")
        pi = [0, 4, 6, 1, 7, 3, 2, 5, 8]
        steps = R.rt_unsigned(pi)
        assert [s["kind"] for s in steps] == ["tr", "rev", "tr", "rev"]
        assert [R.b_count(s["after"]) for s in steps] == [5, 3, 2, 0]
        row = mkrow(pi, False, y=2.3, x=-2.2, size=0.5)
        self.say("Walter, Dias and Meidanis showed that for every unsigned permutation other than the identity, a reversal or a transposition removes at least one breakpoint. The first strip is 0, 1, up to m. The element m plus 1 is either the first element of its strip, then a transposition brings the strip right after m, or the last element of a decreasing strip, then a reversal does it.",
                 FadeIn(row))
        prev = row
        hdr = T(f"b = {R.b_count(pi)}", 20, C_HI).move_to([3.6, row.get_y(), 0], aligned_edge=LEFT)
        self.add(hdr)
        for s in steps:
            m = s["m"]
            pre = box(VGroup(*row.cells[0:0]), GREEN) if False else None
            if s["kind"] == "rev":
                i, j = s["args"]
                self.set_subtitle(f"The prefix 0 to {m} is in place. The element {m + 1} ends a decreasing strip, so we use the reversal rho of {i} and {j}.")
                prev, lab = self.reverse_step(prev, s["after"], i, j, f"rev rho({i},{j})", dy=0.62, label_x=0.1, label_size=18, rt=0.8)
            else:
                i, j, k = s["args"]
                self.set_subtitle(f"The prefix 0 to {m} is in place. The element {m + 1} starts a strip, so we use the transposition rho of {i}, {j}, {k}.")
                prev, lab = self.transp_step(prev, s["after"], i, j, k, f"transp rho({i},{j},{k})", dy=0.62, label_x=0.1, label_size=18, rt=0.8)
            for c in prev.cells[:m + 1]:
                c.sq.set_fill(C_GOOD, 0.35)
            self.add(T(f"b = {R.b_count(s['after'])}", 20, C_HI).move_to([3.6, prev.get_y(), 0], aligned_edge=LEFT))
            self.hold()
        chain = MathT(r"\lceil b/3\rceil\le d_{rt}\le A\le b\ \Rightarrow\ A/d_{rt}\le 3", 34, YELLOW).move_to([0, -1.6, 0])
        self.say("Each iteration removes at least one breakpoint, so the algorithm uses at most b operations. The same recipe for the third time gives a 3-approximation in O of n squared. Here it used 4 operations while the slide says the optimum is 3.",
                 FadeIn(chain))
        self.wipe(title)

        # signed
        sp = [0, -4, 2, 5, -1, 3, 6]
        ops = [("rev", 1, 4), ("rev", 2, 3), ("tr", 3, 4, 5), ("tr", 3, 5, 6)]
        row = mkrow(sp, True, y=2.3, x=-2.2, size=0.5)
        cs = [R.c_count(sp)]
        cur = list(sp)
        self.say("With signs the same authors use the cycle graph. An oriented reversal always increases the number of cycles by one. When all elements are positive, a transposition exists that adds two cycles, or one that adds none followed by one that adds two. Transpositions keep signs, so once the negatives are gone only transpositions are used.",
                 FadeIn(row))
        self.add(T(f"c = {cs[0]}", 20, C_HI).move_to([3.6, row.get_y(), 0], aligned_edge=LEFT))
        prev = row
        for op in ops:
            if op[0] == "rev":
                _, i, j = op
                cur = R.revx(cur, i, j, True)
                prev, _ = self.reverse_step(prev, cur, i, j, f"rev rho({i},{j})", dy=0.62, label_x=0.1, label_size=18, rt=0.8)
            else:
                _, i, j, k = op
                cur = R.transpx(cur, i, j, k)
                prev, _ = self.transp_step(prev, cur, i, j, k, f"transp rho({i},{j},{k})", dy=0.62, label_x=0.1, label_size=18, rt=0.8)
            cs.append(R.c_count(cur))
            self.add(T(f"c = {cs[-1]}", 20, C_HI).move_to([3.6, prev.get_y(), 0], aligned_edge=LEFT))
            self.hold()
        assert cs == [2, 3, 4, 4, 6] and cur == list(range(7)), (cs, cur)
        self.say("The third step is a transposition that adds no cycle, and it prepares the fourth, which adds two. On average each operation adds at least one cycle, so the algorithm uses at most n plus 1 minus c operations, a 2-approximation. Here it used 4, exactly n plus 1 minus c, while the optimum is 3.")
        self.say("Later work: for unsigned permutations the factor is about 2.84, and in 2019 the problem was proved NP-hard with and without signs.")


# =====================================================================
class A6_Summary(RearrScene):
    """Slide 'whole lecture on one page': the common script and the table."""

    def construct(self):
        title = self.corner_title("The whole lecture on one page")
        self.add(title)
        script = [("1", "define a measure of disorder"), ("2", "show how much one operation can improve it"),
                  ("3", "obtain a lower bound"), ("4", "design an algorithm never too far from it")]
        sc = VGroup(*[T(f"{a}. {b}", 26, YELLOW) for a, b in script]).arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([-3.0, 1.2, 0])
        self.say("Almost every part of this lecture followed the same script: define a measure of disorder, show how much one operation can improve it, obtain a lower bound, and design an algorithm that is never too far from that bound.",
                 FadeIn(sc, lag_ratio=0.3), run_time=2.0)
        self.play(FadeOut(sc))
        rows = [
            ("Problem", "Lower bound", "Best known", "Complexity"),
            ("Reversals, unsigned", "b/2, b - cb", "factor ~1.375", "NP-hard"),
            ("Pancakes", "bp - 1", "factor 2", "NP-hard"),
            ("Reversals, signed", "n+1-c", "exact, + h + f", "polynomial"),
            ("DCJ", "exact", "n-(c+i/2)", "O(n)"),
            ("Transpositions", "bt/3, (n+1-c)/2", "factor 1.375", "NP-hard"),
            ("Rev+transp, unsigned", "b/3", "factor ~2.84", "NP-hard"),
            ("Rev+transp, signed", "(n+1-c)/2", "factor 2", "NP-hard"),
        ]
        xs = [-6.6, -2.8, 0.9, 4.6]
        tbl = VGroup()
        for r, row in enumerate(rows):
            for c, txt in enumerate(row):
                col = GREY_B if r == 0 else (C_GOOD if (r in (3, 4) and c == 3) else WHITE)
                t = T(txt, 21 if r else 22, col).move_to([xs[c], 2.3 - 0.52 * r, 0], aligned_edge=LEFT)
                tbl.add(t)
        self.say("Here is the table. The one row that is polynomial among the single chromosome problems is signed reversals. Knowing the orientation of the genes is what turns an NP-hard problem into one with an exact formula.",
                 FadeIn(tbl, lag_ratio=0.05), run_time=2.5)
        self.hold()
