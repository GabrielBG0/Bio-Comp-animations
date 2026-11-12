"""Aula 4: fragment assembly.

Scenes: A4_Shotgun, A4_OverlapGraph, A4_GreedySCS, A4_GreedyVsTopological,
        A4_DeBruijn, A4_KmerAssembly, A4_Consensus
"""
import random
from itertools import product

from common import *


# ---------------------------------------------------------------- pure helpers
def ov(a, b):
    """Largest t with suffix(a, t) == prefix(b, t), t < min(len)."""
    for t in range(min(len(a), len(b)) - 1, 0, -1):
        if a[-t:] == b[:t]:
            return t
    return 0


def all_ov(a, b):
    return [t for t in range(1, min(len(a), len(b))) if a[-t:] == b[:t]]


def glue(frags, order):
    s = frags[order[0]]
    ws = []
    for u, v in zip(order, order[1:]):
        t = ov(frags[u], frags[v])
        ws.append(t)
        s += frags[v][t:]
    return s, ws


def graph_positions(n, radius=2.1, center=ORIGIN, start=90):
    return [center + radius * np.array([np.cos((start - 360 * i / n) * DEGREES),
                                        np.sin((start - 360 * i / n) * DEGREES), 0]) for i in range(n)]


def node(label, pos, r=0.32, color=BLUE_C, font=26):
    c = Circle(radius=r).set_stroke(color, 3).set_fill(BLACK, 1)
    t = Text(str(label), font_size=font)
    g = VGroup(c, t).move_to(pos)
    g.circ = c
    return g


def edge(a, b, pos, r=0.32, arc=0.0, color=WHITE, width=3):
    p, q = pos[a], pos[b]
    d = normalize(q - p)
    e = Arrow(p + d * r, q - d * r, buff=0, path_arc=arc, thickness=width, tip_width_ratio=3).set_color(color)
    e.geom = (p + d * r, q - d * r, arc)
    return e


def edge_label(e, text, color=YELLOW, size=20, push=0.22):
    """Place a label beside the middle of an edge, on the bulging side if it is an arc."""
    a, b, arc = e.geom
    chord = b - a
    L = np.linalg.norm(chord)
    d = chord / L
    right = np.array([d[1], -d[0], 0.0])
    sag = (L / 2) * np.tan(arc / 4) if abs(arc) > 1e-6 else 0.0
    m = (a + b) / 2 + right * sag
    side = right * (1 if sag >= 0 else -1)
    if abs(arc) < 1e-6:
        side = -right
    t = Text(str(text), font_size=size, color=color)
    t.move_to(m + side * push)
    return t


def seq_text(s, size=26, color=WHITE, hl=None):
    """Monospace sequence; hl = {index: color}."""
    t = Text(s, font_size=size, color=color, t2c={}, )
    if hl:
        for i, c in hl.items():
            t[i].set_color(c)
    return t


# =====================================================================
class A4_Shotgun(NarratedScene):
    """Slides 1 to 4 and 12 to 15: why we have fragments, and how many we need."""

    def construct(self):
        title = self.corner_title("Fragment assembly: the problem and the coverage")
        self.add(title)
        self.say("Standard sequencing cannot read more than about 1000 base pairs in one piece, but the organisms we care about have millions. So the lab breaks many copies of the DNA at random positions and reads the small pieces.")

        T_len = 120
        f_len = 20
        sx = 0.095
        ox = -5.7

        def X(x):
            return ox + x * sx

        target = Line([X(0), 2.3, 0], [X(T_len), 2.3, 0]).set_stroke(WHITE, 8)
        tl = T("target sequence", 22, GREY_B).next_to(target, UP, buff=0.15)
        self.play(ShowCreation(target), FadeIn(tl))
        self.say("The pieces are the fragments. Fragment assembly rebuilds the target using the overlaps between fragments: if the end of one equals the beginning of another, they were probably neighbors.")

        palette = [BLUE_C, TEAL, YELLOW, ORANGE, PURPLE_C, GREEN]
        tmin = 3

        def run(n, seed, label):
            rng = random.Random(seed)
            starts = sorted(rng.randint(0, T_len - f_len) for _ in range(n))
            lanes = []
            laneof = []
            for s in starts:
                for k, last in enumerate(lanes):
                    if s > last + 1:
                        lanes[k] = s + f_len
                        laneof.append(k)
                        break
                else:
                    lanes.append(s + f_len)
                    laneof.append(len(lanes) - 1)
            # contigs: sweep with minimum overlap
            groups = []
            for s in starts:
                if groups and s <= groups[-1][1] - tmin:
                    groups[-1][1] = max(groups[-1][1], s + f_len)
                else:
                    groups.append([s, s + f_len])
            return starts, laneof, groups

        summary = VGroup()
        for n, seed in ((6, 5), (12, 4), (24, 9)):
            starts, laneof, groups = run(n, seed, "")
            c = n * f_len / T_len
            bars = VGroup()
            for s, ln in zip(starts, laneof):
                g = next(i for i, (a, b) in enumerate(groups) if a <= s < b)
                bar = Line([X(s), 1.6 - 0.28 * ln, 0], [X(s + f_len), 1.6 - 0.28 * ln, 0]).set_stroke(palette[g % 6], 6)
                bars.add(bar)
            gaps = VGroup()
            cov = [False] * T_len
            for s in starts:
                for x in range(s, s + f_len):
                    cov[x] = True
            x = 0
            while x < T_len:
                if not cov[x]:
                    y = x
                    while y < T_len and not cov[y]:
                        y += 1
                    gaps.add(Line([X(x), 2.3, 0], [X(y), 2.3, 0]).set_stroke(C_BAD, 8))
                    x = y
                else:
                    x += 1
            frac = sum(cov) / T_len
            txt = MathT(rf"n={n},\ c=\frac{{nf}}{{T}}={c:g}", 30, YELLOW).move_to([-4.2, -1.9, 0])
            res = T(f"covered: {100 * frac:.0f}%   contigs: {len(groups)}", 24, WHITE).move_to([0.2, -1.9, 0])
            exp = MathT(rf"1-e^{{-c}}={100 * (1 - np.exp(-c)):.1f}\%", 28, C_GOOD).move_to([4.4, -1.9, 0])
            if n == 6:
                self.set_subtitle("Average coverage c is the total amount read, n times f, divided by the target size T. Here we have 6 fragments of size 20 on a target of 120, so every base is read once on average.")
            elif n == 12:
                self.set_subtitle("Double the fragments and c becomes 2. The gaps shrink, but sampling is random, so some regions still have no fragment and the assembly breaks there.")
            else:
                self.set_subtitle("At coverage 4 almost everything is covered. A base is missed by all n fragments with probability about e to the minus c, so the covered fraction is 1 minus e to the minus c.")
            self.play(FadeIn(bars, lag_ratio=0.05), run_time=1.5)
            self.play(FadeIn(gaps), FadeIn(txt), FadeIn(res), FadeIn(exp))
            self.hold()
            self.play(FadeOut(bars), FadeOut(gaps), FadeOut(txt), FadeOut(res), FadeOut(exp))

        self.say("With 20,000 fragments of 600 base pairs on a 3 Mbp target, c is 4. At c equal to 1 about 63 percent of the target is covered, at 2 about 86.5, and at 4 about 98 percent.")
        f1 = MathT(r"p=n\,e^{-n(f-t)/T}", 38).move_to([0, 0.4, 0])
        self.say("The expected number of contigs is n times e to the minus n times f minus t over T. A fragment ends its contig when no other fragment starts in the next f minus t positions, and each contig has exactly one last fragment.",
                 Write(f1), run_time=1.5)


# =====================================================================
class A4_OverlapGraph(NarratedScene):
    """Slides 16 to 26: overlap multigraph, paths are superstrings."""

    F = ["ACTACAC", "CACTCAGGCA", "GCATTCACTA", "ACTAGAAATATA", "TATACCAGC"]

    def construct(self):
        title = self.corner_title("The overlap graph")
        self.add(title)
        F = self.F
        n = len(F)
        pos = graph_positions(n, radius=2.0, center=np.array([-2.6, -0.1, 0]))
        nodes = [node(i + 1, pos[i]) for i in range(n)]

        # fragment list on the right
        lst = VGroup(*[Text(f"{i + 1}  {f}", font_size=24, color=WHITE) for i, f in enumerate(F)])
        lst.arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to([4.2, 1.6, 0])
        self.say("To attack the shortest common superstring we store all overlap information in one structure. One vertex per fragment, and an edge from a to b with weight t when the last t characters of a equal the first t characters of b.",
                 FadeIn(lst), *[FadeIn(nd) for nd in nodes], run_time=1.5)

        edges = {}
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                ts = all_ov(F[i], F[j])
                if ts and max(ts) >= 2:
                    edges[(i, j)] = ts
        group = VGroup()
        pair_set = set(edges)
        eobjs = {}
        for (i, j), ts in edges.items():
            arc = 0.28 if (j, i) in pair_set else 0.0
            e = edge(i, j, pos, arc=arc, color=GREY_B)
            lab = edge_label(e, "/".join(map(str, sorted(ts, reverse=True))), size=20)
            eobjs[(i, j)] = (e, lab)
            group.add(e, lab)
        self.set_subtitle("Check one edge by hand. Fragment 3 ends with ACTA and fragment 1 starts with ACTA, so there is an edge from 3 to 1 with weight 4. Fragment 4 to 5 has two edges, weights 4 and 2, which is why this is a multigraph.")
        self.play(FadeIn(group, lag_ratio=0.05), run_time=2)
        e31 = eobjs[(2, 0)][0]
        self.play(e31.animate.set_color(YELLOW))
        self.hold()
        self.play(e31.animate.set_color(GREY_B))
        note = T("(edges of weight 1 and the weight-0 edges are not drawn)", 20, GREY_B).move_to([-2.6, -2.75, 0])
        self.play(FadeIn(note))

        # paths are superstrings
        self.say("A path in this graph is a recipe for gluing fragments: follow the path, and each time you cross an edge of weight t, write the next fragment without its first t characters.")
        paths = [([2, 0, 3, 4, 1], "greedy-like path 3, 1, 4, 5, 2"), ([0, 1, 2, 3, 4], "path 1, 2, 3, 4, 5")]
        out = VGroup()
        for order, name in paths:
            s, ws = glue(F, order)
            w = sum(ws)
            hl = [e for e in order]
            # highlight path edges
            hi = []
            for u, v in zip(order, order[1:]):
                if (u, v) in eobjs:
                    hi.append(eobjs[(u, v)][0].animate.set_color(YELLOW))
            self.set_subtitle(f"Take the path {', '.join(str(k + 1) for k in order)}. The edge weights are {', '.join(map(str, ws))}, so w(P) is {w}.")
            if hi:
                self.play(*hi, run_time=0.8)
            sline = T(s, 20, C_GOOD).move_to([4.0, -0.6, 0])
            if sline.get_width() > 6.5:
                sline.set_width(6.5)
            ident = MathT(rf"\lVert F(P)\rVert={sum(map(len, F))},\ w(P)={w},\ |S(P)|={len(s)}", 28, YELLOW).move_to([4.0, -1.6, 0])
            self.play(FadeIn(sline), FadeIn(ident))
            self.hold()
            self.play(FadeOut(sline), FadeOut(ident), *[eobjs[(u, v)][0].animate.set_color(GREY_B) for u, v in zip(order, order[1:]) if (u, v) in eobjs])
        ident2 = MathT(r"\lVert F(P)\rVert=w(P)+|S(P)|", 40, YELLOW).move_to([4.0, -0.6, 0])
        self.say("Always the same identity: the total fragment length equals the total overlap plus the length of the result. If P visits all the vertices, the left side is fixed, so a short superstring is the same as a heavy path. SCS is a maximum weight Hamiltonian path.",
                 Write(ident2), run_time=1.5)


# =====================================================================
class A4_GreedySCS(NarratedScene):
    """Slides 27 to 35: the greedy algorithm, with its trace."""

    F = A4_OverlapGraph.F

    def construct(self):
        title = self.corner_title("Greedy algorithm for SCS")
        self.add(title)
        F = self.F
        n = len(F)
        pos = graph_positions(n, radius=1.8, center=np.array([-4.3, 0.1, 0]))
        nodes = [node(i + 1, pos[i], r=0.28, font=22) for i in range(n)]

        E = sorted([(ov(F[i], F[j]), i, j) for i in range(n) for j in range(n) if i != j and ov(F[i], F[j]) > 0],
                   key=lambda e: (-e[0], e[1], e[2]))
        eobjs = {}
        group = VGroup()
        for w, i, j in E:
            if w < 2:
                continue
            arc = 0.25 if any((b, a) == (i, j) for _, a, b in E if _ >= 2) else 0.0
            e = edge(i, j, pos, r=0.28, arc=arc, color=GREY_D, width=2)
            eobjs[(i, j)] = e
            group.add(e)
        self.say("Between two vertices only the heaviest edge can help, so we keep just that one. This is the overlap graph. The greedy idea: at each step take the heaviest edge that is still valid.",
                 *[FadeIn(nd) for nd in nodes], FadeIn(group), run_time=1.5)
        rules = VGroup(T("valid edge:", 24, YELLOW),
                       T("no two chosen edges leave the same vertex", 22),
                       T("no two chosen edges enter the same vertex", 22),
                       T("no cycle is formed", 22)).arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([2.5, 2.4, 0])
        self.say("An edge is valid when no two chosen edges leave the same vertex, no two enter the same vertex, and no cycle is formed. That is exactly what makes the chosen edges a set of paths.",
                 FadeIn(rules), run_time=1.0)

        # state: pieces
        pieces = {i: [i] for i in range(n)}      # path through vertex -> list
        out_used, in_used = set(), set()
        first_of = {i: i for i in range(n)}
        last_of = {i: i for i in range(n)}
        strip = VGroup().move_to([2.5, -1.2, 0])
        log_pos = np.array([2.5, 0.6, 0])
        cur_pieces = [[i] for i in range(n)]

        def piece_strings():
            return [glue(F, p)[0] for p in cur_pieces]

        def make_strip():
            g = VGroup(*[Text(s, font_size=19, color=C_QUERY) for s in piece_strings()])
            g.arrange(DOWN, aligned_edge=LEFT, buff=0.1)
            for t in g:
                if t.get_width() > 7.5:
                    t.set_width(7.5)
            g.move_to([2.7, -1.35, 0])
            g.align_to([-0.6, 0, 0], LEFT)
            return g

        strip = make_strip()
        cap = T("current pieces", 20, GREY_B).next_to(strip, UP, buff=0.15).align_to(strip, LEFT)
        self.play(FadeIn(strip), FadeIn(cap))
        self.play(FadeOut(rules))

        chosen = []
        msg = VGroup()
        for w, i, j in E:
            if len(cur_pieces) == 1:
                break
            # find pieces
            pi = next(p for p in cur_pieces if p[-1] == i) if any(p[-1] == i for p in cur_pieces) else None
            pj = next(p for p in cur_pieces if p[0] == j) if any(p[0] == j for p in cur_pieces) else None
            ok = True
            why = ""
            if pi is None:
                ok, why = False, f"vertex {i + 1} already has an outgoing edge"
            elif pj is None:
                ok, why = False, f"vertex {j + 1} already has an incoming edge"
            elif pi is pj:
                ok, why = False, "it would close a cycle"
            if w < 2 and not ok:
                continue
            lbl = T(f"edge {i + 1} -> {j + 1}, weight {w}", 26, YELLOW).move_to([2.5, 1.6, 0])
            e = eobjs.get((i, j))
            anims = [FadeIn(lbl)]
            if e is not None:
                anims.append(e.animate.set_color(YELLOW))
            self.play(*anims, run_time=0.5)
            if not ok:
                r = T("invalid: " + why, 22, C_BAD).move_to([2.5, 1.0, 0])
                if e is not None:
                    self.play(FadeIn(r), e.animate.set_color(C_BAD), run_time=0.5)
                self.wait(0.9)
                self.play(FadeOut(r), FadeOut(lbl), *( [FadeOut(e)] if e is not None else []), run_time=0.4)
                continue
            cur_pieces.remove(pi)
            cur_pieces.remove(pj)
            cur_pieces.append(pi + pj)
            chosen.append((i, j, w))
            r = T(f"valid, overlap {F[i][-w:]}", 22, C_GOOD).move_to([2.5, 1.0, 0])
            newstrip = make_strip()
            newcap = cap
            self.play(FadeIn(r), *( [e.animate.set_color(C_GOOD)] if e is not None else []), run_time=0.4)
            self.play(Transform(strip, newstrip), run_time=0.8)
            self.wait(0.8)
            self.play(FadeOut(r), FadeOut(lbl), run_time=0.3)
        w_total = sum(c[2] for c in chosen)
        s_final = piece_strings()[0]
        self.set_subtitle("The result has w(P) equal to 4 plus 4 plus 3 plus 2, that is 13, so the superstring has 48 minus 13 equals 35 characters.")
        res = MathT(rf"w(P)={w_total},\quad |S(P)|=48-{w_total}={len(s_final)}", 34, C_GOOD).move_to([2.5, 1.4, 0])
        self.play(Write(res))
        self.hold()
        better = T("but the path 1,2,3,4,5 with weights 3,3,4,4 gives 34", 24, C_HI).move_to([2.5, -2.9, 0])
        better.set_y(-2.75)
        self.say("The greedy algorithm is not optimal. A shorter string of 34 characters exists. Taking the edge 3 to 1 first blocked both 3 to 4 and 2 to 3. A locally best choice destroyed a globally better solution, as expected for an NP-complete problem.",
                 FadeIn(better), run_time=1.0)
        cplx = MathT(r"\text{total cost: } O(\lVert F\rVert + n^{2})", 30, GREY_B).move_to([2.5, 2.4, 0])
        self.say("With prefix trees to build the graph, counting sort for the edges, and auxiliary vectors for the cycle test, the whole algorithm costs order of the total size plus n squared.",
                 FadeIn(cplx))


# =====================================================================
class A4_GreedyVsTopological(NarratedScene):
    """Slides 39 to 51: acyclic overlap graphs are solved exactly by a topological sort."""

    G = ["AGTATTGGCAATC", "AATCGATG", "ATGCAAACCT", "TTGGCAATCACT", "CCTTTTGG"]

    def construct(self):
        title = self.corner_title("Greedy against topological sort")
        self.add(title)
        G = self.G
        n = len(G)
        pos = [np.array(p) for p in ([-5.0, 1.0, 0], [-3.1, 1.9, 0], [-1.3, 1.0, 0], [-3.2, -1.2, 0], [-1.3, -0.4, 0])]
        pos = [np.array([-5.5, 0.9, 0]), np.array([-3.8, 1.9, 0]), np.array([-2.0, 0.9, 0]),
               np.array([-3.8, -1.4, 0]), np.array([-2.0, -0.4, 0])]
        nodes = [node(i + 1, pos[i], r=0.3, font=22) for i in range(n)]
        lst = VGroup(*[Text(f"{i + 1}  {f}", font_size=20) for i, f in enumerate(G)])
        lst.arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([3.6, 2.3, 0])
        self.say("Five fragments and a minimum overlap t equal to 3, so only edges of weight at least 3 survive in the graph. Small overlaps happen by chance, so we do not trust them.",
                 *[FadeIn(x) for x in nodes], FadeIn(lst), run_time=1.2)
        E = [(0, 1, 4), (1, 2, 3), (2, 4, 3), (4, 3, 4), (0, 3, 9)]
        for u, v, w in E:
            assert ov(G[u], G[v]) == w, (u, v, ov(G[u], G[v]))
        eo = {}
        grp = VGroup()
        for u, v, w in E:
            e = edge(u, v, pos, r=0.3, arc=0.0 if (u, v) != (0, 3) else 0.0, color=GREY_B)
            lab = edge_label(e, w, size=20)
            eo[(u, v)] = e
            grp.add(e, lab)
        self.play(FadeIn(grp, lag_ratio=0.2), run_time=1.5)

        # greedy
        self.set_subtitle("Greedy first grabs the heaviest edge, 1 to 4 with weight 9. Now 1 to 2 is invalid, a second edge leaving 1, and 5 to 4 is invalid, a second edge entering 4.")
        self.play(eo[(0, 3)].animate.set_color(C_GOOD))
        s14, _ = glue(G, [0, 3])
        t14 = T(s14, 20, C_GOOD).move_to([3.6, 0.5, 0])
        self.play(FadeIn(t14))
        self.play(eo[(0, 1)].animate.set_color(C_BAD), eo[(4, 3)].animate.set_color(C_BAD))
        self.hold()
        self.set_subtitle("It can still take 2 to 3 and 3 to 5. Then it is stuck with two pieces and no valid edge of weight 3 or more to join them.")
        self.play(eo[(1, 2)].animate.set_color(C_GOOD), eo[(2, 4)].animate.set_color(C_GOOD))
        s235, _ = glue(G, [1, 2, 4])
        t235 = T(s235, 20, C_GOOD).move_to([3.6, 0.0, 0])
        self.play(FadeIn(t235))
        stuck = T("two pieces, no valid edge left", 24, C_BAD).move_to([3.6, -0.7, 0])
        self.play(FadeIn(stuck))
        self.hold()
        self.play(FadeOut(VGroup(t14, t235, stuck)), *[eo[k].animate.set_color(GREY_B) for k in eo])

        # topological sort
        self.set_subtitle("Now the topological sort. The graph is acyclic. Vertex 1 is the only one with no incoming edge, so it comes first. Remove it and vertex 2 has no incoming edge. Continue, and the order is 1, 2, 3, 5, 4.")
        order = [0, 1, 2, 4, 3]
        removed = []
        for k, v in enumerate(order):
            ring = Circle(radius=0.38).set_stroke(YELLOW, 4).move_to(pos[v])
            num = T(str(k + 1), 24, YELLOW).next_to(ring, UP, buff=0.05)
            self.play(ShowCreation(ring), FadeIn(num), run_time=0.5)
            self.play(nodes[v].animate.set_opacity(0.35), run_time=0.3)
            for (a, b), e in eo.items():
                if a == v:
                    self.play(e.animate.set_opacity(0.2), run_time=0.15)
        self.hold()

        lines = VGroup()
        cur, _ = glue(G, order[:1])
        stack = []
        for k in range(1, len(order) + 1):
            s, ws = glue(G, order[:k])
            stack.append(s)
        text_lines = VGroup(*[Text(s, font_size=18, color=C_QUERY) for s in stack]).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        text_lines.move_to([3.5, -0.8, 0]).align_to([0.4, 0, 0], LEFT)
        self.set_subtitle("Gluing along that order gives the final sequence of 37 characters.")
        for t in text_lines:
            self.play(FadeIn(t), run_time=0.5)
        self.hold()
        final = stack[-1]
        rep = "TTGGCAATC"
        i1 = final.find(rep)
        i2 = final.find(rep, i1 + 1)
        assert i1 >= 0 and i2 > i1
        hl = Text(final, font_size=26)
        hl.move_to([0.5, -2.6, 0])
        for i in list(range(i1, i1 + 9)) + list(range(i2, i2 + 9)):
            hl[i].set_color(C_HI)
        self.play(FadeOut(text_lines), FadeIn(hl))
        self.say("Look at the final sequence. The string TTGGCAATC appears twice. It is a repeat of size 9. Fragment 1 ends inside the first copy and fragment 4 starts at the second, so they overlap without being neighbors. The heaviest edge was a lie caused by a repeat, and greedy believed it.",
                 run_time=0.5)


# =====================================================================
class A4_DeBruijn(NarratedScene):
    """Slides 57 and 58: the de Bruijn graph for the binary alphabet and k = 3."""

    def construct(self):
        title = self.corner_title("The de Bruijn graph, k = 3, binary alphabet")
        self.add(title)
        verts = ["".join(p) for p in product("01", repeat=3)]
        pos = {v: 2.3 * np.array([np.cos((90 - 45 * i) * DEGREES), np.sin((90 - 45 * i) * DEGREES), 0]) + np.array([-1.5, -0.5, 0])
               for i, v in enumerate(verts)}
        self.say("A k-mer is a string of k characters. In the de Bruijn graph the vertices are all the k-mers, and there is an edge from x to y when the last k minus 1 characters of x equal the first k minus 1 characters of y.")
        nodes = {v: node(v, pos[v], r=0.38, font=20) for v in verts}
        self.play(*[FadeIn(nd) for nd in nodes.values()])
        edges = {}
        grp = VGroup()
        for x in verts:
            for c in "01":
                y = x[1:] + c
                if x == y:
                    loop = Circle(radius=0.22).set_stroke(GREY_B, 3)
                    out = normalize(pos[x] - np.array([-1.5, -0.5, 0]))
                    loop.move_to(pos[x] + out * 0.62)
                    lab = T(c, 18, YELLOW).move_to(pos[x] + out * 1.0)
                    edges[(x, y)] = loop
                    grp.add(loop, lab)
                else:
                    rev = y[1:] == x[:2]
                    e = edge_to(pos, nodes, x, y, arc=0.3 if rev else 0.0)
                    lab = edge_label(e, c, size=18)
                    edges[(x, y)] = e
                    grp.add(e, lab)
        self.set_subtitle("Each vertex has one outgoing edge per character, so there are 2 to the 3 equals 8 vertices and 2 to the 4 equals 16 edges. The label on an edge is the appended character. For example 001 goes to 010 with label 0, and 000 has a loop.")
        self.play(FadeIn(grp, lag_ratio=0.03), run_time=3)
        self.hold()
        self.play(edges[("001", "010")].animate.set_color(YELLOW))
        info = VGroup(T("each edge = one string of size k + 1", 24, YELLOW),
                      T("it is Hamiltonian and Eulerian", 24, YELLOW)).arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to([4.2, 1.2, 0])
        self.say("Another way to see an edge: it is one string of size k plus 1, whose first k characters are x and whose last k characters are y. This graph is both Hamiltonian and Eulerian, and it is the setting of new technologies that produce short fragments of equal size.",
                 FadeIn(info), run_time=1.0)


def edge_to(pos, nodes, x, y, arc=0.0):
    p, q = pos[x], pos[y]
    d = normalize(q - p)
    e = Arrow(p + d * 0.38, q - d * 0.38, buff=0, path_arc=arc, thickness=2.5, tip_width_ratio=3).set_color(GREY_B)
    e.geom = (p + d * 0.38, q - d * 0.38, arc)
    return e


# =====================================================================
class A4_KmerAssembly(NarratedScene):
    """Slides 59 to 64: assembly as a walk through all edges of the k-mer graph (exercise 2)."""

    def construct(self):
        title = self.corner_title("Assembly with a k-mer graph")
        self.add(title)
        F = ["CTGCT", "CTCGAC", "CTCTCG", "ACTGC", "GCTCTC"]
        k = 3
        self.say("Now fragments become edges. Each fragment of size k plus 1 is an edge from its first k characters to its last k characters. An assembly is a walk through all the edges, which is polynomial, unlike a Hamiltonian path.")
        fl = VGroup(*[Text(f, font_size=26) for f in F]).arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([-5.2, 1.0, 0])
        self.play(FadeIn(fl))
        self.set_subtitle("Our fragments have 5 and 6 characters, so first cut each one into its windows of size 4.")
        windows = []
        wtxt = VGroup()
        for f in F:
            ws = [f[i:i + k + 1] for i in range(len(f) - k)]
            windows.append(ws)
        col = VGroup()
        for i, ws in enumerate(windows):
            t = Text("  ".join(ws), font_size=24, color=C_QUERY)
            t.next_to(fl[i], RIGHT, buff=0.5)
            col.add(t)
        self.play(FadeIn(col, lag_ratio=0.2), run_time=1.5)
        uniq = []
        for ws in windows:
            for w in ws:
                if w not in uniq:
                    uniq.append(w)
        self.set_subtitle("Removing duplicates leaves 9 different windows. Each becomes an edge from its first 3 characters to its last 3.")
        self.hold()
        self.play(FadeOut(fl), FadeOut(col))

        # graph layout (hand placed): main line + loop
        order = ["ACT", "CTG", "TGC", "GCT", "CTC", "TCG", "CGA", "GAC"]
        pos = {}
        for i, v in enumerate(order):
            pos[v] = np.array([-5.4 + 1.55 * i, 0.9, 0])
        pos["TCT"] = np.array([(pos["CTC"][0]), -0.9, 0])
        pos["TCT"] = np.array([pos["CTC"][0] + 0.78, -0.9, 0])
        nodes = {v: node(v, p, r=0.4, font=17) for v, p in pos.items()}
        self.play(*[FadeIn(nodes[v]) for v in nodes], run_time=0.8)
        eo = {}
        for w in uniq:
            a, b = w[:3], w[1:]
            arc = 0.0
            if {a, b} == {"CTC", "TCT"}:
                arc = 0.5 if a == "CTC" else 0.5
            e = Arrow(pos[a] + normalize(pos[b] - pos[a]) * 0.4, pos[b] - normalize(pos[b] - pos[a]) * 0.4, buff=0,
                      path_arc=arc, thickness=3, tip_width_ratio=3).set_color(GREY_B)
            eo[w] = e
        self.set_subtitle("Draw the graph. It is a line with one small loop. ACT has no incoming edge, so the walk starts there. GAC has no outgoing edge, so it ends there. CTC has two incoming and two outgoing edges.")
        self.play(*[ShowCreation(e) for e in eo.values()], run_time=2)
        self.hold()

        walk = ["ACT", "CTG", "TGC", "GCT", "CTC", "TCT", "CTC", "TCG", "CGA", "GAC"]
        self.set_subtitle("The walk is forced. From ACT we reach CTC. If we went on to TCG now we could never come back for the loop, so we go around the loop first.")
        spelled = Text(walk[0], font_size=38, color=YELLOW).move_to([0, -2.2, 0])
        cur_str = walk[0]
        self.play(FadeIn(spelled), nodes[walk[0]].circ.animate.set_stroke(YELLOW, 6))
        cur = walk[0]
        for nxt in walk[1:]:
            w = cur + nxt[-1]
            e = eo[w]
            cur_str += nxt[-1]
            new = Text(cur_str, font_size=38, color=YELLOW).move_to([0, -2.2, 0])
            self.play(e.animate.set_color(YELLOW), nodes[nxt].circ.animate.set_stroke(YELLOW, 6),
                      Transform(spelled, new), run_time=0.8)
            cur = nxt
        assert cur_str == "ACTGCTCTCGAC"
        self.set_subtitle("To spell the sequence, write the first vertex and then the last character of each following vertex. The assembly is ACTGCTCTCGAC, and every fragment is inside it. CTC is visited twice because the target has a small repeat there.")
        self.hold()
        small = T("with a single loop there is no choice, so this assembly is unique", 24, C_HI).move_to([0, -2.85, 0])
        small.set_y(-2.75)
        self.play(FadeIn(small))

        # second example: two loops give two assemblies
        self.play(FadeOut(VGroup(*nodes.values(), *eo.values(), spelled, small)))
        ex = VGroup(T("AAGA AGAC GACT ACTG CTGG TGGG GGGA GGAC", 22), T("ACTC CTCC TCCG CCGA CGAC ACTT CTTT", 22)).arrange(DOWN, buff=0.15).move_to([0, 2.5, 0])
        self.say("The slide 60 example has 15 fragments of size 4. The graph has a main line and two loops that leave ACT and come back to GAC. A walk through all edges must go around both loops, in either order.",
                 FadeIn(ex), run_time=1)
        a1 = "AAGACTCCGACTGGGACTTT"
        a2 = "AAGACTGGGACTCCGACTTT"
        t1 = Text(a1, font_size=30, color=C_QUERY).move_to([0, 0.6, 0])
        t2 = Text(a2, font_size=30, color=C_HI).move_to([0, -0.4, 0])
        for t, a in ((t1, a1), (t2, a2)):
            for i in [m for m in range(len(a) - 3) if a[m:m + 4] == "GACT"]:
                for j in range(i, i + 4):
                    t[j].set_color(YELLOW)
        self.say("Upper loop first, or lower loop first. Both assemblies have 20 characters and contain all 15 fragments. GACT appears three times in each. A repeat is a part of the graph that the walk crosses more than once, and the data cannot decide between the two.",
                 FadeIn(t1), FadeIn(t2), run_time=1.2)


# =====================================================================
class A4_Consensus(NarratedScene):
    """Slides 67, 68 and 75 to 77: Phred quality and the CAP3 weighted consensus."""

    def construct(self):
        title = self.corner_title("Phred quality and the CAP3 consensus")
        self.add(title)
        q = MathT(r"Q=-10\log_{10}P_e", 44).move_to([0, 2.3, 0])
        self.say("Phred attaches to each base the number Q, minus 10 times log base 10 of the probability that the base is wrong.", Write(q))
        rows = VGroup(*[VGroup(T(f"Q = {Q}", 28, C_QUERY), T(f"error {100 * 10 ** (-Q / 10):g}%", 28, C_HI), T(f"1 wrong base in {10 ** (Q // 10):,}", 28, GREY_B)).arrange(RIGHT, buff=0.8)
                        for Q in (10, 20, 30)]).arrange(DOWN, buff=0.25).move_to([0, 0.4, 0])
        self.say("Every 10 points divide the error by 10. Q 10 is 1 error in 10, Q 20 is 1 in 100, and Q 30 is 1 in 1,000.",
                 FadeIn(rows, lag_ratio=0.3), run_time=1.5)
        self.play(FadeOut(rows), FadeOut(q))

        self.say("CAP3 builds the consensus column by column. For each base type, add the qualities of the reads that show it. In each reading direction the highest quality gets weight 1 and the others get weight 0.5.")

        def column(fwd, rev, ypos, name):
            """fwd, rev: list of (base, quality)"""
            reads = VGroup()
            x0 = -5.8
            for k, (b, qv) in enumerate(fwd):
                reads.add(Text(f"{b}({qv}) ->", font_size=24, color=C_QUERY).move_to([x0 + 1.5 * k, ypos, 0]))
            for k, (b, qv) in enumerate(rev):
                reads.add(Text(f"<- {b}({qv})", font_size=24, color=C_HI).move_to([x0 + 1.5 * k, ypos - 0.55, 0]))
            return reads

        def weighted(reads):
            tot = {}
            parts = {}
            bases = sorted({b for b, _ in reads})
            for b in bases:
                qs = sorted([qv for bb, qv in reads if bb == b], reverse=True)
                tot[b] = qs[0] + sum(0.5 * x for x in qs[1:])
            return tot

        cols = [([("A", 30), ("A", 20), ("C", 10), ("A", 20)], [("A", 20), ("A", 10), ("T", 30)]),
                ([("A", 30), ("T", 30), ("T", 20), ("A", 20)], [("A", 20), ("A", 10), ("T", 30)])]
        base_y = 1.7
        for ci, (fwd, rev) in enumerate(cols):
            reads = column(fwd, rev, base_y, "")
            hdr = T(f"column {ci + 1}", 26, YELLOW).move_to([-6.0, base_y + 0.7, 0], aligned_edge=LEFT)
            self.set_subtitle(f"Column {ci + 1}. Arrows pointing right are forward reads, arrows pointing left are reverse reads, with the quality in parentheses.")
            self.play(FadeIn(hdr), FadeIn(reads, lag_ratio=0.1), run_time=1.5)
            self.hold()
            fw = weighted(fwd)
            rv = weighted(rev)
            bases = sorted(set(fw) | set(rv))
            tab = VGroup()
            for k, b in enumerate(bases):
                f_ = fw.get(b, 0)
                r_ = rv.get(b, 0)
                tab.add(VGroup(T(b, 28, YELLOW), T(f"forward {f_:g}", 24, C_QUERY), T(f"reverse {r_:g}", 24, C_HI),
                               T(f"total {f_ + r_:g}", 28, WHITE)).arrange(RIGHT, buff=0.6))
            tab.arrange(DOWN, aligned_edge=LEFT, buff=0.25).move_to([-1.5, -0.6, 0])
            tot = {b: fw.get(b, 0) + rv.get(b, 0) for b in bases}
            best = max(tot, key=tot.get)
            qual = tot[best] - sum(v for b, v in tot.items() if b != best)
            self.set_subtitle("Sum the weighted qualities per base. The base with the largest total goes to the consensus, and its quality is that total minus the totals of all the other bases.")
            self.play(FadeIn(tab, lag_ratio=0.3), run_time=1.5)
            res = T(f"consensus {best}, quality {tot[best]:g} - " + " - ".join(f"{v:g}" for b, v in tot.items() if b != best) + f" = {qual:g}", 28, C_GOOD).move_to([-1.0, -2.2, 0])
            self.play(FadeIn(res))
            if ci == 0:
                self.set_subtitle("In the first column almost everyone agrees, and the consensus A has quality 35.")
            else:
                self.set_subtitle("In the second column the reads are divided, T wins by little, and the quality 5 warns whoever uses the consensus that this base is doubtful.")
            self.hold()
            self.play(FadeOut(VGroup(hdr, reads, tab, res)))
        self.say("Several reads in the same direction are not fully independent evidence, since they may share a systematic error. So the best counts in full and the others count half. Agreement between the two directions is stronger, so each direction gets its own full weight.")
