"""Aula 5: phylogenetic trees.

Scenes: A5_AdditiveBasics, A5_AdditiveInsertion, A5_UPGMA, A5_NeighborJoining,
        A5_Sandwich, A5_PerfectPhylogeny, A5_Fitch
"""
from fractions import Fraction

from common import *


# ---------------------------------------------------------------- helpers
def fnum(x):
    """Format a Fraction / number without trailing zeros."""
    x = Fraction(x).limit_denominator(1000)
    if x.denominator == 1:
        return str(x.numerator)
    return f"{float(x):g}"


def table(rows, col_w=0.8, row_h=0.42, font=22, header=None, colors=None, row_labels=None):
    """Grid of Text cells. rows: list of lists of str. Returns VGroup with .cells[r][c]."""
    cells = []
    g = VGroup()
    nrows = len(rows)
    ncols = max(len(r) for r in rows)
    for r, row in enumerate(rows):
        cr = []
        for c, val in enumerate(row):
            col = WHITE
            if colors and (r, c) in colors:
                col = colors[(r, c)]
            t = Text(str(val), font_size=font, color=col)
            t.move_to([c * col_w, -r * row_h, 0])
            cr.append(t)
            g.add(t)
        cells.append(cr)
    g.cells = cells
    g.col_w, g.row_h = col_w, row_h
    return g


def matrix_view(names, D, font=22, cw=0.75, rh=0.42, colors=None, fill_diag=True):
    rows = [[""] + list(names)]
    for a in names:
        rows.append([a] + [("0" if (a == b and fill_diag) else ("" if a == b else fnum(D[a][b]))) for b in names])
    return table(rows, col_w=cw, row_h=rh, font=font, colors=colors)


def tnode(label, pos, r=0.28, color=BLUE_C, font=22, fill=True):
    c = Circle(radius=r).set_stroke(color, 3).set_fill(BLACK, 1 if fill else 0)
    t = Text(str(label), font_size=font)
    g = VGroup(c, t).move_to(pos)
    return g


def seg(p, q, color=WHITE, width=4):
    return Line(p, q).set_stroke(color, width)


def wlabel(p, q, text, color=YELLOW, size=20, off=0.22):
    m = (np.array(p) + np.array(q)) / 2
    d = np.array(q) - np.array(p)
    n = normalize(rotate_vector(d, PI / 2)) if np.linalg.norm(d) > 1e-6 else UP
    if n[1] < 0 or (abs(n[1]) < 1e-6 and n[0] < 0):
        n = -n
    t = Text(str(text), font_size=size, color=color)
    t.move_to(m + n * off)
    return t


# =====================================================================
class A5_AdditiveBasics(NarratedScene):
    """Slides 12 to 25: three objects always have a tree; the four-point condition."""

    def construct(self):
        title = self.corner_title("Additive trees: three and four objects")
        self.add(title)

        self.say("A distance matrix is additive if there is a tree whose leaf to leaf distances are exactly the matrix. With three objects there is always one: a star with a center.")
        a = tnode("a", [-5.0, 1.4, 0]); b = tnode("b", [-5.0, -1.4, 0]); c = tnode("c", [-0.5, 0.0, 0]); u = Dot([-2.5, 0.0, 0], radius=0.1, color=YELLOW)
        self.play(FadeIn(VGroup(a, b, c)), FadeIn(u))
        la = seg(a.get_center(), u.get_center()); lb = seg(b.get_center(), u.get_center()); lc = seg(c.get_center(), u.get_center())
        da = T("?", 26, YELLOW).move_to((a.get_center() + u.get_center()) / 2 + 0.35 * UP + 0.2 * LEFT)
        db = T("?", 26, YELLOW).move_to((b.get_center() + u.get_center()) / 2 + 0.35 * DOWN + 0.2 * LEFT)
        dc = T("?", 26, YELLOW).move_to((c.get_center() + u.get_center()) / 2 + 0.35 * UP)
        self.play(ShowCreation(la), ShowCreation(lb), ShowCreation(lc), FadeIn(VGroup(da, db, dc)))
        data = MathT(r"\delta(a,b)=63\quad \delta(a,c)=94\quad \delta(b,c)=79", 30, GREY_B).move_to([2.2, 2.6, 0])
        self.say("The three unknown branch lengths must add up to the three given distances. Here a to b is 63, a to c is 94, and b to c is 79.", FadeIn(data))
        f = MathT(r"d(x,u)=\frac{\delta(x,y)+\delta(x,z)-\delta(y,z)}{2}", 38).move_to([3.0, 1.0, 0])
        self.say("Solving the three equations gives a pattern: add the two distances that leave x, subtract the one that does not touch x, and halve. The parts that do not belong to the branch of x cancel.",
                 Write(f), run_time=1.5)
        vals = [("a", (63 + 94 - 79) // 2, da), ("b", (63 + 79 - 94) // 2, db), ("c", (94 + 79 - 63) // 2, dc)]
        calc = VGroup(MathT(r"d(a,u)=\frac{63+94-79}{2}=39", 30), MathT(r"d(b,u)=\frac{63+79-94}{2}=24", 30), MathT(r"d(c,u)=\frac{94+79-63}{2}=55", 30))
        calc.arrange(DOWN, aligned_edge=LEFT, buff=0.3).move_to([3.0, -0.8, 0])
        anims = []
        for (nm, v, lab), cl in zip(vals, calc):
            self.set_subtitle(f"For {nm}: the result is {v}.")
            new = T(str(v), 28, C_GOOD).move_to(lab)
            self.play(FadeIn(cl), Transform(lab, new), run_time=1.0)
            self.wait(0.6)
        self.say("Check: 39 plus 24 is 63, 39 plus 55 is 94, and 24 plus 55 is 79. All three distances are matched.")
        self.play(FadeOut(VGroup(a, b, c, u, la, lb, lc, da, db, dc, data, f, calc)))

        # four point condition
        names = list("abcd")
        D = {"a": {"b": 5, "c": 7, "d": 13}, "b": {"a": 5, "c": 8, "d": 14}, "c": {"a": 7, "b": 8, "d": 8}, "d": {"a": 13, "b": 14, "c": 8}}
        mv = matrix_view(names, D, font=26, cw=0.8, rh=0.5).move_to([-4.6, 0.6, 0])
        self.say("For four objects we need a test. The four point condition says: for every group of four, compute the three sums of opposite pairs. The two largest sums must be equal.",
                 FadeIn(mv))
        sums = [("ab + cd", 5 + 8 if False else 5 + 8), ("ac + bd", 7 + 14), ("ad + bc", 13 + 8)]
        sums = [("ab + cd", D["a"]["b"] + D["c"]["d"]), ("ac + bd", D["a"]["c"] + D["b"]["d"]), ("ad + bc", D["a"]["d"] + D["b"]["c"])]
        lines = VGroup()
        for k, (nm, v) in enumerate(sums):
            lines.add(T(f"{nm} = {v}", 30, WHITE).move_to([2.0, 1.5 - 0.7 * k, 0], aligned_edge=LEFT))
        self.set_subtitle("Here ab plus cd is 13, ac plus bd is 21, and ad plus bc is 21. The two largest are equal, so the group is fine.")
        self.play(FadeIn(lines, lag_ratio=0.4), run_time=1.5)
        box = SurroundingRectangle(VGroup(lines[1], lines[2]), buff=0.12).set_stroke(C_GOOD, 3)
        self.play(ShowCreation(box))
        self.hold()
        # the tree behind it
        tr = VGroup()
        pa, pb, pc, pd = [-0.8, -1.0, 0], [-0.8, -2.6, 0], [3.6, -1.0, 0], [3.6, -2.6, 0]
        pu, pv = [0.4, -1.8, 0], [2.4, -1.8, 0]
        tr.add(seg(pa, pu), seg(pb, pu), seg(pu, pv), seg(pv, pc), seg(pv, pd))
        tr.add(T("a", 22).move_to(pa), T("b", 22).move_to(pb), T("c", 22).move_to(pc), T("d", 22).move_to(pd))
        mid = T("(21 - 13) / 2 = 4", 22, YELLOW).move_to([1.4, -1.45, 0])
        mid2 = T("(sum of the two largest - smallest) / 2 = middle edge", 20, GREY_B)
        self.say("The smallest sum tells us who the neighbors are: a pairs with b, and c with d. And the middle edge of the tree is the largest sum minus the smallest, divided by 2.",
                 FadeIn(tr), FadeIn(mid), run_time=1.0)


# =====================================================================
def three_point_t(D, x, y, z):
    t = Fraction(D[x][z] + D[x][y] - D[y][z], 2)
    l = Fraction(D[x][z] + D[y][z] - D[x][y], 2)
    return t, l


def tree_path(adj, s, t):
    stack = [(s, [s])]
    seen = {s}
    while stack:
        v, p = stack.pop()
        if v == t:
            return p
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append((w, p + [w]))
    return None


class A5_AdditiveInsertion(NarratedScene):
    """Slides 26 to 42: build the additive tree inserting one object at a time."""

    D = {"a": {"b": 5, "c": 7, "d": 13, "e": 9}, "b": {"a": 5, "c": 8, "d": 14, "e": 10},
         "c": {"a": 7, "b": 8, "d": 8, "e": 4}, "d": {"a": 13, "b": 14, "c": 8, "e": 6},
         "e": {"a": 9, "b": 10, "c": 4, "d": 6}}
    POS = {"a": [-6.0, 1.0, 0], "b": [-6.0, -1.2, 0], "u1": [-4.4, -0.1, 0], "u2": [-2.6, -0.1, 0], "c": [-2.6, -1.9, 0],
           "u3": [-0.8, -0.1, 0], "d": [0.8, 1.0, 0], "e": [0.8, -1.2, 0]}

    def construct(self):
        title = self.corner_title("Additive tree: insert one object at a time")
        self.add(title)
        D = self.D
        names = list("abcde")
        mv = matrix_view(names, D, font=18, cw=0.55, rh=0.32).move_to([4.8, 2.35, 0])
        pos = {k: np.array(v, dtype=float) for k, v in self.POS.items()}

        # ---- run the algorithm for real, recording a trace
        adj = {"a": {"b": Fraction(5)}, "b": {"a": Fraction(5)}}
        leaves = ["a", "b"]
        counter = [0]
        trace = []

        def snapshot():
            es = {}
            for u in adj:
                for v, w in adj[u].items():
                    es[frozenset((u, v))] = w
            return es

        trace.append(("start", snapshot(), []))
        for z in names[2:]:
            x, y = "a", "b"
            tried = []
            ycands = [l for l in leaves if l != x]
            idx = 0
            while True:
                y = ycands[idx]
                t, l = three_point_t(D, x, y, z)
                path = tree_path(adj, x, y)
                acc = Fraction(0)
                landing = None
                for u, v in zip(path, path[1:]):
                    w = adj[u][v]
                    if acc < t < acc + w:
                        landing = ("edge", u, v, t - acc, acc + w - t)
                        break
                    if acc + w == t:
                        landing = ("vertex", v)
                        break
                    acc += w
                tried.append((x, y, t, l, landing))
                if landing[0] == "vertex":
                    c = landing[1]
                    nxt = None
                    for j in range(idx + 1, len(ycands)):
                        pth = tree_path(adj, x, ycands[j])
                        if c in pth and ycands[j] != c:
                            nxt = j
                            break
                    if nxt is None:
                        adj.setdefault(z, {})
                        adj[c][z] = l
                        adj[z][c] = l
                        trace.append((z, snapshot(), tried))
                        break
                    idx = nxt
                    continue
                _, u, v, w1, w2 = landing
                counter[0] += 1
                nn = f"u{counter[0]}"
                del adj[u][v]
                del adj[v][u]
                adj[nn] = {u: w1, v: w2}
                adj[u][nn] = w1
                adj[v][nn] = w2
                adj[z] = {nn: l}
                adj[nn][z] = l
                trace.append((z, snapshot(), tried))
                break
            leaves.append(z)
        # verify: leaf distances equal the matrix
        for p in names:
            for q in names:
                if p < q:
                    pth = tree_path(adj, p, q)
                    assert sum(adj[u][v] for u, v in zip(pth, pth[1:])) == D[p][q], (p, q)
        assert counter[0] == 3

        # ---- animate it
        node_m, edge_m, lab_m = {}, {}, {}

        def mk_node(k):
            if k.startswith("u"):
                return Dot(pos[k], radius=0.09, color=YELLOW)
            return tnode(k, pos[k], r=0.26, font=22)

        def sync(es, run=0.8):
            new_nodes = set()
            for key in es:
                new_nodes |= set(key)
            outs, ins = [], []
            for key in list(edge_m):
                if key not in es:
                    outs += [FadeOut(edge_m.pop(key)), FadeOut(lab_m.pop(key))]
            for k in new_nodes:
                if k not in node_m:
                    node_m[k] = mk_node(k)
                    ins.append(FadeIn(node_m[k]))
            for key, w in es.items():
                u, v = tuple(key)
                if key not in edge_m:
                    pu, pv = pos[u], pos[v]
                    # shorten a bit so the line does not cross the node circles
                    e = seg(pu, pv, GREY_B, 4)
                    lab = wlabel(pu, pv, fnum(w), size=20)
                    edge_m[key], lab_m[key] = e, lab
                    lab._txt = fnum(w)
                    ins += [ShowCreation(e), FadeIn(lab)]
            if outs or ins:
                self.play(*outs, *ins, run_time=run)
            self.bring_to_front(*node_m.values())

        self.say("The algorithm starts with the tree of three objects and inserts the others one at a time. To insert z, pick two leaves x and y already in the tree and apply the three object formula.",
                 FadeIn(mv), run_time=1.0)
        f = MathT(r"t=\frac{\delta(x,z)+\delta(x,y)-\delta(y,z)}{2},\qquad \ell=\frac{\delta(x,z)+\delta(y,z)-\delta(x,y)}{2}", 26, YELLOW)
        f.move_to([-1.6, 2.65, 0])
        self.say("The number t says how far from x, along the path to y, the branch of z leaves. The number l is the length of that branch.", FadeIn(f), run_time=1.0)

        sync(trace[0][1])
        note = T("a and b: one edge of weight 5", 24, WHITE).move_to([4.0, -0.2, 0])
        self.say("Start with a and b, joined by one edge of weight 5.", FadeIn(note))
        self.play(FadeOut(note))

        for z, snap, tried in trace[1:]:
            self.set_subtitle(f"Insert {z}.")
            for (x, y, t, l, landing) in tried:
                txt = MathT(rf"x={x},\ y={y}:\ t=\frac{{{D[x][z]}+{D[x][y]}-{D[y][z]}}}{{2}}={fnum(t)}", 26)
                txt.move_to([4.0, -0.3, 0])
                path = tree_path({k: dict(v) for k, v in self._adj_before(trace, z).items()}, x, y)
                # highlight the path x -> y
                hl = []
                for p, q in zip(path, path[1:]):
                    key = frozenset((p, q))
                    if key in edge_m:
                        hl.append(edge_m[key].animate.set_color(YELLOW))
                if landing[0] == "vertex":
                    res = T(f"on vertex {landing[1]}: try another y", 22, C_HI)
                    self.set_subtitle(f"With y equal to {y}, t is {fnum(t)}. That is exactly the vertex {landing[1]}. The branch of {z} may leave from there or from a subtree hanging there, so we try another leaf.")
                else:
                    res = T(f"inside an edge: split {fnum(landing[3])} + {fnum(landing[4])}", 22, C_GOOD)
                    self.set_subtitle(f"With y equal to {y}, t is {fnum(t)}. That falls in the middle of an edge, so we split the edge there and hang {z} with weight l equal to {fnum(l)}.")
                res.move_to([4.0, -1.0, 0])
                self.play(FadeIn(txt), *hl, run_time=0.8)
                self.play(FadeIn(res), run_time=0.5)
                self.hold()
                self.play(FadeOut(txt), FadeOut(res), *[edge_m[frozenset((p, q))].animate.set_color(GREY_B) for p, q in zip(path, path[1:]) if frozenset((p, q)) in edge_m], run_time=0.4)
            # update edge labels and structure
            old_keys = set(edge_m)
            sync(snap, run=1.0)
            # refresh weights on surviving edges
            self._refresh_labels(snap, edge_m, lab_m, pos)
            self.wait(0.5)

        final = T("a to e = 2 + 4 + 2 + 1 = 9,   b to d = 3 + 4 + 2 + 5 = 14,   c to e = 1 + 2 + 1 = 4", 22, C_GOOD).move_to([0.5, -2.7, 0])
        self.say("Check the final tree. From a to e the path weighs 2 plus 4 plus 2 plus 1, which is 9. From b to d it is 14, and from c to e it is 4. All ten distances match the matrix.",
                 FadeIn(final), run_time=1.0)
        cx = MathT(r"\Theta(n^{3})", 40, GREY_B).move_to([5.0, -1.3, 0])
        self.say("Inserting the k-th object may try up to k minus 2 leaves, and each try walks a path of size k. Summing over k gives order n cubed. The additive tree of an additive matrix is unique.", FadeIn(cx))

    def _adj_before(self, trace, z):
        """Adjacency of the tree just before inserting z."""
        prev = None
        for zz, snap, _ in trace:
            if zz == z:
                break
            prev = snap
        adj = {}
        for key, w in prev.items():
            u, v = tuple(key)
            adj.setdefault(u, {})[v] = w
            adj.setdefault(v, {})[u] = w
        return adj

    def _refresh_labels(self, snap, edge_m, lab_m, pos):
        anims = []
        for key, w in snap.items():
            u, v = tuple(key)
            new = wlabel(pos[u], pos[v], fnum(w), size=20)
            if lab_m[key].get_center()[0] != new.get_center()[0] or True:
                old_txt = getattr(lab_m[key], "_txt", None)
                if old_txt != fnum(w):
                    anims.append(Transform(lab_m[key], new))
                    lab_m[key]._txt = fnum(w)
        if anims:
            self.play(*anims, run_time=0.6)


# =====================================================================
class A5_UPGMA(NarratedScene):
    """Slides 60 to 72: UPGMA on an ultrametric matrix."""

    def construct(self):
        title = self.corner_title("UPGMA")
        self.add(title)
        names = list("abcde")
        D0 = {"a": {"b": 12, "c": 12, "d": 12, "e": 12}, "b": {"a": 12, "c": 4, "d": 6, "e": 6},
              "c": {"a": 12, "b": 4, "d": 6, "e": 6}, "d": {"a": 12, "b": 6, "c": 6, "e": 2},
              "e": {"a": 12, "b": 6, "c": 6, "d": 2}}
        self.say("UPGMA clusters from the bottom up. At the start each object is a group. At each iteration join the two closest groups, where the distance between groups is the average of the distances between their objects.")

        # dendrogram geometry
        leaf_x = {"a": 1.8, "b": 2.9, "c": 3.7, "d": 4.8, "e": 5.6}
        base_y = -2.3
        sc = 0.42

        def hy(h):
            return base_y + h * sc

        # state
        groups = {n: (n,) for n in names}               # name -> members
        dist = {frozenset((x, y)): Fraction(D0[x][y]) for x in names for y in names if x < y}
        height = {n: Fraction(0) for n in names}
        xpos = dict(leaf_x)
        cur_names = list(names)
        gm = {}

        def nm(g):
            return "".join(g)

        # leaf marks
        leaf_mobs = {n: Text(n, font_size=26).move_to([leaf_x[n], base_y - 0.3, 0]) for n in names}
        self.play(*[FadeIn(m) for m in leaf_mobs.values()])
        axis = Line([1.2, base_y, 0], [6.3, base_y, 0]).set_stroke(GREY_D, 2)
        self.add(axis)

        def build_matrix():
            ns = cur_names
            Dm = {a: {b: dist[frozenset((a, b))] for b in ns if b != a} for a in ns}
            labels = [nm(groups[n]) for n in ns]
            rows = [[""] + labels]
            for a in ns:
                rows.append([nm(groups[a])] + [("0" if a == b else fnum(dist[frozenset((a, b))])) for b in ns])
            m = table(rows, col_w=0.85, row_h=0.42, font=22)
            m.move_to([-4.2, 1.4, 0])
            return m

        mat = build_matrix()
        self.play(FadeIn(mat))
        self.say("The matrix of slide 64 is ultrametric: any three objects have two equal largest distances.")
        edges_drawn = VGroup()
        counter = 0
        while len(cur_names) > 1:
            pairs = sorted(((dist[frozenset((x, y))], x, y) for i, x in enumerate(cur_names) for y in cur_names[i + 1:]), key=lambda p: (p[0], p[1], p[2]))
            dmin, x, y = pairs[0]
            # highlight entry
            i, j = cur_names.index(x) + 1, cur_names.index(y) + 1
            r1 = SurroundingRectangle(mat.cells[i][j + 0], buff=0.08).set_stroke(YELLOW, 3)
            r2 = SurroundingRectangle(mat.cells[j][i + 0], buff=0.08).set_stroke(YELLOW, 3)
            self.set_subtitle(f"The smallest entry is dist({nm(groups[x])}, {nm(groups[y])}) equal to {fnum(dmin)}. Join these two groups.")
            self.play(ShowCreation(r1), ShowCreation(r2), run_time=0.6)
            self.hold()
            h = dmin / 2
            counter += 1
            newn = f"G{counter}"
            members = groups[x] + groups[y]
            cx = (xpos[x] + xpos[y]) / 2
            # edges (vertical + horizontal)
            hx, hyv = height[x], height[y]
            e1 = VGroup(seg([xpos[x], hy(hx), 0], [xpos[x], hy(h), 0]), seg([xpos[y], hy(hyv), 0], [xpos[y], hy(h), 0]), seg([xpos[x], hy(h), 0], [xpos[y], hy(h), 0]))
            w1 = T(fnum(h - hx), 20, YELLOW).move_to([xpos[x] - 0.25, (hy(hx) + hy(h)) / 2, 0])
            w2 = T(fnum(h - hyv), 20, YELLOW).move_to([xpos[y] + 0.25, (hy(hyv) + hy(h)) / 2, 0])
            hlab = MathT(rf"\text{{height}}={fnum(dmin)}/2={fnum(h)}", 24, C_GOOD).move_to([-4.2, -1.3, 0])
            self.set_subtitle(f"The new vertex gets height equal to the distance divided by 2, here {fnum(h)}. In an ultrametric tree two leaves at distance d have their common ancestor at height d over 2. Each new edge weighs the difference of heights.")
            self.play(FadeIn(e1), FadeIn(w1), FadeIn(w2), FadeIn(hlab), run_time=1.2)
            self.hold()
            # update table
            others = [o for o in cur_names if o not in (x, y)]
            sx, sy = len(groups[x]), len(groups[y])
            for o in others:
                dist[frozenset((newn, o))] = (sx * dist[frozenset((x, o))] + sy * dist[frozenset((y, o))]) / (sx + sy)
            for k in [k for k in dist if x in k or y in k]:
                dist.pop(k, None)
            groups[newn] = members
            height[newn] = h
            xpos[newn] = cx
            cur_names = [o for o in cur_names if o not in (x, y)] + [newn]
            self.play(FadeOut(r1), FadeOut(r2), FadeOut(hlab))
            if len(cur_names) > 1:
                newmat = build_matrix()
                self.set_subtitle("Update the table. The new distance to each other group is the average over all pairs of original objects, one in each group, so each group is weighted by its size.")
                self.play(Transform(mat, newmat), run_time=1.0)
                self.hold()
        root = T(f"root at height {fnum(height[cur_names[0]])}", 24, C_GOOD).move_to([6.0, hy(height[cur_names[0]]) + 0.3, 0])
        root.set_x(4.0)
        self.play(FadeIn(root))
        self.say("The final tree: every leaf is at distance 6 from the root. Check two distances. d to e is 1 plus 1 equals 2, and b to e is 2 plus 1 plus 2 plus 1 equals 6. Both match the matrix.")
        wp = MathT(r"\text{WPGMA: }\ \mathrm{dist}(D,C)=\frac{\mathrm{dist}(C_1,C)+\mathrm{dist}(C_2,C)}{2}", 28, GREY_B).move_to([-3.0, -2.5, 0])
        self.say("WPGMA differs in one line: it takes the plain mean of the two table entries, so a lone object counts as much as three. For an ultrametric matrix both give the same tree.",
                 FadeIn(wp))


# =====================================================================
class A5_NeighborJoining(NarratedScene):
    """Slides 73 to 104: Neighbor-Joining, and why UPGMA fails when rates differ."""

    def construct(self):
        title = self.corner_title("Neighbor-Joining")
        self.add(title)
        names = ["a", "b", "c", "d"]
        D0 = {"a": {"b": 6, "c": 3, "d": 8}, "b": {"a": 6, "c": 7, "d": 12}, "c": {"a": 3, "b": 7, "d": 7}, "d": {"a": 8, "b": 12, "c": 7}}
        mv = matrix_view(names, D0, font=24, cw=0.75, rh=0.45).move_to([-5.4, 2.0, 0])
        self.say("Why not simply join the closest pair like UPGMA? Without a molecular clock, the closest pair may not be neighbors in the tree. Two objects on short branches can be close without being siblings.", FadeIn(mv))

        u = {x: Fraction(sum(D0[x].values()), len(names) - 2) for x in names}
        ut = VGroup(*[MathT(rf"u({x})=\frac{{{sum(D0[x].values())}}}{{{len(names) - 2}}}={fnum(u[x])}", 26) for x in names])
        ut.arrange(DOWN, aligned_edge=LEFT, buff=0.18).move_to([-4.6, -0.6, 0])
        self.say("NJ computes for each object a separation u, the sum of its distances divided by the number of nodes minus 2. It measures how far x is from everybody else.", FadeIn(ut, lag_ratio=0.3), run_time=1.5)
        sform = MathT(r"S(i,j)=\delta(i,j)-u(i)-u(j)", 34, YELLOW).move_to([2.2, 2.7, 0])
        self.say("Then it joins the pair that minimizes S. So NJ looks for two objects that are close to each other and at the same time far from all the rest.", Write(sform))
        pairs = [("a", "b"), ("a", "c"), ("a", "d"), ("b", "c"), ("b", "d"), ("c", "d")]
        S = {p: D0[p[0]][p[1]] - u[p[0]] - u[p[1]] for p in pairs}
        rows = VGroup()
        for k, p in enumerate(pairs):
            rows.add(MathT(rf"S({p[0]},{p[1]})={D0[p[0]][p[1]]}-{fnum(u[p[0]])}-{fnum(u[p[1]])}={fnum(S[p])}", 24))
        rows.arrange(DOWN, aligned_edge=LEFT, buff=0.18).move_to([2.6, 0.5, 0])
        self.set_subtitle("Compute S for the six pairs. The minimum, minus 15, is shared by a b and c d. Ties may be broken either way.")
        self.play(FadeIn(rows, lag_ratio=0.3), run_time=2.0)
        mn = min(S.values())
        marks = VGroup(*[SurroundingRectangle(rows[k], buff=0.06).set_stroke(C_GOOD, 3) for k, p in enumerate(pairs) if S[p] == mn])
        self.play(ShowCreation(marks))
        self.hold()
        self.play(FadeOut(rows), FadeOut(marks), FadeOut(sform), FadeOut(ut))

        # join a,b
        d_ab_a = (D0["a"]["b"] + u["a"] - u["b"]) / 2
        d_ab_b = D0["a"]["b"] - d_ab_a
        f1 = MathT(rf"d(N_1,a)=\frac{{6+{fnum(u['a'])}-{fnum(u['b'])}}}{{2}}={fnum(d_ab_a)},\quad d(N_1,b)=6-{fnum(d_ab_a)}={fnum(d_ab_b)}", 26).move_to([2.2, 2.0, 0])
        self.say("Join a and b into a new node N1. The node that is farther from everyone gets the longer edge: a gets 1 and b gets 5.", FadeIn(f1))
        dN1c = (D0["a"]["c"] + D0["b"]["c"] - D0["a"]["b"]) / 2
        dN1d = (D0["a"]["d"] + D0["b"]["d"] - D0["a"]["b"]) / 2
        f2 = MathT(rf"\mathrm{{dist}}(N_1,c)=\frac{{3+7-6}}{{2}}={fnum(dN1c)},\quad \mathrm{{dist}}(N_1,d)=\frac{{8+12-6}}{{2}}={fnum(dN1d)}", 26).move_to([2.2, 1.0, 0])
        self.say("The new distances come from the three object formula: the distance from c to the point where the branches of a and b meet.", FadeIn(f2))
        # tree so far
        pa, pb, pN1 = [-1.0, -0.6, 0], [-1.0, -2.2, 0], [0.6, -1.4, 0]
        tr1 = VGroup(seg(pa, pN1), seg(pb, pN1), tnode("a", pa, r=0.25), tnode("b", pb, r=0.25), Dot(pN1, radius=0.09, color=YELLOW))
        tr1.add(wlabel(pa, pN1, "1"), wlabel(pb, pN1, "5"))
        self.play(FadeIn(tr1))
        self.hold()
        self.play(FadeOut(f1), FadeOut(f2))

        # second round: N1, c, d
        D1 = {"N1": {"c": dN1c, "d": dN1d}, "c": {"N1": dN1c, "d": Fraction(7)}, "d": {"N1": dN1d, "c": Fraction(7)}}
        u1 = {x: sum(D1[x].values()) / 1 for x in D1}
        S1 = {("N1", "c"): D1["N1"]["c"] - u1["N1"] - u1["c"], ("N1", "d"): D1["N1"]["d"] - u1["N1"] - u1["d"], ("c", "d"): D1["c"]["d"] - u1["c"] - u1["d"]}
        assert all(v == -16 for v in S1.values()), S1
        t2 = T("three nodes left: u(N1) = 9, u(c) = 9, u(d) = 14\nevery pair has S = -16", 22, WHITE).move_to([3.2, 2.2, 0])
        self.say("Three nodes remain, with u of 9, 9 and 14, and every pair has S equal to minus 16. Join c and d into N2.", FadeIn(t2))
        f3 = MathT(r"d(N_2,c)=\frac{7+9-14}{2}=1,\quad d(N_2,d)=7-1=6", 26).move_to([2.6, 1.3, 0])
        f4 = MathT(r"\mathrm{dist}(N_2,N_1)=\frac{2+7-7}{2}=1", 26).move_to([2.6, 0.5, 0])
        self.say("c gets 1 and d gets 6. The last edge, between N1 and N2, has weight 1, the distance of the two remaining nodes.", FadeIn(f3), FadeIn(f4))
        pc, pd, pN2 = [3.2, -0.6, 0], [3.2, -2.2, 0], [1.9, -1.4, 0]
        tr2 = VGroup(seg(pN1, pN2), seg(pc, pN2), seg(pd, pN2), tnode("c", pc, r=0.25), tnode("d", pd, r=0.25), Dot(pN2, radius=0.09, color=YELLOW))
        tr2.add(wlabel(pN1, pN2, "1", off=0.22), wlabel(pc, pN2, "1"), wlabel(pd, pN2, "6"))
        self.play(FadeIn(tr2))
        chk = T("ab=6, ac=3, ad=8, bc=7, bd=12, cd=7: all six match, NJ recovered the additive tree", 20, C_GOOD).move_to([0, -3.0, 0])
        chk.set_y(-2.85)
        self.say("Check all six distances: they all match the matrix. NJ recovered the additive tree.", FadeIn(chk))
        self.play(FadeOut(VGroup(tr1, tr2, f3, f4, t2, chk, mv)))

        # UPGMA comparison
        self.say("Now compare with UPGMA on the same matrix. It joins a and c, the closest pair at distance 3, then adds b and finally d.")
        ptsA = {"a": [-4.4, -2.2, 0], "c": [-3.2, -2.2, 0], "b": [-2.0, -2.2, 0], "d": [-0.8, -2.2, 0]}
        # UPGMA heights: ac 1.5, acb 3.25, root 4.5  (scale 0.5)
        sc = 0.5
        base = -2.2
        def yy(h): return base + h * sc
        up = VGroup()
        for n, p in ptsA.items():
            up.add(Text(n, font_size=24).move_to([p[0], base - 0.3, 0]))
        xac = (ptsA["a"][0] + ptsA["c"][0]) / 2
        xacb = (xac + ptsA["b"][0]) / 2
        xroot = (xacb + ptsA["d"][0]) / 2
        up.add(seg(ptsA["a"], [ptsA["a"][0], yy(1.5), 0]), seg(ptsA["c"], [ptsA["c"][0], yy(1.5), 0]), seg([ptsA["a"][0], yy(1.5), 0], [ptsA["c"][0], yy(1.5), 0]))
        up.add(seg([xac, yy(1.5), 0], [xac, yy(3.25), 0]), seg(ptsA["b"], [ptsA["b"][0], yy(3.25), 0]), seg([xac, yy(3.25), 0], [ptsA["b"][0], yy(3.25), 0]))
        up.add(seg([xacb, yy(3.25), 0], [xacb, yy(4.5), 0]), seg(ptsA["d"], [ptsA["d"][0], yy(4.5), 0]), seg([xacb, yy(4.5), 0], [ptsA["d"][0], yy(4.5), 0]))
        lab_up = T("UPGMA: (((a,c),b),d)", 24, C_BAD).move_to([-2.6, 0.9, 0])
        self.play(FadeIn(up), FadeIn(lab_up), run_time=1.5)
        nj = VGroup(seg([2.4, 0.2, 0], [3.4, -0.5, 0]), seg([2.4, -1.4, 0], [3.4, -0.5, 0]), seg([3.4, -0.5, 0], [4.9, -0.5, 0]), seg([4.9, -0.5, 0], [5.9, 0.2, 0]), seg([4.9, -0.5, 0], [5.9, -1.4, 0]))
        nj.add(Text("a", font_size=24).move_to([2.1, 0.2, 0]), Text("b", font_size=24).move_to([2.1, -1.4, 0]), Text("c", font_size=24).move_to([6.2, 0.2, 0]), Text("d", font_size=24).move_to([6.2, -1.4, 0]))
        lab_nj = T("NJ: ((a,b),(c,d))", 24, C_GOOD).move_to([4.2, 0.9, 0])
        self.set_subtitle("The topologies differ. The four point condition already told us a pairs with b, because the smallest sum was ab plus cd. NJ is right.")
        self.play(FadeIn(nj), FadeIn(lab_nj), run_time=1.5)
        self.hold()
        why = T("a and c are close only because both changed little (branches 1 and 1);\nUPGMA trusts the raw distance, NJ subtracts the separations u", 20, YELLOW).move_to([0, -1.75, 0])
        why.set_y(2.5)
        self.say("Why did UPGMA fail? The branches of b and d are long, those lineages changed a lot, while a and c changed little. So a and c are close without being siblings. UPGMA trusts the raw distance, and NJ discounts it with the separations.",
                 FadeIn(why))


# =====================================================================
class A5_Sandwich(NarratedScene):
    """Slides 105 to 134: the ultrametric sandwich."""

    PAIRS = {("A", "B"): (3, 7), ("A", "C"): (2, 3), ("A", "D"): (4, 6), ("A", "E"): (3, 5), ("B", "C"): (4, 10),
             ("B", "D"): (1, 4), ("B", "E"): (1, 8), ("C", "D"): (3, 8), ("C", "E"): (3, 5), ("D", "E"): (1, 7)}

    def construct(self):
        title = self.corner_title("The ultrametric sandwich")
        self.add(title)
        P = self.PAIRS
        names = list("ABCDE")
        self.say("Sometimes the honest input is an interval: the distance between i and j is at least Ml and at most Mh. We want an ultrametric tree that fits between the two matrices, like the filling of a sandwich.")
        rows = [["pair", "Ml", "Mh"]] + [[a + b, str(P[(a, b)][0]), str(P[(a, b)][1])] for (a, b) in P]
        tb = table(rows, col_w=0.85, row_h=0.32, font=18).move_to([-5.4, 1.2, 0])
        self.play(FadeIn(tb))
        self.say("Mh will decide the shape of the tree, and Ml will decide the heights. First build a minimum spanning tree of the graph that has the Mh weights, using Kruskal.")

        ang = {n: 90 - 72 * i for i, n in enumerate(names)}
        gpos = {n: np.array([-0.3, 0.2, 0]) + 1.6 * np.array([np.cos(ang[n] * DEGREES), np.sin(ang[n] * DEGREES), 0]) for n in names}
        gnodes = {n: tnode(n, gpos[n], r=0.26) for n in names}
        self.play(*[FadeIn(g) for g in gnodes.values()])
        order = sorted(P, key=lambda p: (P[p][1], p))
        parent = {n: n for n in names}

        def find(x):
            while parent[x] != x:
                x = parent[x]
            return x

        mst = []
        elist = VGroup()
        elab = T("", 20)
        tree_edges = {}
        for k, (a, b) in enumerate(order):
            w = P[(a, b)][1]
            ln = seg(gpos[a], gpos[b], YELLOW, 4)
            lb = T(f"{a}{b}: {w}", 22, YELLOW).move_to([3.4, 2.6, 0])
            self.set_subtitle(f"Next lightest edge: {a} {b} with weight {w}.")
            self.play(ShowCreation(ln), FadeIn(lb), run_time=0.4)
            if find(a) != find(b):
                parent[find(a)] = find(b)
                mst.append((a, b, w))
                self.play(ln.animate.set_color(C_GOOD), run_time=0.3)
                tree_edges[(a, b)] = ln
                self.play(FadeOut(lb), run_time=0.2)
            else:
                self.play(ln.animate.set_color(C_BAD), run_time=0.3)
                self.wait(0.3)
                self.play(FadeOut(ln), FadeOut(lb), run_time=0.3)
            if len(mst) == 4:
                break
        assert sorted((a + b) for a, b, _ in mst) == ["AC", "AD", "AE", "BD"], mst
        self.say("The tree T has the edges A C with 3, B D with 4, A E with 5, and A D with 6. The edge C E also weighs 5, but it would close the cycle A, C, E, so it is skipped.")

        # cut weights
        adj = {n: {} for n in names}
        for a, b, w in mst:
            adj[a][b] = (w, (a, b))
            adj[b][a] = (w, (a, b))

        def heaviest(s, t):
            stack = [(s, [])]
            seen = {s}
            while stack:
                v, es = stack.pop()
                if v == t:
                    return max(es, key=lambda e: (e[0], e[1]))
                for nb, (w, e) in adj[v].items():
                    if nb not in seen:
                        seen.add(nb)
                        stack.append((nb, es + [(w, e)]))

        cw = {}
        assigned = {}
        for (a, b), (lo, hi) in P.items():
            w, e = heaviest(a, b)
            assigned.setdefault(e, []).append(((a, b), lo))
            cw[e] = max(cw.get(e, 0), lo)
        self.play(FadeOut(tb))
        self.say("Now the cut weight. Each pair of objects is assigned to the heaviest edge on its path in T. All pairs of one edge are separated at the same node of the final tree, so they get the same distance, which must be at least the largest Ml among them.")
        cwt = []
        for e in sorted(cw, key=lambda e: (cw[e], e)):
            prs = ", ".join(f"{p[0]}{p[1]}" for p, _ in assigned[e])
            lows = ", ".join(str(l) for _, l in assigned[e])
            cwt.append([f"{e[0]}{e[1]}", prs, lows, str(cw[e])])
        tbl = table([["edge", "pairs", "Ml", "CW"]] + cwt, col_w=2.0, row_h=0.45, font=19)
        # custom widths: make the pairs column wider by pushing the last two columns right
        for r_ in range(len(cwt) + 1):
            tbl.cells[r_][2].shift(RIGHT * 0.9)
            tbl.cells[r_][3].shift(RIGHT * 0.9)
        tbl.move_to([-3.1, 0.4, 0])
        gshift = VGroup(*tree_edges.values(), *gnodes.values())
        self.play(FadeIn(tbl, lag_ratio=0.1), gshift.animate.shift(RIGHT * 3.9), run_time=1.5)
        self.hold()
        ex = T("example: C to E goes C, A, E with weights 3 and 5, so (C,E) is assigned to (A,E)", 20, GREY_B).move_to([-0.5, -2.45, 0])
        self.say("Take C and E. The path is C, A, E with weights 3 and 5, so the pair is assigned to the edge A E. And B to E goes B, D, A, E with weights 4, 6 and 5, so it is assigned to A D.", FadeIn(ex))
        self.play(FadeOut(tbl), FadeOut(ex), *[FadeOut(v) for v in tree_edges.values()], *[FadeOut(g) for g in gnodes.values()])

        # build the tree (UPGMA-like driven by CW)
        self.say("Now the algorithm looks like UPGMA driven by the cut weights. Each object starts as a tree of height 0. Process the edges of T by non-decreasing cut weight, joining the trees that contain the two ends, with height CW over 2.")
        leaf_order = ["B", "D", "A", "C", "E"]
        lx = {n: -2.0 + 1.1 * i for i, n in enumerate(leaf_order)}
        base_y, sc = -2.1, 0.7

        def yy(h):
            return base_y + h * sc

        comp = {n: n for n in names}
        members = {n: [n] for n in names}
        top = {n: (lx[n], Fraction(0)) for n in names}
        for n in names:
            self.add(Text(n, font_size=24).move_to([lx[n], base_y - 0.3, 0]))
        self.add(Line([-2.6, base_y, 0], [4.0, base_y, 0]).set_stroke(GREY_D, 2))
        steps_rows = [["edge", "CW", "height", "new edges"]]
        sorted_edges = sorted(cw, key=lambda e: (cw[e], e))
        for e in sorted_edges:
            a, b = e
            ra, rb = comp[a], comp[b]
            h = Fraction(cw[e], 2)
            (xa, ha), (xb, hb) = top[ra], top[rb]
            grp = VGroup(seg([xa, yy(ha), 0], [xa, yy(h), 0]), seg([xb, yy(hb), 0], [xb, yy(h), 0]), seg([xa, yy(h), 0], [xb, yy(h), 0]))
            w1 = T(fnum(h - ha), 18, YELLOW).move_to([xa - 0.25, (yy(ha) + yy(h)) / 2, 0])
            w2 = T(fnum(h - hb), 18, YELLOW).move_to([xb + 0.25, (yy(hb) + yy(h)) / 2, 0])
            info = T(f"edge {a}{b}: CW = {cw[e]}, height = {fnum(h)}", 24, C_GOOD).move_to([4.2, 1.6, 0])
            self.set_subtitle(f"Edge {a} {b} has cut weight {cw[e]}. Create a node of height {fnum(h)}. Its edges weigh {fnum(h - ha)} and {fnum(h - hb)}.")
            self.play(FadeIn(grp), FadeIn(w1), FadeIn(w2), FadeIn(info), run_time=1.0)
            self.hold()
            self.play(FadeOut(info), run_time=0.2)
            newx = (xa + xb) / 2
            for m in members[rb]:
                comp[m] = ra
            members[ra] += members[rb]
            top[ra] = (newx, h)

        # sandwich check
        # tree distance = 2 * height of the lowest common node; recompute from the construction
        def treedist(a, b):
            w, e = heaviest(a, b)
            return cw[e]
        lines = ["pair  Ml  tree  Mh"] + [f"{a}{b}     {lo}    {treedist(a, b)}     {hi}" for (a, b), (lo, hi) in P.items()]
        for (a, b), (lo, hi) in P.items():
            assert lo <= treedist(a, b) <= hi
        chk = T("every tree distance lies between its two limits: the sandwich exists", 24, C_GOOD).move_to([0.5, 3.0, 0])
        self.say("The distance between two leaves is twice the height of the node where they meet. Compare with the input: every tree distance, for example 4 for A to B between 3 and 7, lies between its two limits. The sandwich exists.", FadeIn(chk))
        cx = MathT(r"O(n^{2}\log n)", 34, GREY_B).move_to([4.8, -1.2, 0])
        self.say("With Kruskal for the spanning tree and a union-find structure for the sets, the whole construction is cheap: sorting the edges dominates.", FadeIn(cx))


# =====================================================================
def layout_tree(children, root, xgap=1.5, ygap=1.4, top=2.4, cx=0.0):
    """children: dict node -> list of children (ordered). Returns {node: position}."""
    xs = {}
    counter = [0]

    def place(v):
        ch = children.get(v, [])
        if not ch:
            xs[v] = counter[0]
            counter[0] += 1
        else:
            for c in ch:
                place(c)
            xs[v] = (xs[ch[0]] + xs[ch[-1]]) / 2

    place(root)
    depth = {root: 0}
    stack = [root]
    while stack:
        v = stack.pop()
        for c in children.get(v, []):
            depth[c] = depth[v] + 1
            stack.append(c)
    mid = (max(xs.values()) + min(xs.values())) / 2
    return {v: np.array([cx + (xs[v] - mid) * xgap, top - depth[v] * ygap, 0]) for v in xs}


class A5_PerfectPhylogeny(NarratedScene):
    """Slides 136 to 165: perfect phylogeny, the matrix L, and the construction."""

    M = {"a1": [1, 1, 0, 1, 0], "a2": [1, 1, 0, 0, 0], "a3": [1, 0, 0, 0, 0], "a4": [0, 0, 1, 0, 0], "a5": [0, 0, 1, 0, 1]}
    CH = ["c1", "c2", "c3", "c4", "c5"]

    def construct(self):
        title = self.corner_title("Perfect phylogeny and the matrix L")
        self.add(title)
        objs = list(self.M)
        self.say("Now the input is a table of binary characters. Each character is absent or present, evolves from 0 to 1, and once it reaches 1 it never goes back. So each character is gained exactly once in history.")
        rows = [[""] + self.CH] + [[o] + [str(v) for v in self.M[o]] for o in objs]
        Mt = table(rows, col_w=0.75, row_h=0.45, font=24).move_to([-4.8, 1.6, 0])
        Mlabel = T("M", 26, GREY_B).next_to(Mt, UP, buff=0.2)
        self.play(FadeIn(Mt), FadeIn(Mlabel))
        self.say("A perfect phylogeny is a rooted tree whose leaves are the objects, where each character labels exactly one edge, and for each object the labels on the path to the root are exactly its characters. It exists if and only if, for every two characters, the sets of objects are disjoint or one contains the other.")

        # sort columns
        counts = [sum(self.M[o][j] for o in objs) for j in range(5)]
        order = sorted(range(5), key=lambda j: (-counts[j], j))
        self.set_subtitle("Gusfield's test: first sort the columns by number of ones, the columns with more ones on the left. Here c1 has 3 ones, c2 and c3 have 2, c4 and c5 have 1, and the columns are already sorted.")
        cnt = T("ones per column: " + ", ".join(f"{self.CH[j]}:{counts[j]}" for j in range(5)), 22, YELLOW).move_to([-4.0, -1.1, 0])
        self.play(FadeIn(cnt))
        self.hold()
        self.play(FadeOut(cnt))

        # build L
        L = {}
        for o in objs:
            row = [self.M[o][j] for j in order]
            Lr = []
            last = None
            for pos_, v in enumerate(row, start=1):
                if v == 0:
                    Lr.append(0)
                else:
                    Lr.append(-1 if last is None else last)
                    last = pos_
            L[o] = Lr
        Lrows = [[""] + [self.CH[j] for j in order]] + [[o] + [str(v) for v in L[o]] for o in objs]
        Lt = table(Lrows, col_w=0.75, row_h=0.45, font=24).move_to([1.2, 1.6, 0])
        Llabel = T("L", 26, GREY_B).next_to(Lt, UP, buff=0.2)
        self.set_subtitle("Build the matrix L. Where M has a 0, L has a 0. Where M has a 1, L stores the position of the nearest 1 to its left in the same row, or minus 1 if there is none.")
        self.play(FadeIn(Llabel))
        for ri, o in enumerate(objs):
            self.play(FadeIn(VGroup(*Lt.cells[ri + 1])), run_time=0.7)
        self.play(FadeIn(VGroup(*Lt.cells[0])), run_time=0.3)
        self.hold()
        # check columns
        self.set_subtitle("Check each column of L. In a perfect phylogeny these are the characters on the path of each object, so the nearest 1 to the left is the parent of the character. A character has only one parent, so all nonzero values in a column must agree.")
        marks = VGroup()
        ok = True
        for j in range(5):
            vals = [L[o][j] for o in objs if L[o][j] != 0]
            col_ok = len(set(vals)) <= 1
            ok = ok and col_ok
            cells = [Lt.cells[ri + 1][j + 1] for ri, o in enumerate(objs) if L[o][j] != 0]
            if cells:
                r = SurroundingRectangle(VGroup(*cells), buff=0.08).set_stroke(C_GOOD if col_ok else C_BAD, 3)
                marks.add(r)
        self.play(ShowCreation(marks), run_time=1.5)
        self.hold()
        assert ok
        verdict = T("all nonzero values agree in every column: a perfect phylogeny exists", 24, C_GOOD).move_to([0.5, -1.3, 0])
        self.play(FadeIn(verdict))
        self.say("In every column the nonzero values are equal, so a perfect phylogeny exists.")
        self.play(FadeOut(VGroup(marks, verdict, Lt, Llabel)))

        # construction
        sets = {self.CH[j]: frozenset(o for o in objs if self.M[o][j]) for j in range(5)}
        self.set_subtitle("Now build the tree. Create the root, the set of all objects. For each character in sorted order, find the smallest superset already in the tree, hang a new vertex for its set, and label the edge with the character.")
        children = {"root": []}
        edge_label = {}
        nodes_set = {"root": frozenset(objs)}
        name_of = {frozenset(objs): "root"}
        cur = VGroup()

        def setname(s):
            return "{" + ",".join(sorted(s)) + "}"

        def draw(children, edge_label, labels, base=None):
            pos = layout_tree(children, "root", xgap=1.35, ygap=1.15, top=2.2, cx=3.0)
            g = VGroup()
            for p, cs in children.items():
                for c in cs:
                    g.add(seg(pos[p], pos[c], GREY_B, 3))
                    g.add(wlabel(pos[p], pos[c], edge_label[c], size=18, off=0.2))
            for v, p in pos.items():
                if v == "root":
                    g.add(Dot(p, radius=0.09, color=YELLOW))
                else:
                    g.add(Text(labels[v], font_size=16, color=C_QUERY).move_to(p))
            return g

        labels = {}
        cur = draw(children, edge_label, labels)
        self.add(cur)
        order_chars = [self.CH[j] for j in order]
        for ch in order_chars:
            s = sets[ch]
            # smallest superset already placed
            parent = min((v for v, vs in nodes_set.items() if s < vs), key=lambda v: len(nodes_set[v]))
            nid = ch
            children[parent].append(nid)
            children[nid] = []
            edge_label[nid] = ch
            nodes_set[nid] = s
            labels[nid] = setname(s) if len(s) > 1 else next(iter(s))
            new = draw(children, edge_label, labels)
            self.set_subtitle(f"Character {ch} has the set {setname(s)}. The smallest superset already in the tree is {setname(nodes_set[parent]) if parent != 'root' else 'the root'}, so the new vertex goes below it, labeled {ch}.")
            self.play(ReplacementTransform(cur, new), run_time=1.0)
            cur = new
            self.hold()
        # missing leaves
        for o in objs:
            if frozenset([o]) not in nodes_set.values():
                parent = min((v for v, vs in nodes_set.items() if o in vs), key=lambda v: len(nodes_set[v]))
                nid = o
                children[parent].append(nid)
                children[nid] = []
                edge_label[nid] = "λ"
                nodes_set[nid] = frozenset([o])
                labels[nid] = o
                new = draw(children, edge_label, labels)
                self.set_subtitle(f"The object {o} has no vertex of its own, so it becomes a leaf below the smallest set that contains it, with the empty label lambda.")
                self.play(ReplacementTransform(cur, new), run_time=1.0)
                cur = new
                self.hold()
        self.say("Check: a1 has c4, c2 and c1 on its path, and its row is 1 1 0 1 0. a5 has c5 and c3, and its row is 0 0 1 0 1. The tree is a perfect phylogeny, with exactly one change per character.")
        self.play(FadeOut(VGroup(cur, Mt, Mlabel)))

        # failing example
        M2 = {"a1": [1, 1, 0, 1, 0], "a2": [1, 1, 0, 0, 0], "a3": [1, 0, 1, 0, 0], "a4": [0, 0, 1, 0, 0], "a5": [0, 0, 1, 0, 1]}
        chs = ["c1", "c3", "c2", "c4", "c5"]
        cols = [0, 2, 1, 3, 4]
        Ms = {o: [M2[o][j] for j in cols] for o in objs}
        Lr2 = {}
        for o in objs:
            last = None
            r = []
            for pos_, v in enumerate(Ms[o], start=1):
                if v == 0:
                    r.append(0)
                else:
                    r.append(-1 if last is None else last)
                    last = pos_
            Lr2[o] = r
        A = VGroup(table([[""] + chs] + [[o] + [str(v) for v in Ms[o]] for o in objs], col_w=0.75, row_h=0.45, font=22).move_to([-4.2, 1.2, 0]))
        B = VGroup(table([[""] + chs] + [[o] + [str(v) for v in Lr2[o]] for o in objs], col_w=0.75, row_h=0.45, font=22).move_to([1.6, 1.2, 0]))
        self.set_subtitle("Now change one cell: a3 also has c3. After sorting the columns as c1, c3, c2, c4, c5, look at the column of c3 in L.")
        self.play(FadeIn(A), FadeIn(B))
        bt = B[0]
        col = [bt.cells[ri + 1][2] for ri, o in enumerate(objs) if Lr2[o][1] != 0]
        bad = SurroundingRectangle(VGroup(*col), buff=0.08).set_stroke(C_BAD, 3)
        self.play(ShowCreation(bad))
        self.hold()
        msg = T("column c3 holds 1, -1, -1: they disagree, no perfect phylogeny", 24, C_BAD).move_to([0, -1.4, 0])
        self.say("The nonzero values are 1, minus 1 and minus 1. They are not all equal. The incompatible sets are A1 equal to a1, a2, a3 and A3 equal to a3, a4, a5: they share a3 and neither contains the other. A single origin of c3 cannot satisfy a3, which says it appeared below c1, and a4 and a5, which say it appeared below the root.",
                 FadeIn(msg))
        cx = MathT(r"O(nm)", 40, GREY_B).move_to([0, -2.4, 0])
        self.say("The whole test is order n m, the size of the input: counting ones, sorting columns with counting sort, building L, and checking its columns.", FadeIn(cx))


# =====================================================================
class A5_Fitch(NarratedScene):
    """Slide 166: small parsimony with Fitch's algorithm (a worked example that is not on the slides)."""

    def construct(self):
        title = self.corner_title("Small parsimony: Fitch")
        self.add(title)
        self.say("Real data rarely allow a perfect phylogeny, so we use parsimony: the best tree explains the data with the fewest state changes. Small parsimony gives the tree and asks for the minimum number of changes.")
        note = T("(example built for this animation, not taken from the slides)", 20, GREY_B).move_to([0, 3.1, 0])
        self.add(note)
        children = {"r": ["X", "Y"], "X": ["l1", "l2"], "Y": ["l3", "W"], "W": ["l4", "l5"]}
        leaf = {"l1": "A", "l2": "C", "l3": "C", "l4": "G", "l5": "G"}
        pos = layout_tree(children, "r", xgap=1.7, ygap=1.35, top=2.0, cx=0.0)
        lines = VGroup()
        for p, cs in children.items():
            for c in cs:
                lines.add(seg(pos[p], pos[c], GREY_B, 3))
        self.play(ShowCreation(lines))
        leaf_m = {}
        for k, s in leaf.items():
            m = VGroup(Circle(radius=0.3).set_stroke(C_QUERY, 3).set_fill(BLACK, 1), Text(s, font_size=26)).move_to(pos[k])
            leaf_m[k] = m
        self.play(*[FadeIn(m) for m in leaf_m.values()])

        sets = {k: {v} for k, v in leaf.items()}
        cost = 0
        shown = {}
        self.set_subtitle("One pass from the leaves to the root. At each internal node take the intersection of the children's sets if it is not empty. Otherwise take the union and count one change.")
        post = ["X", "W", "Y", "r"]
        for v in post:
            a, b = children[v]
            inter = sets[a] & sets[b]
            if inter:
                sets[v] = inter
                note_txt = f"{{{''.join(sorted(sets[a]))}}} and {{{''.join(sorted(sets[b]))}}} intersect: keep {{{''.join(sorted(inter))}}}"
            else:
                sets[v] = sets[a] | sets[b]
                cost += 1
                note_txt = f"{{{''.join(sorted(sets[a]))}}} and {{{''.join(sorted(sets[b]))}}} are disjoint: union, one change"
            m = VGroup(Circle(radius=0.42).set_stroke(YELLOW if not inter else C_GOOD, 3).set_fill(BLACK, 1),
                       Text("{" + "".join(sorted(sets[v])) + "}", font_size=22)).move_to(pos[v])
            msg = T(note_txt, 24, WHITE).move_to([4.6, 2.4 - 0.6 * post.index(v), 0])
            if msg.get_width() > 4.6:
                msg.set_width(4.6)
            self.play(FadeIn(m), FadeIn(msg), run_time=0.9)
            self.wait(0.9)
            shown[v] = m
        assert cost == 2
        total = T(f"minimum changes = {cost}", 26, C_GOOD).move_to([5.2, -0.6, 0])
        self.say("The number of unions is the minimum number of changes, here 2. This takes one pass over the tree, order n k for k states. Sankoff handles different costs for each change in order n k squared.", FadeIn(total))
        self.say("Large parsimony, when the tree is unknown, is NP-complete. Programs use local search, moving to a neighbor tree with NNI or SPR and keeping the moves that reduce the mutations. A perfect phylogeny is the ideal case: exactly one change per character.")
