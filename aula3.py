"""Aula 3: BLAST.

Render one scene:   uv run manimgl aula3.py A3_SeedList -w
Render them all:    ./render_all.sh aula3
Scenes: A3_Overview, A3_SeedList, A3_FindHits, A3_ExtendHSP, A3_TwoHit,
        A3_GappedBLAST, A3_Significance
"""
import random

from common import *

QUERY = "VHREMAARTSPLRPLVATAGPALSPVPPCVHLTLR"   # slide 14, 35 residues


# =====================================================================
class A3_Overview(NarratedScene):
    """Slide 10 (why exact alignment fails) and slides 11 and 12 (BLAST in one idea)."""

    def construct(self):
        title = self.corner_title("BLAST: why a heuristic?")
        self.add(title)

        self.say("Smith-Waterman fills one cell for every pair of a query position and a database position. That is n times M cells.")
        grid = VGroup(*[Square(0.28).set_stroke(GREY_B, 1) for _ in range(30)])
        grid.arrange_in_grid(5, 6, buff=0)
        grid.shift(3.8 * LEFT + 0.4 * UP)
        nlab = MathT("n", 34).next_to(grid, UP, buff=0.15)
        mlab = MathT("M", 34).next_to(grid, LEFT, buff=0.2)
        self.play(ShowCreation(grid), Write(nlab), Write(mlab))

        formula = MathT(r"\text{time}=\frac{n\times 10^{13}}{3\times 2^{30}}\ \text{s}", 34)
        formula.shift(3.3 * RIGHT + 2.5 * UP)
        self.say("Take M = 10^13 bases, one instruction per cell, and a 3 GHz machine. Now vary the query size n.",
                 Write(formula), run_time=1.5)

        rows = [("100", "about 4 days"), ("1,000", "about 36 days"), ("10,000", "about 1 year"),
                ("100,000", "about 1 decade"), ("1,000,000", "about 1 century"),
                ("10,000,000", "about 1 millennium")]
        table = VGroup()
        for n, t in rows:
            a = T(n, 26, color=C_QUERY)
            b = T(t, 26, color=C_HI)
            table.add(VGroup(a, b))
        for k, r in enumerate(table):
            r[0].move_to([1.9, 1.3 - 0.52 * k, 0], aligned_edge=RIGHT)
            r[1].move_to([2.4, 1.3 - 0.52 * k, 0], aligned_edge=LEFT)
        head = VGroup(T("n", 26, GREY_B).move_to([1.9, 1.85, 0], aligned_edge=RIGHT),
                      T("time", 26, GREY_B).move_to([2.4, 1.85, 0], aligned_edge=LEFT))
        self.set_subtitle("Each line is ten times the previous one. A single exact search takes days to years, and the database doubles every 18 months.")
        self.play(FadeIn(head))
        for r in table:
            self.play(FadeIn(r, shift=0.2 * RIGHT), run_time=0.5)
        self.hold()

        self.play(FadeOut(VGroup(grid, nlab, mlab, formula, table, head)))

        # BLAST in one idea
        idea = T("Good local alignments contain a short, almost identical stretch", 30, YELLOW)
        idea.shift(2.4 * UP)
        self.say("Here is the idea behind BLAST. A good local alignment almost always contains a short stretch that is identical or nearly identical.",
                 Write(idea), run_time=1.5)

        steps = VGroup(
            VGroup(T("1", 40, YELLOW), T("build a list\nof seeds\nfrom the query", 26)),
            VGroup(T("2", 40, YELLOW), T("find hits of\nthe seeds in\nthe database", 26)),
            VGroup(T("3", 40, YELLOW), T("extend each hit\nto obtain HSPs", 26)),
        )
        boxes = VGroup()
        for s in steps:
            s[1].next_to(s[0], DOWN, buff=0.2)
            b = SurroundingRectangle(s, buff=0.3).set_stroke(C_SEED, 2)
            boxes.add(VGroup(b, s))
        boxes.arrange(RIGHT, buff=0.8).shift(0.2 * DOWN)
        arrows = VGroup(*[Arrow(boxes[i].get_right(), boxes[i + 1].get_left(), buff=0.1) for i in range(2)])
        self.set_subtitle("Finding short stretches is cheap, so BLAST looks for them first and only spends alignment effort around the places where it found one. These are the three basic steps.")
        for i, b in enumerate(boxes):
            self.play(FadeIn(b, shift=0.3 * UP), *( [ShowCreation(arrows[i - 1])] if i else []), run_time=0.8)
        self.hold()

        vocab = T("query   database   seed   hit   HSP (high-scoring segment pair)", 24, GREY_B)
        vocab.next_to(boxes, DOWN, buff=0.6)
        self.say("Five words to learn: query, database, seed, hit, and HSP, the pair of similar segments that BLAST finally reports.",
                 FadeIn(vocab))


# =====================================================================
class A3_SeedList(NarratedScene):
    """Slides 13 to 29: words of the query and neighbor words above the threshold T."""

    def construct(self):
        title = self.corner_title("Step 1: the seed list")
        self.add(title)

        self.say("Seeds are words of length w. In the slides w is 4 and the query has 35 residues.")
        q = row_of_cells(list(QUERY), size=0.385, font=17)
        q.move_to(2.6 * UP)
        self.play(FadeIn(q, lag_ratio=0.02), run_time=1.5)

        win = SurroundingRectangle(VGroup(*q.cells[0:4]), buff=0.0).set_stroke(YELLOW, 4)
        self.add(win)
        self.set_subtitle("A window of size 4 slides one position at a time: VHRE, HREM, REMA, EMAA, MAAR, and so on. That is the staircase you see on the slides.")
        words = VGroup()
        for k in range(8):
            if k:
                self.play(win.animate.move_to(VGroup(*q.cells[k:k + 4])), run_time=0.35)
            w = T(QUERY[k:k + 4], 24, YELLOW)
            w.move_to(q.cells[k].get_center() + 0.55 * DOWN * (k + 1) * 0.55 + 0.25 * DOWN + 0.55 * RIGHT)
            words.add(w)
            self.play(FadeIn(w, shift=0.15 * DOWN), run_time=0.3)
        self.hold()

        count = MathT(r"n-w+1=35-4+1=32\ \text{words}", 32, C_HI)
        count.to_edge(RIGHT, buff=0.8).shift(0.4 * DOWN)
        self.say("A sequence of length n has n minus w plus 1 words. Here that is 32.", Write(count))
        self.play(FadeOut(win), FadeOut(words), FadeOut(count))

        # neighbors of EMAA
        self.play(q.animate.scale(0.75).to_edge(UP, buff=0.9))
        em = T("EMAA", 44, YELLOW)
        em.move_to(3.6 * LEFT + 0.6 * UP)
        self.say("But the seeds are not only the words that appear in the query. They are also the neighbors, any word whose score against a query word is at least T.",
                 FadeIn(em))
        self.say("Pick the word EMAA and score candidate words against it, position by position, with no gaps, adding the substitution matrix values.")

        rows = [("EMAA", 18), ("EHAI", 14), ("EHAK", 9), ("AAAC", 6), ("AAAD", 3), ("YYYY", -15)]
        T_thr = 11
        hdr = VGroup(T("word", 24, GREY_B), T("score", 24, GREY_B), T("score >= T ?", 24, GREY_B))
        hdr[0].move_to([-0.5, 1.6, 0]); hdr[1].move_to([1.6, 1.6, 0]); hdr[2].move_to([4.2, 1.6, 0])
        self.play(FadeIn(hdr))
        tline = T("T = 11", 34, C_HI).move_to([4.3, 2.4, 0])
        self.play(Write(tline))
        shown = VGroup()
        for k, (wd, sc) in enumerate(rows):
            y = 1.0 - 0.55 * k
            ok = sc >= T_thr
            col = C_GOOD if ok else C_BAD
            a = T(wd, 28, WHITE).move_to([-0.5, y, 0])
            b = T(str(sc), 28, col).move_to([1.6, y, 0])
            c = T("yes, a seed" if ok else "no", 26, col).move_to([4.2, y, 0])
            self.play(FadeIn(a), FadeIn(b), run_time=0.4)
            self.play(FadeIn(c), run_time=0.3)
            shown.add(a, b, c)
        self.hold()

        note = T("EHAI is not in the query, yet it is a seed:  M->H and A->I cost only 4 points", 22, YELLOW)
        note.move_to([1.0, -2.55, 0])
        self.say("Notice EHAI. It does not occur in the query, but replacing M by H and A by I costs only a few points, so a database protein with EHAI there is a plausible relative. Searching only for EMAA would miss it.",
                 FadeIn(note))
        self.play(FadeOut(VGroup(shown, hdr, tline, em, note, q)))

        # role of T
        up = T("higher T: fewer neighbors, faster, less sensitive", 30, C_BAD)
        dn = T("lower T: more neighbors, slower, more sensitive", 30, C_GOOD)
        up.shift(0.8 * UP); dn.shift(0.4 * DOWN)
        self.say("The role of T is a trade-off. Higher T keeps fewer neighbors, so the search is faster and less sensitive. Lower T keeps more, so it is slower and more sensitive. BLAST 2.0 will play with this.",
                 FadeIn(up), FadeIn(dn), run_time=1.2)
        tune = T("w and T are tuned so that the seed list is about 50 times larger than the query", 24, GREY_B)
        tune.shift(1.8 * DOWN)
        self.say("The parameters w and T are tuned so that the seed list is about 50 times the size of the query.", FadeIn(tune))


# =====================================================================
class A3_FindHits(NarratedScene):
    """Slide 30: one pass over each database sequence, lookup of every word in the seed table."""

    def construct(self):
        title = self.corner_title("Step 2: finding hits")
        self.add(title)
        DB = "KWKEHAIQADPIAS"
        seeds = {"VHRE": [0], "HREM": [1], "REMA": [2], "EMAA": [3], "EHAI": [3], "MAAR": [4], "AART": [5],
                 "ARTS": [6], "RTSP": [7], "TSPL": [8], "SPLR": [9]}

        self.say("We have the seed list. Now we must find every place in the database where a seed occurs. We want one single pass over each database sequence.")
        # table on the left
        t_title = T("seed table (word -> query positions)", 24, GREY_B).move_to([-4.2, 2.7, 0])
        items = VGroup()
        for k, (w, ps) in enumerate(seeds.items()):
            col = C_SEED if w == "EHAI" else WHITE
            items.add(T(f"{w}  ->  query {', '.join(map(str, ps))}", 22, col))
        items.arrange(DOWN, aligned_edge=LEFT, buff=0.12).next_to(t_title, DOWN, buff=0.25)
        self.play(FadeIn(t_title), FadeIn(items, lag_ratio=0.1), run_time=1.5)

        dbrow = row_of_cells(list(DB), size=0.5, font=24, color=C_DB)
        dbrow.move_to([2.5, 1.2, 0])
        lab = T("database sequence", 22, GREY_B).next_to(dbrow, UP, buff=0.2)
        self.play(FadeIn(dbrow), FadeIn(lab))

        win = SurroundingRectangle(VGroup(*dbrow.cells[0:4]), buff=0).set_stroke(YELLOW, 4)
        cur = T("", 30)
        self.add(win)
        self.set_subtitle("A table has one position for every possible word. We take the word under the window, jump to its position, and read all hits in constant time.")
        hits = []
        for k in range(len(DB) - 3):
            word = DB[k:k + 4]
            if k:
                self.play(win.animate.move_to(VGroup(*dbrow.cells[k:k + 4])), run_time=0.4)
            newcur = T(word, 34, YELLOW).move_to([2.5, 0.0, 0])
            self.play(Transform(cur, newcur) if k else FadeIn(cur.become(newcur)), run_time=0.15)
            found = word in seeds
            if found:
                idx = list(seeds).index(word)
                hl = SurroundingRectangle(items[idx], buff=0.05).set_stroke(C_GOOD, 3)
                flash = SurroundingRectangle(VGroup(*dbrow.cells[k:k + 4]), buff=0.02).set_stroke(C_GOOD, 5)
                msg = T(f"HIT: database position {k} matches query position {seeds[word][0]}", 26, C_GOOD)
                msg.move_to([2.5, -1.0, 0])
                self.play(ShowCreation(hl), ShowCreation(flash), FadeIn(msg), run_time=0.8)
                self.wait(1.2)
                hits.append((word, k))
                self.play(FadeOut(hl), FadeOut(msg))
                self.add(flash)
            else:
                miss = T("no entry", 22, GREY_B).move_to([2.5, -1.0, 0])
                self.play(FadeIn(miss), run_time=0.1)
                self.play(FadeOut(miss), run_time=0.1)
        self.hold()

        auto = T("In practice BLAST runs a finite automaton instead:\none transition per letter, seed states announce a hit.", 24, YELLOW)
        auto.move_to([2.5, -1.8, 0])
        self.say("The table has 20 to the 4, that is 160,000 positions, and few are used, so it is stored as a hash table. BLAST in practice uses a finite automaton, one transition per database letter.",
                 FadeIn(auto), run_time=1.0)
        cost = MathT(r"\text{cost of this step: one pass over the database}= O(M)", 30, C_HI)
        cost.move_to([2.5, -2.7, 0])
        self.say("The cost of this step is a single pass over the database, order M.", FadeIn(cost))


# =====================================================================
class A3_ExtendHSP(NarratedScene):
    """Slides 31 to 56: gapless extension in both directions, with drop-off."""

    def construct(self):
        title = self.corner_title("Step 3: from a hit to an HSP")
        self.add(title)
        Q = "VHREMAARTSPLRP"
        D = "KWKEHAIQADPIAS"
        seed_i = (3, 4, 5, 6)
        seed_score = 14
        left = [("R", "K", 2), ("H", "W", -2), ("V", "K", -2)]
        right = [("R", "Q", 1), ("T", "A", 0), ("S", "D", 0), ("P", "P", 7), ("L", "I", 2), ("R", "A", -1), ("P", "S", -1)]
        col_score = {}
        for (qa, da, s), idx in zip(left, (2, 1, 0)):
            col_score[idx] = s
        for (qa, da, s), idx in zip(right, range(7, 14)):
            col_score[idx] = s

        self.say("Each hit is extended in both directions without gaps. In the dynamic programming matrix this means walking along one diagonal.")
        size = 0.62
        qrow = row_of_cells(list(Q), size=size, font=26, color=C_QUERY)
        drow = row_of_cells(list(D), size=size, font=26, color=C_DB)
        qrow.move_to([0, 2.5, 0]); drow.next_to(qrow, DOWN, buff=0)
        ql = T("query", 22, C_QUERY).next_to(qrow, LEFT, buff=0.3)
        dl = T("database", 22, C_DB).next_to(drow, LEFT, buff=0.3)
        self.play(FadeIn(VGroup(qrow, drow, ql, dl)))
        self.set_subtitle("The seed is EHAI in the database against EMAA in the query. Remember, EHAI scored 14 against EMAA, so it was kept as a neighbor word.")
        seedbox = SurroundingRectangle(VGroup(*qrow.cells[3:7], *drow.cells[3:7]), buff=0.03).set_stroke(C_SEED, 4)
        self.play(ShowCreation(seedbox))
        self.hold()

        # score rows
        sc_label = T("column score", 20, GREY_B).move_to([qrow.get_left()[0] - 1.0, 1.1, 0])
        run_label = T("running score", 20, GREY_B).move_to([qrow.get_left()[0] - 1.0, 0.6, 0])
        best_label = T("best so far", 20, GREY_B).move_to([qrow.get_left()[0] - 1.0, 0.1, 0])
        self.play(FadeIn(VGroup(sc_label, run_label, best_label)))
        xs = [c.get_center()[0] for c in qrow.cells]
        seed_txt = T("14", 26, C_SEED).move_to([(xs[3] + xs[6]) / 2, 1.1, 0])
        run_txt = T("14", 26, WHITE).move_to([(xs[3] + xs[6]) / 2, 0.6, 0])
        best_txt = T("14", 26, C_GOOD).move_to([(xs[3] + xs[6]) / 2, 0.1, 0])
        self.play(FadeIn(VGroup(seed_txt, run_txt, best_txt)))

        # mini plot of the running score
        base_y = -2.35
        sc = 0.05
        axis = Line([xs[0] - 0.4, base_y, 0], [xs[-1] + 0.4, base_y, 0]).set_stroke(GREY_B, 1)
        self.add(axis)
        pts = {}

        def px(i, v):
            return np.array([xs[i], base_y + v * sc, 0])

        marks = VGroup()
        dropline = DashedLine([xs[0] - 0.4, base_y, 0], [xs[-1] + 0.4, base_y, 0], color=C_BAD)
        # seed point
        mid = 4.5
        d_seed = Dot(np.array([(xs[3] + xs[6]) / 2, base_y + 14 * sc, 0]), radius=0.06, color=C_SEED)
        self.add(d_seed)
        running = {}
        # --- extend left
        self.set_subtitle("First to the left, one column at a time. R against K adds 2, so the running score becomes 16 and a new best. Then H against W and V against K both cost 2.")
        run_val = seed_score
        best = seed_score
        best_idx_left = 3
        prev_pt = d_seed.get_center()
        rows_left = []
        for idx in (2, 1, 0):
            s = col_score[idx]
            run_val += s
            box = SurroundingRectangle(VGroup(qrow.cells[idx], drow.cells[idx]), buff=0.03).set_stroke(YELLOW, 4)
            ct = T(f"{s:+d}", 26, C_GOOD if s > 0 else (WHITE if s == 0 else C_BAD)).move_to([xs[idx], 1.1, 0])
            rt = T(str(run_val), 26).move_to([xs[idx], 0.6, 0])
            if run_val > best:
                best = run_val
                best_idx_left = idx
            bt = T(str(best), 26, C_GOOD).move_to([xs[idx], 0.1, 0])
            dot = Dot(px(idx, run_val), radius=0.06, color=YELLOW)
            seg = Line(prev_pt, dot.get_center()).set_stroke(YELLOW, 3)
            prev_pt = dot.get_center()
            self.play(ShowCreation(box), FadeIn(ct), FadeIn(rt), FadeIn(bt), ShowCreation(seg), FadeIn(dot), run_time=0.9)
            self.wait(0.5)
            self.play(FadeOut(box), run_time=0.2)
        lbound = best_idx_left
        lend = SurroundingRectangle(VGroup(qrow.cells[lbound], drow.cells[lbound]), buff=0.03).set_stroke(C_GOOD, 4)
        self.set_subtitle("Going left only got worse after R against K. So the left end of the HSP is fixed there, with score 16.")
        self.play(ShowCreation(lend))
        self.hold()

        # --- extend right
        self.set_subtitle("Now to the right, starting from 16. The extension does not stop at the zeros. It walks through T/A and S/D and is rewarded with P/P, plus 7.")
        run_val = 16
        best = 16
        best_idx_right = 6
        prev_pt = px(lbound, 16)
        # restart plot line from the best-left point
        for idx in range(7, 14):
            s = col_score[idx]
            run_val += s
            if run_val > best:
                best = run_val
                best_idx_right = idx
            box = SurroundingRectangle(VGroup(qrow.cells[idx], drow.cells[idx]), buff=0.03).set_stroke(YELLOW, 4)
            ct = T(f"{s:+d}" if s else "0", 26, C_GOOD if s > 0 else (WHITE if s == 0 else C_BAD)).move_to([xs[idx], 1.1, 0])
            rt = T(str(run_val), 26).move_to([xs[idx], 0.6, 0])
            bt = T(str(best), 26, C_GOOD).move_to([xs[idx], 0.1, 0])
            dot = Dot(px(idx, run_val), radius=0.06, color=ORANGE)
            seg = Line(prev_pt if idx > 7 else px(6, 16), dot.get_center()).set_stroke(ORANGE, 3)
            prev_pt = dot.get_center()
            self.play(ShowCreation(box), FadeIn(ct), FadeIn(rt), FadeIn(bt), ShowCreation(seg), FadeIn(dot), run_time=0.8)
            self.wait(0.35)
            self.play(FadeOut(box), run_time=0.2)
        self.hold()

        # drop-off remark and final HSP
        hsp_cells = VGroup(*[VGroup(qrow.cells[i], drow.cells[i]) for i in range(lbound, best_idx_right + 1)])
        hsp = SurroundingRectangle(hsp_cells, buff=0.06).set_stroke(C_GOOD, 6)
        final = T("HSP score = 26", 34, C_GOOD).move_to([5.0, 0.55, 0])
        self.say("The reported HSP ends where the score was maximum, not where the extension stopped. The last two columns, minus 1 and minus 1, are thrown away. The HSP scores 26.",
                 ShowCreation(hsp), FadeIn(final), run_time=1.2)
        rule = T("stop when running score < best - 20 (proteins), report the best point", 24, YELLOW)
        rule.move_to([0, -2.85, 0])
        self.say("The drop limit is 20 for proteins. It exists so that the extension can be patient through zero columns, and with this value the chance of missing a better alignment is about 0.1 percent. Note there are no gaps anywhere.",
                 FadeIn(rule))


# =====================================================================
class A3_TwoHit(NarratedScene):
    """Slides 58 to 62: the two-hit method."""

    A_WIN = 20
    W = 3

    def construct(self):
        title = self.corner_title("BLAST 2.0: the two-hit method")
        self.add(title)
        rng = random.Random(11)
        XM, YM = 100, 60
        sx, sy = 0.07, 0.062
        ox, oy = -6.6, -2.2   # origin of the dot plot

        def P(x, y):
            return np.array([ox + x * sx, oy + y * sy, 0])

        frame = Rectangle(width=XM * sx, height=YM * sy).set_stroke(GREY_B, 2)
        frame.move_to(P(XM / 2, YM / 2))
        xl = T("query position", 20, GREY_B).next_to(frame, DOWN, buff=0.1)
        yl = T("database position", 20, GREY_B).rotate(PI / 2).next_to(frame, LEFT, buff=0.1)

        def diag(h):
            return h[0] - h[1]

        def conflict(h, others, win):
            return any(diag(h) == diag(o) and abs(h[0] - o[0]) <= win for o in others)

        planted = [((14, 8), (27, 21)), ((52, 30), (66, 44)), ((70, 3), (85, 18))]
        t13, t11 = [], []
        for a, b in planted:
            t13.append(a); t11.append(b)
        allh = list(t13) + list(t11)
        while len(t13) < 34:
            h = (rng.randint(2, 97), rng.randint(2, 57))
            if not conflict(h, allh, 24) and diag(h) not in [diag(o) for o in allh]:
                t13.append(h); allh.append(h)
        n13 = len(t13)
        while len(t11) < 3 + 82:
            h = (rng.randint(2, 97), rng.randint(2, 57))
            if not conflict(h, allh, 24) and diag(h) not in [diag(o) for o in allh]:
                t11.append(h); allh.append(h)

        def sq(h, color):
            s = Square(0.1).set_fill(color, 1).set_stroke(width=0)
            s.move_to(P(*h))
            return s

        self.say("Picture the query on one axis and a database sequence on the other. Every hit is a small square.")
        self.play(ShowCreation(frame), FadeIn(xl), FadeIn(yl))
        reds = VGroup(*[sq(h, C_BAD) for h in t13])
        self.say("These are the hits with T at least 13. The original one-hit method would extend every single one of them.",
                 FadeIn(reds, lag_ratio=0.05), run_time=1.5)
        c1 = T(f"one-hit: {n13} hits, {n13} extensions", 26, C_BAD).move_to([1.2, 2.2, 0], aligned_edge=LEFT)
        self.play(FadeIn(c1))

        blues = VGroup(*[sq(h, BLUE_C) for h in t11])
        self.say("Now lower T to 11. The blue squares are the new hits. There are clearly many more of them, so we will have many more chances to find a pair.",
                 FadeIn(blues, lag_ratio=0.03), run_time=2.0)
        n11 = n13 + len(t11)
        c2 = T(f"T = 11: {n11} hits (about 3 times more)", 26, BLUE_C).move_to([1.2, 1.5, 0], aligned_edge=LEFT)
        self.play(FadeIn(c2))

        # two-hit rule
        everything = t13 + t11
        pairs = []
        for i, a in enumerate(everything):
            for b in everything[i + 1:]:
                if diag(a) == diag(b):
                    lo, hi = (a, b) if a[0] < b[0] else (b, a)
                    if hi[0] - lo[0] >= self.W and hi[0] - lo[0] <= self.A_WIN:
                        pairs.append((lo, hi))
        assert len(pairs) == 3, len(pairs)
        self.set_subtitle("The rule: only start an extension when two hits sit on the same diagonal, closer than a distance A. Two overlapping hits count as one.")
        # draw diagonals through each pair
        marks = VGroup()
        for lo, hi in pairs:
            ln = Line(P(*lo), P(*hi)).set_stroke(WHITE, 4)
            ring = SurroundingRectangle(VGroup(sq(lo, WHITE), sq(hi, WHITE)), buff=0.12).set_stroke(YELLOW, 3)
            marks.add(VGroup(ln, ring))
        # fade everything else
        others = VGroup()
        for m in list(reds) + list(blues):
            others.add(m)
        self.play(others.animate.set_opacity(0.18), run_time=1.0)
        for m in marks:
            self.play(ShowCreation(m), run_time=0.9)
        self.hold()
        c3 = T("two-hit: 3 extensions", 26, YELLOW).move_to([1.2, 0.8, 0], aligned_edge=LEFT)
        self.say("Here only three places have two hits on the same diagonal within the window, so only those three are extended. All the other hits are discarded with no extension at all.",
                 FadeIn(c3), run_time=1.0)
        trade = T("finding a hit = table lookup (cheap)\nan extension = real work (expensive)", 24, GREY_B).move_to([1.2, -0.4, 0], aligned_edge=LEFT)
        self.play(FadeIn(trade))
        self.say("With the default protein parameters, about 3.2 times more hits, but about 86 percent fewer extensions. If one-hit found 100 hits and 100 extensions, two-hit finds 320 hits and only 14 extensions. Many cheap operations traded for few expensive ones.",
                 extra=0.5)

    @property
    def A(self):
        return self.A_WIN


# =====================================================================
class A3_GappedBLAST(NarratedScene):
    """Slides 63 and 64: gapped extension by dynamic programming with an X_G drop-off."""

    def construct(self):
        title = self.corner_title("Gapped BLAST: dynamic programming around a hit")
        self.add(title)
        q = "GATTACAGGCTA"
        d = "GATTACCAGGCTA"
        M_, MM_, G_, X = 2, -2, -3, 6
        n, m = len(q), len(d)

        self.say("In the original BLAST the extension lives on one diagonal. In gapped BLAST the extension is done with dynamic programming around the hit, so the alignment can contain gaps.")

        sz = 0.37
        ox, oy = -6.0, 2.3
        cells = {}
        grid = VGroup()
        for i in range(n + 1):
            for j in range(m + 1):
                sq = Square(sz).set_stroke(GREY_D, 1)
                sq.move_to([ox + j * sz, oy - i * sz, 0])
                cells[(i, j)] = sq
                grid.add(sq)
        qlabels = VGroup(*[T(c, 15, C_QUERY).move_to([ox - sz, oy - (i + 1) * sz, 0]) for i, c in enumerate(q)])
        dlabels = VGroup(*[T(c, 15, C_DB).move_to([ox + (j + 1) * sz, oy + sz, 0]) for j, c in enumerate(d)])
        self.play(FadeIn(grid), FadeIn(qlabels), FadeIn(dlabels), run_time=1.2)
        seed = T("seed", 20, C_SEED).move_to(cells[(0, 0)].get_center() + 0.5 * UP + 0.35 * LEFT)
        self.play(cells[(0, 0)].animate.set_fill(C_SEED, 0.8))
        self.set_subtitle(f"We start at the seed in the corner. Scores: match +{M_}, mismatch {MM_}, gap {G_}. The extension stops exploring a cell when its score falls more than X below the best score found, here X = {X}.")
        self.hold()

        S = {(0, 0): 0}
        kept = {(0, 0)}
        pruned = set()
        best, bestc = 0, (0, 0)
        order = []
        for k in range(1, n + m + 1):
            for i in range(0, n + 1):
                j = k - i
                if j < 0 or j > m:
                    continue
                cands = []
                if i > 0 and j > 0 and (i - 1, j - 1) in kept:
                    cands.append(S[(i - 1, j - 1)] + (M_ if q[i - 1] == d[j - 1] else MM_))
                if i > 0 and (i - 1, j) in kept:
                    cands.append(S[(i - 1, j)] + G_)
                if j > 0 and (i, j - 1) in kept:
                    cands.append(S[(i, j - 1)] + G_)
                if not cands:
                    continue
                v = max(cands)
                S[(i, j)] = v
                if v >= best - X:
                    kept.add((i, j))
                    if v > best:
                        best, bestc = v, (i, j)
                    order.append(((i, j), True))
                else:
                    pruned.add((i, j))
                    order.append(((i, j), False))

        # animate anti-diagonal by anti-diagonal
        self.set_subtitle("We fill the matrix one anti-diagonal at a time. Gray cells are computed and kept. Red crosses are cells pruned because they fell too far below the best score.")
        bydiag = {}
        for (c, ok) in order:
            bydiag.setdefault(c[0] + c[1], []).append((c, ok))
        numbers = VGroup()
        for k in sorted(bydiag):
            anims = []
            for (c, ok) in bydiag[k]:
                if ok:
                    anims.append(cells[c].animate.set_fill(GREY_B, 0.35 + 0.02 * max(S[c], 0)))
                    t = T(str(S[c]), 11, WHITE).move_to(cells[c])
                    numbers.add(t)
                    anims.append(FadeIn(t))
                else:
                    anims.append(cells[c].animate.set_fill(C_BAD, 0.5))
            self.play(*anims, run_time=0.25)
        self.hold()

        # traceback
        i, j = bestc
        path = [(i, j)]
        while (i, j) != (0, 0):
            for a, b, s in [(i - 1, j - 1, (M_ if i > 0 and j > 0 and q[i - 1] == d[j - 1] else MM_)), (i - 1, j, G_), (i, j - 1, G_)]:
                if (a, b) in kept and S[(a, b)] + s == S[(i, j)]:
                    i, j = a, b
                    path.append((i, j))
                    break
        path = path[::-1]
        line = VMobject().set_stroke(YELLOW, 5)
        line.set_points_as_corners([cells[c].get_center() for c in path])
        self.say(f"The best score is {best}. Tracing back from that cell gives the alignment. The path steps off the diagonal once, and that step is a gap. A one-diagonal extension would have stopped right there.",
                 ShowCreation(line), run_time=2.0)
        total = n * m
        frac = T(f"cells computed: {len(kept) + len(pruned)} of {total}", 26, C_HI).move_to([3.2, 1.2, 0])
        self.play(FadeIn(frac))
        self.say("Only a small region near each promising pair of hits is explored, not the whole matrix, so the dynamic programming stays cheap. BLAST 2.0 is about 3 times faster than the original and finds gapped alignments.")


# =====================================================================
class A3_Significance(NarratedScene):
    """Slides 65 to 67: complexity, bit score, E-value and P-value."""

    def construct(self):
        title = self.corner_title("How surprising is an HSP?")
        self.add(title)

        c = MathT(r"O\big(n\,w\,|A|^{w}\;+\;M\;+\;h\,n^{2}\big)", 42)
        c.shift(2.3 * UP)
        self.say("First the complexity. Building seeds costs n w times the alphabet to the w. Finding hits costs M. Extending costs h times n squared. Unlike Smith-Waterman, M is not multiplied by n.",
                 Write(c), run_time=2)
        self.play(FadeOut(c))

        f1 = MathT(r"S'=\frac{\lambda S-\ln K}{\ln 2}", 40)
        f2 = MathT(r"E=\frac{K\,M\,n}{e^{\lambda S}}", 40)
        f1.shift(2.2 * UP + 3 * LEFT)
        f2.shift(2.2 * UP + 3 * RIGHT)
        self.say("For each HSP BLAST reports two numbers. The bit score S prime normalizes the raw score so searches with different matrices can be compared.",
                 Write(f1), run_time=1.5)
        self.say("The E-value is the expected number of HSPs with score at least S that appear by pure chance. It grows with the database size M and with the query size n, and it falls exponentially with S.",
                 Write(f2), run_time=1.5)
        f3 = MathT(r"E=M\,n\,2^{-S'}", 52, YELLOW)
        f3.shift(0.4 * UP)
        self.say("Combine them: e to the lambda S equals K times 2 to the S prime. Substituting gives one easy formula. Each extra bit halves the E-value.",
                 ReplacementTransform(VGroup(f1.copy(), f2.copy()), f3), run_time=1.5)
        self.play(FadeOut(f1), FadeOut(f2))

        # table
        mn = 3e11
        scene_note = MathT(r"M\,n=3\times 10^{11}\approx 2^{38}", 30, GREY_B).next_to(f3, DOWN, buff=0.4)
        self.say("For a feeling of the scale, take a query of 300 residues and a database of 10^9 residues. M n is about 2 to the 38.", FadeIn(scene_note))
        rows = VGroup()
        for k, s in enumerate((30, 38, 40, 50, 60)):
            e = mn * 2.0 ** (-s)
            rows.add(VGroup(T(f"S' = {s}", 28, C_QUERY), T(f"E = {e:.2g}", 28, C_GOOD if e < 0.01 else C_HI)).arrange(RIGHT, buff=1.2))
        rows.arrange(DOWN, buff=0.22).next_to(scene_note, DOWN, buff=0.35)
        self.say("A bit score of 40 gives E of about 0.27. A bit score of 50 gives about 0.0003. Ten more bits divides E by 1024.",
                 FadeIn(rows, lag_ratio=0.2), run_time=2)
        self.play(FadeOut(VGroup(rows, scene_note, f3)))

        p = MathT(r"P=1-e^{-E}", 50, YELLOW).shift(2.0 * UP)
        self.say("The E-value is not a P-value. E can be larger than 1. P is the probability of at least one chance alignment this good, and the count behaves like a Poisson variable with mean E.",
                 Write(p), run_time=1.5)
        tab = VGroup()
        for E in (0.01, 0.1, 1, 5):
            tab.add(VGroup(T(f"E = {E}", 28, C_HI), T(f"P = {1 - np.exp(-E):.4f}", 28, C_GOOD)).arrange(RIGHT, buff=1.2))
        tab.arrange(DOWN, buff=0.25).shift(0.2 * DOWN)
        self.say("For small E the two are almost equal. For E equal to 5, P is already 0.993, so the hit tells you nothing. An E-value near zero means the HSP is very unlikely to be a coincidence.",
                 FadeIn(tab, lag_ratio=0.2), run_time=2)
