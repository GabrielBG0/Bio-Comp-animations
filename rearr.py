"""Pure-Python helpers for aula 6 (genome rearrangements).

Everything here is plain Python, no manim. The animations call these functions,
so what you see on screen is computed, not typed in.

Conventions
-----------
* An "extended" permutation is a list that starts with 0 and ends with n + 1,
  e.g. [0, 3, 1, 2, 4].  Positions are 0-based indices into that list, so the
  slides' rho(i, j) is `revx(ext, i, j)` for 1 <= i <= j <= n.
* Signed permutations use negative numbers for minus signs.
"""
from collections import deque
from itertools import product


def ext_of(p, signed=False):
    return [0] + list(p) + [len(p) + 1]


def plain(e):
    return e[1:-1]


def fmt(v, signed=False):
    if signed:
        return ("+" if v >= 0 else "-") + str(abs(v))
    return str(v)


# --------------------------------------------------------------- operations
def revx(e, i, j, signed=False):
    seg = e[i:j + 1][::-1]
    if signed:
        seg = [-x for x in seg]
    return e[:i] + seg + e[j + 1:]


def transpx(e, i, j, k):
    """rho(i, j, k): exchange the blocks e[i..j-1] and e[j..k-1]."""
    return e[:i] + e[j:k] + e[i:j] + e[k:]


def flip(e, i):
    """Pancake flip rho(1, i)."""
    return revx(e, 1, i)


# ------------------------------------------------------------- breakpoints
def adj_unsigned(a, b):
    return abs(a - b) == 1


def breakpoints(e):
    return [x for x in range(len(e) - 1) if not adj_unsigned(e[x], e[x + 1])]


def b_count(e):
    return len(breakpoints(e))


def bp_pancake(e):
    """Breakpoints for prefix reversals: the pair (e0, e1) always counts."""
    bps = set(breakpoints(e))
    bps.add(0)
    return len(bps)


def bt_count(e):
    """Breakpoints for transpositions: only +1 steps are adjacencies."""
    return sum(1 for x in range(len(e) - 1) if e[x + 1] - e[x] != 1)


def strips(e):
    """Maximal runs of adjacencies. Returns (start, end, 'inc'|'dec')."""
    out = []
    s = 0
    N = len(e)
    for x in range(1, N + 1):
        if x == N or not adj_unsigned(e[x - 1], e[x]):
            end = x - 1
            if end > s:
                d = "inc" if e[s + 1] > e[s] else "dec"
            else:
                d = "inc" if (s == 0 or s == N - 1) else "dec"
            out.append((s, end, d))
            s = x
    return out


def greedy_bp(e):
    """Kececioglu-Sankoff greedy with breakpoints. Returns steps."""
    steps = []
    cur = list(e)
    guard = 0
    while cur != sorted(cur) and guard < 100:
        guard += 1
        st = strips(cur)
        dec = [x for (s, t, d) in st if d == "dec" for x in cur[s:t + 1]]
        before = b_count(cur)
        info = {}
        if not dec:
            bps = breakpoints(cur)
            i, j = bps[0] + 1, bps[1]
            kind = "none"
            nxt = revx(cur, i, j)
        else:
            k = min(dec)
            pa, pb = cur.index(k), cur.index(k - 1)
            lo, hi = min(pa, pb), max(pa, pb)
            ik, jk = lo + 1, hi
            cand_k = revx(cur, ik, jk)
            has_dec = any(d == "dec" for (_, _, d) in strips(cand_k))
            if has_dec:
                i, j, kind, nxt = ik, jk, "k", cand_k
                info = dict(k=k)
            else:
                l = max(dec)
                pa, pb = cur.index(l), cur.index(l + 1)
                lo, hi = min(pa, pb), max(pa, pb)
                i, j = lo, hi - 1
                kind, nxt = "l", revx(cur, i, j)
                info = dict(k=k, l=l)
        steps.append(dict(i=i, j=j, kind=kind, before=cur, after=nxt,
                          db=b_count(nxt) - before, info=info))
        cur = nxt
    return steps


# ------------------------------------------------------- signed reversals
def oriented_pairs(e):
    """Pairs of positions (p, q), p < q, with consecutive absolute values and
    different signs. e is a signed extended permutation (ends are positive)."""
    pos = {abs(v): x for x, v in enumerate(e)}
    out = []
    for v in range(len(e) - 1):
        p, q = pos[v], pos[v + 1]
        if (e[p] >= 0) != (e[q] >= 0):
            a, b = min(p, q), max(p, q)
            out.append((a, b))
    return out


def pair_reversal(e, p, q):
    s = e[p] + e[q]
    return (p, q - 1) if s == 1 else (p + 1, q)


def score(e):
    return len(oriented_pairs(e))


def bergeron(e, prefer=None):
    """Oriented reversals of maximum score. Returns steps with candidates."""
    steps = []
    cur = list(e)
    t = 0
    while oriented_pairs(cur):
        assert t < 200, 'bergeron did not terminate'
        cands = []
        for (p, q) in oriented_pairs(cur):
            i, j = pair_reversal(cur, p, q)
            res = revx(cur, i, j, True)
            cands.append(dict(pair=(p, q), rev=(i, j), result=res, score=score(res)))
        best = max(c["score"] for c in cands)
        top = [c for c in cands if c["score"] == best]
        pick = top[0]
        if prefer is not None and t < len(prefer):
            for c in top:
                if c["rev"] == prefer[t]:
                    pick = c
                    break
        steps.append(dict(before=cur, cands=cands, pick=pick, best=best))
        cur = pick["result"]
        t += 1
    return steps, cur


# ----------------------------------------------------------- cycle graph
def cg_vertices(e):
    """Vertex sequence of the cycle graph: list of blocks of extremities."""
    n = len(e) - 2
    blocks = []
    for v in e:
        if v == 0:
            blocks.append([0])
        elif v == n + 1:
            blocks.append([2 * n + 1])
        elif v > 0:
            blocks.append([2 * v - 1, 2 * v])
        else:
            blocks.append([2 * -v, 2 * -v - 1])
    return blocks


def cg_edges(e):
    blocks = cg_vertices(e)
    n = len(e) - 2
    black = [(blocks[t][-1], blocks[t + 1][0]) for t in range(len(blocks) - 1)]
    gray = [(2 * m, 2 * m + 1) for m in range(n + 1)]
    return blocks, black, gray


def cg_cycles(e):
    """Cycles as lists of vertex ids; each alternates black, gray."""
    blocks, black, gray = cg_edges(e)
    bmap, gmap = {}, {}
    for a, b in black:
        bmap[a], bmap[b] = b, a
    for a, b in gray:
        gmap[a], gmap[b] = b, a
    seen = set()
    cycles = []
    for v0 in sorted(bmap):
        if v0 in seen:
            continue
        cyc = []
        v = v0
        while True:
            seen.add(v)
            w = bmap[v]
            seen.add(w)
            cyc.append((v, w))  # a black edge, traversed v -> w
            v = gmap[w]
            if v == v0:
                break
        cycles.append(cyc)
    return cycles


def c_count(e):
    return len(cg_cycles(e))


# ---------------------------------------------------------------- hurdles
def framed(e):
    """Shortest framed interval from each start (circular). e is a positive
    extended permutation. Returns list of (start, length) in positions."""
    N = len(e)
    out = []
    for p in range(N):
        i = e[p]
        c = m = 0
        for step in range(1, N):
            x = (p + step) % N
            r = (e[x] - i) % N
            c += 1
            m = max(m, r)
            if c >= 2 and r == m and r == c:
                out.append((p, step + 1))
                break
    return out


def _cover(p, L, N):
    return {(p + t) % N for t in range(L)}


def hurdles(e):
    N = len(e)
    cands = framed(e)
    res = []
    for (p, L) in cands:
        cs = _cover(p, L, N)
        inside = False
        for (q, M) in cands:
            if (q, M) == (p, L):
                continue
            if _cover(q, M, N) < cs:
                inside = True
                break
        if not inside:
            res.append((p, L))
    return res, cands


# -------------------------------------------------------- exact distances
def signed_distance(p):
    """BFS over signed permutations (use only for n <= 6)."""
    n = len(p)
    start = tuple(p)
    goal = tuple(range(1, n + 1))
    if start == goal:
        return 0
    seen = {start}
    dq = deque([(start, 0)])
    while dq:
        cur, d = dq.popleft()
        for i in range(n):
            for j in range(i, n):
                nxt = cur[:i] + tuple(-x for x in reversed(cur[i:j + 1])) + cur[j + 1:]
                if nxt == goal:
                    return d + 1
                if nxt not in seen:
                    seen.add(nxt)
                    dq.append((nxt, d + 1))
    return None


def unsigned_dist(p, ops="r"):
    n = len(p)
    start = tuple(p)
    goal = tuple(range(1, n + 1))
    if start == goal:
        return 0
    seen = {start}
    dq = deque([(start, 0)])
    while dq:
        cur, d = dq.popleft()
        nxts = []
        if "r" in ops:
            for i in range(n):
                for j in range(i + 1, n):
                    nxts.append(cur[:i] + cur[i:j + 1][::-1] + cur[j + 1:])
        if "t" in ops:
            for i in range(n + 1):
                for j in range(i + 1, n + 1):
                    for k in range(j + 1, n + 2):
                        nxts.append(cur[:i] + cur[j:k] + cur[i:j] + cur[k:])
        for nx in nxts:
            if nx == goal:
                return d + 1
            if nx not in seen:
                seen.add(nx)
                dq.append((nx, d + 1))
    return None


# ----------------------------------------------------------- transpositions
def naive_transp(e):
    steps = []
    cur = list(e)
    while cur != sorted(cur):
        i = next(x for x in range(1, len(cur)) if cur[x] != cur[x - 1] + 1)
        j = cur.index(cur[i - 1] + 1)
        k = j + 1
        while k < len(cur) and cur[k] == cur[k - 1] + 1:
            k += 1
        nxt = transpx(cur, i, j, k)
        steps.append(dict(i=i, j=j, k=k, before=cur, after=nxt))
        cur = nxt
    return steps


def rt_unsigned(e):
    """Walter-Dias-Meidanis 3-approximation for reversals + transpositions."""
    steps = []
    cur = list(e)
    while cur != sorted(cur):
        m = 0
        while m + 1 < len(cur) and cur[m + 1] == cur[m] + 1:
            m += 1
        p = cur.index(m + 1)
        st = next((s, t, d) for (s, t, d) in strips(cur) if s <= p <= t)
        s, t, d = st
        if d == "dec" and t > s and p == t:
            nxt = revx(cur, m + 1, p)
            steps.append(dict(kind="rev", args=(m + 1, p), before=cur, after=nxt, m=m))
        else:
            assert p == s, (cur, m, p, st)
            nxt = transpx(cur, m + 1, p, t + 1)
            steps.append(dict(kind="tr", args=(m + 1, p, t + 1), before=cur, after=nxt, m=m))
        cur = nxt
    return steps


# -------------------------------------------------------------------- DCJ
def dcj_adjacency_graph(chrs_a, chrs_b):
    """Each genome is a list of linear chromosomes (lists of signed genes).
    Returns (vertices_a, vertices_b, edges, components)."""
    def exts(g):
        return [(g, "t"), (g, "h")] if g > 0 else [(-g, "h"), (-g, "t")]

    def verts(chrs):
        out = []
        for ch in chrs:
            seq = []
            for g in ch:
                seq += exts(g)
            out.append(frozenset([seq[0]]))
            for t in range(1, len(seq) - 1, 2):
                out.append(frozenset([seq[t], seq[t + 1]]))
            out.append(frozenset([seq[-1]]))
        return out

    va, vb = verts(chrs_a), verts(chrs_b)
    edges = []
    for x, a in enumerate(va):
        for e in a:
            y = next(t for t, b in enumerate(vb) if e in b)
            edges.append((x, y, e))
    # components
    adj = {}
    for x, y, e in edges:
        adj.setdefault(("a", x), []).append(("b", y))
        adj.setdefault(("b", y), []).append(("a", x))
    seen, comps = set(), []
    for v in adj:
        if v in seen:
            continue
        comp, dq = [], deque([v])
        seen.add(v)
        while dq:
            u = dq.popleft()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    dq.append(w)
        ne = sum(1 for x, y, e in edges if ("a", x) in comp)
        degs = [len(adj[u]) for u in comp]
        kind = "cycle" if all(d == 2 for d in degs) else ("odd path" if ne % 2 == 1 else "even path")
        comps.append(dict(vertices=comp, edges=ne, kind=kind))
    return va, vb, edges, comps


def dcj_distance(chrs_a, chrs_b):
    n = sum(len(c) for c in chrs_a)
    _, _, _, comps = dcj_adjacency_graph(chrs_a, chrs_b)
    c = sum(1 for k in comps if k["kind"] == "cycle")
    i = sum(1 for k in comps if k["kind"] == "odd path")
    return n - (c + i / 2), c, i


if __name__ == "__main__":
    # ------------------------------------------------ self tests = slide numbers
    e = [0, 4, 3, 2, 5, 6, 7, 1, 8, 9]
    assert b_count(e) == 4
    seq = [(2, 4), (4, 5), (4, 7), (3, 4), (5, 6)]
    bs = [4, 5, 7, 6, 6, 4]
    cur = e
    for (i, j), b in zip(seq, bs[1:]):
        cur = revx(cur, i, j)
        assert b_count(cur) == b, (cur, b)
    # greedy with breakpoints
    st = greedy_bp([0, 3, 4, 1, 2, 7, 8, 5, 6, 9])
    assert [(s["i"], s["j"]) for s in st] == [(1, 2), (3, 4), (1, 6), (3, 8), (1, 8)], [(s["i"], s["j"]) for s in st]
    assert [s["db"] for s in st] == [0, -1, -1, -1, -2]
    st = greedy_bp([0, 4, 1, 2, 6, 5, 3, 7])
    assert [(s["i"], s["j"]) for s in st] == [(4, 6), (2, 4), (1, 4)], [(s["i"], s["j"]) for s in st]
    # signed
    steps, fin = bergeron([0, 2, -1, 4, -6, -3, 7, -5, 8],
                          prefer=[(1, 2), (4, 7), (5, 7), (2, 2), (3, 5), (3, 6), (6, 6)])
    assert len(steps) == 7 and fin == list(range(9)), (len(steps), fin)
    assert [s["pick"]["rev"] for s in steps] == [(1, 2), (4, 7), (5, 7), (2, 2), (3, 5), (3, 6), (6, 6)]
    # cycle graph
    e = [0, 1, -2, -5, -4, -6, 3, 7]
    assert c_count(e) == 4
    assert c_count([0, 1, 2, -5, -4, -6, 3, 7]) == 5
    # hurdles
    e = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8, 10, 12, 11, 13]
    hs, cands = hurdles(e)
    assert sorted(hs) == [(0, 4), (3, 5), (7, 4), (10, 4)], hs
    hs, _ = hurdles([0, 2, 1, 3, 5, 4, 6])
    assert sorted(hs) == [(0, 4), (3, 4)], hs
    # DCJ
    d, c, i = dcj_distance([[1, 4], [-2, -3, -5]], [[1, 2, 3, 4, 5]])
    assert (d, c, i) == (3, 1, 2)
    # transpositions
    st = naive_transp([0, 5, 2, 1, 4, 3, 6])
    assert len(st) == 4
    st = rt_unsigned([0, 4, 6, 1, 7, 3, 2, 5, 8])
    assert [s["kind"] for s in st] == ["tr", "rev", "tr", "rev"], [s["kind"] for s in st]
    # exact distances
    assert signed_distance((-3, -1, 2)) == 3
    assert signed_distance((4, -2, -5, 3, 1)) == 3
    assert signed_distance((2, 1, 3, 5, 4)) == 6
    assert signed_distance((1, -2, -5, -4, -6, 3)) == 3
    assert signed_distance((1, 2, -5, -4, -6, 3)) == 2
    print("rearr self tests ok")
