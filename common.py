"""Shared helpers for the Bio-Comp animations (manimgl / 3b1b manim).

NarratedScene adds burned-in subtitles (manimgl has no subtitle feature of its
own) and also writes an .srt file per scene into ./subtitles/ so the same text
can be loaded as a soft subtitle track in a video player.

Usage inside construct():

    self.say("Text of the lecture.", FadeIn(x), run_time=1)   # show text, play, then wait for reading
    self.say("More text.")                                    # only show text and wait
    self.say("Text", wait=False); self.play(...); self.hold() # text while several animations run
"""
import os
import textwrap

from manimlib import *

# ---------------------------------------------------------------- palette
C_QUERY = BLUE_C
C_DB = YELLOW
C_GOOD = GREEN
C_BAD = RED
C_HI = ORANGE
C_DIM = GREY_B
C_SEED = TEAL
BG_PANEL = "#1b1b1b"


def T(text, size=30, color=WHITE, **kw):
    """Plain text shortcut."""
    return Text(text, font_size=size, color=color, **kw)


def MathT(tex, size=36, color=WHITE, **kw):
    return Tex(tex, font_size=size, fill_color=color, **kw)


def fmt_srt(t):
    h = int(t // 3600)
    m = int((t % 3600) // 60)
    s = t % 60
    return f"{h:02d}:{m:02d}:{int(s):02d},{int(round((s - int(s)) * 1000)):03d}"


class NarratedScene(Scene):
    WORDS_PER_SECOND = 3.0      # reading speed used to decide how long a subtitle stays
    MIN_SHOW = 1.8
    SUB_WIDTH_CHARS = 62        # wrap width of the subtitle
    SUB_FONT = 26

    # ----------------------------------------------------------- life cycle
    def setup(self):
        self._sub = None
        self._sub_bg = None
        self._sub_t0 = 0.0
        self._sub_need = 0.0
        self._srt = []          # [start, end, text]
        self._pending = []

    def tear_down(self):
        self.clear_subtitle()
        os.makedirs("subtitles", exist_ok=True)
        path = os.path.join("subtitles", f"{self.__class__.__name__}.srt")
        with open(path, "w", encoding="utf-8") as f:
            for i, (a, b, txt) in enumerate(self._srt, 1):
                f.write(f"{i}\n{fmt_srt(a)} --> {fmt_srt(max(b, a + 0.5))}\n{txt}\n\n")
        super().tear_down()

    # ------------------------------------------------------------ subtitles
    MAX_CHUNK_CHARS = 124       # about two lines of SUB_WIDTH_CHARS

    @classmethod
    def split_chunks(cls, text):
        """Split a long narration into pieces that fit in two subtitle lines."""
        import re
        sentences = re.split(r"(?<=[.!?:;])\s+", text.strip())
        pieces = []
        for s in sentences:
            while len(s) > cls.MAX_CHUNK_CHARS:
                cut = s.rfind(", ", 0, cls.MAX_CHUNK_CHARS)
                if cut < 40:
                    cut = s.rfind(" ", 0, cls.MAX_CHUNK_CHARS)
                else:
                    cut += 1
                pieces.append(s[:cut].strip())
                s = s[cut:].strip()
            if s:
                pieces.append(s)
        chunks = []
        for p in pieces:
            if chunks and len(chunks[-1]) + 1 + len(p) <= cls.MAX_CHUNK_CHARS:
                chunks[-1] += " " + p
            else:
                chunks.append(p)
        return chunks or [text]

    def _close_current(self):
        if self._srt and self._srt[-1][1] is None:
            self._srt[-1][1] = self.time

    def _make_sub(self, text):
        wrapped = "\n".join(textwrap.wrap(text, self.SUB_WIDTH_CHARS))
        sub = Text(wrapped, font_size=self.SUB_FONT, color=WHITE)
        if sub.get_width() > 12.4:
            sub.set_width(12.4)
        sub.to_edge(DOWN, buff=0.22)
        bg = Rectangle(width=sub.get_width() + 0.5, height=sub.get_height() + 0.3)
        bg.set_fill(BLACK, opacity=0.78).set_stroke(width=0)
        bg.move_to(sub)
        grp = VGroup(bg, sub)
        grp.fix_in_frame()
        return grp

    def _show_chunk(self, text):
        self._close_current()
        new = self._make_sub(text)
        if self._sub is not None:
            self.remove(self._sub)
        self.add(new)
        self._sub = new
        self._sub_t0 = self.time
        self._sub_need = max(self.MIN_SHOW, len(text.split()) / self.WORDS_PER_SECOND)
        self._srt.append([self.time, None, text])

    def set_subtitle(self, text):
        chunks = self.split_chunks(text)
        self._pending = chunks[1:]
        self._show_chunk(chunks[0])

    def clear_subtitle(self):
        self._close_current()
        if self._sub is not None:
            self.remove(self._sub)
            self._sub = None
        self._pending = []

    def hold(self, extra=0.0):
        """Wait until the current subtitle (and any remaining pieces of it) had time to be read."""
        left = self._sub_need - (self.time - self._sub_t0)
        self.wait(max(left, 0.0))
        while self._pending:
            self._show_chunk(self._pending.pop(0))
            self.wait(self._sub_need)
        if extra:
            self.wait(extra)

    def say(self, text, *anims, run_time=1.0, wait=True, extra=0.0, **kw):
        self.set_subtitle(text)
        if anims:
            self.play(*anims, run_time=run_time, **kw)
        if wait:
            self.hold(extra)

    # ----------------------------------------------------------- title card
    def title_card(self, title, subtitle=None, time=2.5):
        t = Text(title, font_size=54)
        group = VGroup(t)
        if subtitle:
            s = Text(subtitle, font_size=28, color=GREY_B)
            s.next_to(t, DOWN, buff=0.4)
            group.add(s)
        group.move_to(ORIGIN)
        self.play(FadeIn(group, shift=0.3 * UP))
        self.wait(time)
        self.play(FadeOut(group))

    def corner_title(self, text, color=BLUE_B):
        t = Text(text, font_size=34, color=color)
        t.to_corner(UL, buff=0.35)
        t.fix_in_frame()
        return t


# ------------------------------------------------------------ small widgets
def cell(label, size=0.62, fill=None, color=WHITE, font=28, stroke=WHITE):
    sq = Square(side_length=size)
    sq.set_stroke(stroke, 2)
    if fill is not None:
        sq.set_fill(fill, opacity=0.35)
    txt = Text(str(label), font_size=font, color=color)
    txt.move_to(sq)
    g = VGroup(sq, txt)
    g.sq, g.txt = sq, txt
    return g


def row_of_cells(labels, size=0.62, font=28, buff=0.0, **kw):
    cells = [cell(l, size=size, font=font, **kw) for l in labels]
    g = VGroup(*cells).arrange(RIGHT, buff=buff)
    g.cells = cells
    return g


def surround(mob, color=YELLOW, buff=0.08, width=3):
    r = SurroundingRectangle(mob, buff=buff)
    r.set_stroke(color, width)
    return r
