from manimlib import *


class Hello(Scene):
    def construct(self):
        circle = Circle().set_fill(BLUE, opacity=0.5)
        self.play(ShowCreation(circle))
        self.play(circle.animate.shift(2 * RIGHT))
        self.play(FadeOut(circle))
