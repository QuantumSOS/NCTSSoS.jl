"""
NCTSSoS.jl — Hierarchy Convergence Animation
==============================================
Animates the moment-SOS hierarchy tightening around the true optimum
of a noncommutative polynomial optimization problem.

Run:  manim -pqh hierarchy_convergence.py HierarchyConvergence
GIF:  manim -pql --format=gif hierarchy_convergence.py HierarchyConvergence

Requires: pip install manim
"""

from manim import *
import numpy as np

# Julia palette
JULIA_PURPLE = "#9558B2"
JULIA_RED = "#CB3C33"
JULIA_GREEN = "#389826"
JULIA_BLUE = "#4063D8"

# Softer versions for fills
JULIA_PURPLE_SOFT = "#B07CC8"
JULIA_RED_SOFT = "#E06060"
JULIA_GREEN_SOFT = "#5BBF4A"
JULIA_BLUE_SOFT = "#6B8BE0"

DARK_BG = "#1a1a2e"


def smooth_polygon(points, tension=0.3):
    """Create a smooth closed curve through polygon vertices."""
    from manim import CubicBezier
    n = len(points)
    curves = []
    for i in range(n):
        p0 = points[i]
        p3 = points[(i + 1) % n]
        # Control points for smooth curve
        prev_dir = points[(i + 1) % n] - points[(i - 1) % n]
        next_dir = points[(i + 2) % n] - points[i]
        p1 = p0 + tension * prev_dir
        p2 = p3 - tension * next_dir
        curves.append(CubicBezier(p0, p1, p2, p3))
    group = VGroup(*curves)
    return group


class HierarchyConvergence(Scene):
    """
    Main animation: shows moment-SOS hierarchy bounds converging
    to the true polynomial minimum.
    """

    def construct(self):
        self.camera.background_color = DARK_BG

        # ── Title ──
        title = Text("Moment-SOS Hierarchy", font_size=36, color=WHITE)
        subtitle = Text(
            "Converging to the true optimum",
            font_size=22,
            color=GRAY_B,
        )
        subtitle.next_to(title, DOWN, buff=0.3)
        title_group = VGroup(title, subtitle).to_edge(UP, buff=0.5)

        self.play(Write(title), FadeIn(subtitle, shift=UP * 0.2), run_time=1.5)
        self.wait(0.5)

        # ── Axes ──
        axes = Axes(
            x_range=[-3, 3, 1],
            y_range=[-1, 6, 1],
            x_length=8,
            y_length=4.5,
            axis_config={"color": GRAY_C, "stroke_width": 1.5},
            tips=False,
        ).shift(DOWN * 0.5)

        x_label = axes.get_x_axis_label("x", direction=RIGHT).set_color(GRAY_B)
        y_label = axes.get_y_axis_label("f(x)", direction=UP).set_color(GRAY_B)

        self.play(Create(axes), FadeIn(x_label), FadeIn(y_label), run_time=1)

        # ── True polynomial: f(x) = x⁴ - 2x² + x + 1 ──
        # Has a global minimum around x ≈ -1.13, f ≈ -0.63
        def f(x):
            return x**4 - 2 * x**2 + 0.5 * x + 1.5

        # Find approximate minimum for marking
        xs = np.linspace(-2.5, 2.5, 1000)
        ys = [f(x) for x in xs]
        min_idx = np.argmin(ys)
        x_min, y_min = xs[min_idx], ys[min_idx]

        poly_curve = axes.plot(f, x_range=[-2.2, 2.0], color=WHITE, stroke_width=3)
        poly_label = MathTex(
            r"f(x) = x^4 - 2x^2 + \tfrac{1}{2}x + \tfrac{3}{2}",
            font_size=28,
            color=WHITE,
        ).next_to(axes, DOWN, buff=0.3)

        self.play(Create(poly_curve), FadeIn(poly_label), run_time=2)
        self.wait(0.5)

        # ── True minimum marker ──
        min_dot = Dot(
            axes.c2p(x_min, y_min), radius=0.06, color=JULIA_GREEN, z_index=10
        )
        min_label = MathTex(r"f^*", font_size=24, color=JULIA_GREEN)
        min_label.next_to(min_dot, DOWN + LEFT, buff=0.15)

        # Dashed line at true minimum
        true_min_line = DashedLine(
            axes.c2p(-3, y_min),
            axes.c2p(3, y_min),
            color=JULIA_GREEN,
            stroke_width=1.5,
            stroke_opacity=0.5,
            dash_length=0.1,
        )

        self.play(
            FadeIn(min_dot, scale=2),
            FadeIn(min_label),
            Create(true_min_line),
            run_time=1,
        )
        self.wait(0.5)

        # ── SOS lower bounds (hierarchy levels) ──
        # Each level gives a tighter lower bound
        # Level 1: very loose, Level 2: tighter, etc.
        bounds = [
            (-2.5, JULIA_RED, r"d=1", 0.7),
            (-0.8, JULIA_PURPLE, r"d=2", 0.7),
            (0.15, JULIA_BLUE, r"d=3", 0.7),
            (y_min - 0.05, JULIA_GREEN, r"d=4", 0.7),  # Nearly exact
        ]

        bound_lines = []
        bound_labels = []
        # Track sidebar info
        sidebar_items = []

        # Sidebar header
        sidebar_header = Text("Relaxation bounds", font_size=18, color=GRAY_B)
        sidebar_header.to_corner(UR, buff=0.5).shift(DOWN * 0.8)
        self.play(FadeIn(sidebar_header), run_time=0.5)

        for i, (bound_val, color, label_tex, opacity) in enumerate(bounds):
            # Horizontal line at bound level
            bound_line = DashedLine(
                axes.c2p(-3, bound_val),
                axes.c2p(3, bound_val),
                color=color,
                stroke_width=2,
                stroke_opacity=opacity,
                dash_length=0.15,
            )

            # Shade the region between this bound and the polynomial
            # (shows the "gap" being closed)
            shade_points = []
            x_pts = np.linspace(-2.2, 2.0, 80)
            for x in x_pts:
                shade_points.append(axes.c2p(x, f(x)))
            shade_points.append(axes.c2p(2.0, bound_val))
            shade_points.append(axes.c2p(-2.2, bound_val))

            shade = Polygon(
                *shade_points,
                fill_color=color,
                fill_opacity=0.08,
                stroke_width=0,
            )

            # Label on the right
            level_label = MathTex(label_tex, font_size=22, color=color)
            level_val = MathTex(
                rf"\geq {bound_val:.2f}", font_size=20, color=color
            )
            level_label_pos = sidebar_header.get_bottom() + DOWN * (0.4 + i * 0.5)
            level_label.move_to(level_label_pos)
            level_val.next_to(level_label, RIGHT, buff=0.2)

            if i < len(bounds) - 1:
                self.play(
                    Create(bound_line),
                    FadeIn(shade),
                    FadeIn(level_label),
                    FadeIn(level_val),
                    run_time=1.2,
                )
                self.wait(0.4)
            else:
                # Final level — dramatic convergence
                flash = Flash(
                    min_dot,
                    color=JULIA_GREEN,
                    line_length=0.3,
                    num_lines=12,
                    flash_radius=0.5,
                )
                converge_text = Text(
                    "Converged!", font_size=24, color=JULIA_GREEN
                )
                converge_text.next_to(min_dot, RIGHT, buff=0.4)

                self.play(
                    Create(bound_line),
                    FadeIn(shade),
                    FadeIn(level_label),
                    FadeIn(level_val),
                    run_time=1.2,
                )
                self.play(flash, FadeIn(converge_text, scale=1.5), run_time=1)

            bound_lines.append(bound_line)
            bound_labels.append(VGroup(level_label, level_val))

        self.wait(1)

        # ── Final: show the gap closing to zero ──
        gap_arrow = DoubleArrow(
            axes.c2p(2.3, bounds[-1][0]),
            axes.c2p(2.3, y_min),
            color=YELLOW,
            stroke_width=2,
            buff=0,
            max_tip_length_to_length_ratio=0.5,
        )
        gap_label = MathTex(
            r"\text{gap} \to 0", font_size=22, color=YELLOW
        ).next_to(gap_arrow, RIGHT, buff=0.15)

        self.play(GrowArrow(gap_arrow), FadeIn(gap_label), run_time=1)
        self.wait(2)

        # ── Fade out ──
        self.play(*[FadeOut(mob) for mob in self.mobjects], run_time=1.5)

        # ── End card ──
        end_logo = Text("NCTSSoS.jl", font_size=48, color=WHITE)
        end_tagline = Text(
            "Sparse Noncommutative Polynomial Optimization",
            font_size=20,
            color=GRAY_B,
        )
        end_tagline.next_to(end_logo, DOWN, buff=0.3)
        end_group = VGroup(end_logo, end_tagline)

        self.play(FadeIn(end_group, scale=0.9), run_time=1.5)
        self.wait(2)


class SparsityAnimation(Scene):
    """
    Bonus animation: shows a dense moment matrix being decomposed
    into sparse blocks via clique decomposition.
    """

    def construct(self):
        self.camera.background_color = DARK_BG

        title = Text("Exploiting Sparsity", font_size=36, color=WHITE)
        title.to_edge(UP, buff=0.5)
        self.play(Write(title), run_time=1)

        # ── Dense moment matrix (8x8) ──
        n = 8
        dense_cells = []
        cell_size = 0.45
        matrix_group = VGroup()

        for i in range(n):
            row = []
            for j in range(n):
                cell = Square(
                    side_length=cell_size,
                    fill_color=JULIA_PURPLE,
                    fill_opacity=0.6,
                    stroke_color=WHITE,
                    stroke_width=0.5,
                )
                cell.move_to(
                    np.array(
                        [
                            (j - n / 2 + 0.5) * cell_size,
                            (n / 2 - 0.5 - i) * cell_size,
                            0,
                        ]
                    )
                )
                row.append(cell)
                matrix_group.add(cell)
            dense_cells.append(row)

        matrix_group.shift(LEFT * 2.5 + DOWN * 0.3)
        dense_label = Text("Dense moment matrix", font_size=18, color=GRAY_B)
        dense_label.next_to(matrix_group, DOWN, buff=0.3)

        # PSD symbol
        psd_sym = MathTex(r"M \succeq 0", font_size=28, color=WHITE)
        psd_sym.next_to(matrix_group, UP, buff=0.3)

        self.play(
            *[FadeIn(cell, scale=0.5) for cell in matrix_group],
            FadeIn(dense_label),
            FadeIn(psd_sym),
            run_time=2,
        )
        self.wait(0.5)

        # ── Arrow ──
        arrow = Arrow(LEFT * 0.3, RIGHT * 0.3, color=WHITE, stroke_width=3)
        arrow.shift(DOWN * 0.3)
        arrow_label = Text("clique\ndecomposition", font_size=14, color=GRAY_B)
        arrow_label.next_to(arrow, UP, buff=0.2)

        self.play(GrowArrow(arrow), FadeIn(arrow_label), run_time=0.8)

        # ── Sparse block-diagonal matrix ──
        # Three blocks: 3x3, 3x3, 2x2 (overlapping by 1)
        blocks = [
            (3, JULIA_RED, 0),
            (3, JULIA_GREEN, 2),  # overlaps block 1 by 1 row/col
            (3, JULIA_BLUE, 5),
        ]

        sparse_group = VGroup()
        block_outlines = []

        # First draw the full 8x8 grid with empty cells
        for i in range(n):
            for j in range(n):
                cell = Square(
                    side_length=cell_size,
                    fill_color=DARK_BG,
                    fill_opacity=1,
                    stroke_color=GRAY_D,
                    stroke_width=0.3,
                )
                cell.move_to(
                    np.array(
                        [
                            (j - n / 2 + 0.5) * cell_size,
                            (n / 2 - 0.5 - i) * cell_size,
                            0,
                        ]
                    )
                )
                sparse_group.add(cell)

        # Then fill in the block regions
        for block_size, color, offset in blocks:
            block_cells = VGroup()
            for i in range(block_size):
                for j in range(block_size):
                    cell = Square(
                        side_length=cell_size,
                        fill_color=color,
                        fill_opacity=0.5,
                        stroke_color=WHITE,
                        stroke_width=0.5,
                    )
                    cell.move_to(
                        np.array(
                            [
                                (offset + j - n / 2 + 0.5) * cell_size,
                                (n / 2 - 0.5 - offset - i) * cell_size,
                                0,
                            ]
                        )
                    )
                    block_cells.add(cell)
                    sparse_group.add(cell)

            # Block outline
            outline = SurroundingRectangle(
                block_cells, color=color, stroke_width=2, buff=0
            )
            block_outlines.append(outline)
            sparse_group.add(outline)

        sparse_group.shift(RIGHT * 2.5 + DOWN * 0.3)
        for outline in block_outlines:
            outline.shift(RIGHT * 2.5 + DOWN * 0.3)

        sparse_label = Text("Sparse blocks (cliques)", font_size=18, color=GRAY_B)
        sparse_label.next_to(sparse_group, DOWN, buff=0.3)

        # Size comparison
        size_text = MathTex(
            r"8 \times 8 \;\to\; 3{\times}3 + 3{\times}3 + 3{\times}3",
            font_size=22,
            color=JULIA_GREEN,
        )
        size_text.next_to(sparse_group, UP, buff=0.3)

        self.play(
            *[FadeIn(cell, scale=0.5) for cell in sparse_group],
            FadeIn(sparse_label),
            FadeIn(size_text),
            run_time=2,
        )
        self.wait(2)

        # Highlight cost reduction
        savings = Text("→ Smaller SDPs, faster solve!", font_size=22, color=YELLOW)
        savings.to_edge(DOWN, buff=0.5)
        self.play(FadeIn(savings, shift=UP * 0.3), run_time=1)
        self.wait(2)


class FullAnimation(Scene):
    """
    Combined animation: hierarchy convergence + sparsity decomposition.
    For a complete package overview video.

    Run: manim -pqh hierarchy_convergence.py FullAnimation
    """

    def construct(self):
        self.camera.background_color = DARK_BG

        # Delegate to sub-scenes by replaying their logic
        # (Manim doesn't support nested Scene calls cleanly,
        #  so for a full video, render separately and concatenate,
        #  or inline the logic here.)

        intro = Text("NCTSSoS.jl", font_size=56, color=WHITE)
        tagline = Text(
            "Sparse Noncommutative Polynomial Optimization",
            font_size=24,
            color=GRAY_B,
        )
        tagline.next_to(intro, DOWN, buff=0.4)

        self.play(FadeIn(intro, scale=0.8), run_time=1.5)
        self.play(FadeIn(tagline, shift=UP * 0.2), run_time=1)
        self.wait(2)
        self.play(FadeOut(intro), FadeOut(tagline), run_time=1)
