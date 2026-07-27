#!/usr/bin/env python3
import csv
import html
import math
import statistics
from collections import defaultdict
from pathlib import Path


RESULT_DIR = Path(__file__).resolve().parent
TABLE_DIR = RESULT_DIR / "tables"
FIGURE_DIR = RESULT_DIR / "figures"


NUMERIC_FIELDS = {
    "case_order",
    "nranks",
    "n1",
    "n2",
    "n3",
    "nsys",
    "nrow_min",
    "nrow_max",
    "iter",
    "iterations",
    "total_s_max",
    "total_s_avg",
    "local_compute_s_max",
    "pack_forward_s_max",
    "mpi_forward_s_max",
    "unpack_forward_s_max",
    "reduced_compute_s_max",
    "pack_backward_s_max",
    "mpi_backward_s_max",
    "unpack_backward_s_max",
    "update_compute_s_max",
    "compute_s_max",
    "communication_s_max",
    "packing_s_max",
    "solution_sum",
    "solution_l2",
    "solution_linf",
    "sample_z0",
    "sample_zmid",
    "sample_zlast",
    "expected_value",
    "max_abs_error_to_expected",
}


def latest(pattern):
    matches = sorted(RESULT_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"missing {pattern} under {RESULT_DIR}")
    return matches[-1]


def load_csv(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key in list(row.keys()):
            if key in NUMERIC_FIELDS and row[key] != "":
                try:
                    if key in {"case_order", "nranks", "n1", "n2", "n3", "nsys", "nrow_min", "nrow_max", "iter", "iterations"}:
                        row[key] = int(float(row[key]))
                    else:
                        row[key] = float(row[key])
                except ValueError:
                    pass
        row["tags"] = set(tag.strip() for tag in row.get("study_tags", "").split(";") if tag.strip())
    return rows


def mean(values):
    values = list(values)
    return sum(values) / len(values) if values else float("nan")


def stdev(values):
    values = list(values)
    if len(values) <= 1:
        return 0.0
    return statistics.stdev(values)


def safe_div(a, b):
    if b == 0 or math.isnan(b):
        return float("nan")
    return a / b


def fmt_float(value, digits=3):
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.{digits}f}"


def fmt_sci(value, digits=3):
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.{digits}e}"


def md_table(headers, rows):
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines) + "\n"


def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def write_csv(path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def esc(text):
    return html.escape(str(text), quote=True)


def svg_header(width, height):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>',
        'text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}',
        '.title{font-size:18px;font-weight:700}',
        '.axis{stroke:#52606d;stroke-width:1}',
        '.grid{stroke:#d9e2ec;stroke-width:1}',
        '.tick{font-size:11px;fill:#52606d}',
        '.label{font-size:12px;fill:#334e68}',
        '.point-label{font-size:11px;fill:#243b53;font-weight:600}',
        '.legend{font-size:12px;fill:#243b53}',
        '</style>',
    ]


def svg_footer():
    return ["</svg>"]


def nice_ticks(ymax, count=5):
    if ymax <= 0 or math.isnan(ymax):
        return [0, 1]
    raw = ymax / count
    exponent = math.floor(math.log10(raw))
    fraction = raw / (10 ** exponent)
    if fraction <= 1:
        step = 1 * 10 ** exponent
    elif fraction <= 2:
        step = 2 * 10 ** exponent
    elif fraction <= 5:
        step = 5 * 10 ** exponent
    else:
        step = 10 * 10 ** exponent
    top = math.ceil(ymax / step) * step
    ticks = []
    current = 0
    while current <= top + step * 0.5:
        ticks.append(current)
        current += step
    return ticks


def draw_axes(lines, width, height, left, right, top, bottom, y_ticks, y_max, title, xlabel, ylabel):
    plot_w = width - left - right
    plot_h = height - top - bottom
    lines.append(f'<text x="{width / 2}" y="28" text-anchor="middle" class="title">{esc(title)}</text>')
    lines.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" class="axis"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>')
    for tick in y_ticks:
        y = top + plot_h - safe_div(tick, y_max) * plot_h if y_max else top + plot_h
        lines.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" class="grid"/>')
        lines.append(f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" class="tick">{fmt_float(tick, 2)}</text>')
    lines.append(f'<text x="{left + plot_w / 2}" y="{height - 12}" text-anchor="middle" class="label">{esc(xlabel)}</text>')
    lines.append(f'<text x="16" y="{top + plot_h / 2}" text-anchor="middle" class="label" transform="rotate(-90 16 {top + plot_h / 2})">{esc(ylabel)}</text>')


def line_chart(path, title, xlabel, ylabel, series, x_tick_formatter=None):
    width, height = 900, 520
    left, right, top, bottom = 86, 34, 58, 86
    plot_w = width - left - right
    plot_h = height - top - bottom
    all_points = [p for _, _, points in series for p in points]
    if not all_points:
        return
    xs = [p[0] for p in all_points]
    ys = [p[1] for p in all_points]
    xmin, xmax = min(xs), max(xs)
    if xmin == xmax:
        xmin -= 1
        xmax += 1
    y_ticks = nice_ticks(max(ys) * 1.08)
    y_max = y_ticks[-1]
    lines = svg_header(width, height)
    draw_axes(lines, width, height, left, right, top, bottom, y_ticks, y_max, title, xlabel, ylabel)

    def sx(x):
        return left + safe_div(x - xmin, xmax - xmin) * plot_w

    def sy(y):
        return top + plot_h - safe_div(y, y_max) * plot_h

    x_ticks = sorted(set(xs))
    for x in x_ticks:
        x_pos = sx(x)
        label = x_tick_formatter(x) if x_tick_formatter else str(x)
        lines.append(f'<line x1="{x_pos:.1f}" y1="{top + plot_h}" x2="{x_pos:.1f}" y2="{top + plot_h + 5}" class="axis"/>')
        lines.append(f'<text x="{x_pos:.1f}" y="{top + plot_h + 22}" text-anchor="middle" class="tick">{esc(label)}</text>')

    legend_x = left + 8
    legend_y = top + 16
    for idx, (name, color, points) in enumerate(series):
        pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in sorted(points))
        lines.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for x, y in points:
            lines.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{color}"/>')
        y0 = legend_y + idx * 20
        lines.append(f'<rect x="{legend_x}" y="{y0 - 9}" width="12" height="12" fill="{color}"/>')
        lines.append(f'<text x="{legend_x + 18}" y="{y0 + 2}" class="legend">{esc(name)}</text>')

    lines += svg_footer()
    write_text(path, "\n".join(lines) + "\n")


def log_ticks(vmin, vmax):
    if vmin <= 0 or vmax <= 0:
        return []
    candidates = []
    start = math.floor(math.log10(vmin))
    end = math.ceil(math.log10(vmax))
    for exponent in range(start, end + 1):
        for multiplier in (1, 2, 5):
            value = multiplier * (10 ** exponent)
            if vmin <= value <= vmax:
                candidates.append(value)
    if not candidates:
        return [vmin, vmax]
    if candidates[0] > vmin:
        candidates.insert(0, vmin)
    if candidates[-1] < vmax:
        candidates.append(vmax)
    return candidates


def fmt_log_tick(value):
    if value >= 1000:
        return f"{value:.0f}"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.0f}"
    if value >= 1:
        return f"{value:g}"
    return f"{value:.3g}"


def series_parts(item):
    if len(item) == 3:
        name, color, points = item
        style = {}
    else:
        name, color, points, style = item
    return name, color, points, style


def log_line_chart(path, title, xlabel, ylabel, series, x_tick_formatter=None):
    width, height = 900, 540
    left, right, top, bottom = 92, 34, 62, 92
    plot_w = width - left - right
    plot_h = height - top - bottom
    all_points = [p for item in series for p in series_parts(item)[2]]
    if not all_points:
        return

    xs = [p[0] for p in all_points if p[0] > 0 and p[1] > 0]
    ys = [p[1] for p in all_points if p[0] > 0 and p[1] > 0]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    x_min_plot = xmin / 1.12
    x_max_plot = xmax * 1.12
    y_min_plot = ymin / 1.18
    y_max_plot = ymax * 1.18
    lxmin, lxmax = math.log10(x_min_plot), math.log10(x_max_plot)
    lymin, lymax = math.log10(y_min_plot), math.log10(y_max_plot)

    lines = svg_header(width, height)
    lines.append(f'<text x="{width / 2}" y="28" text-anchor="middle" class="title">{esc(title)}</text>')
    lines.append(f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" class="axis"/>')
    lines.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" class="axis"/>')

    def sx(x):
        return left + safe_div(math.log10(x) - lxmin, lxmax - lxmin) * plot_w

    def sy(y):
        return top + plot_h - safe_div(math.log10(y) - lymin, lymax - lymin) * plot_h

    x_ticks = sorted(set(xs))
    for x in x_ticks:
        x_pos = sx(x)
        label = x_tick_formatter(x) if x_tick_formatter else fmt_log_tick(x)
        lines.append(f'<line x1="{x_pos:.1f}" y1="{top}" x2="{x_pos:.1f}" y2="{top + plot_h}" class="grid"/>')
        lines.append(f'<line x1="{x_pos:.1f}" y1="{top + plot_h}" x2="{x_pos:.1f}" y2="{top + plot_h + 5}" class="axis"/>')
        lines.append(f'<text x="{x_pos:.1f}" y="{top + plot_h + 23}" text-anchor="middle" class="tick">{esc(label)}</text>')

    for tick in log_ticks(y_min_plot, y_max_plot):
        y_pos = sy(tick)
        lines.append(f'<line x1="{left}" y1="{y_pos:.1f}" x2="{left + plot_w}" y2="{y_pos:.1f}" class="grid"/>')
        lines.append(f'<text x="{left - 8}" y="{y_pos + 4:.1f}" text-anchor="end" class="tick">{esc(fmt_log_tick(tick))}</text>')

    lines.append(f'<text x="{left + plot_w / 2}" y="{height - 14}" text-anchor="middle" class="label">{esc(xlabel)} (log scale)</text>')
    lines.append(f'<text x="16" y="{top + plot_h / 2}" text-anchor="middle" class="label" transform="rotate(-90 16 {top + plot_h / 2})">{esc(ylabel)} (log scale)</text>')

    legend_x = left + 8
    legend_y = top + 16
    for idx, item in enumerate(series):
        name, color, points, style = series_parts(item)
        points = sorted([p for p in points if p[0] > 0 and p[1] > 0])
        pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in points)
        dash = ' stroke-dasharray="7 5"' if style.get("dash") else ""
        opacity = f' opacity="{style.get("opacity")}"' if style.get("opacity") else ""
        stroke_width = style.get("stroke_width", 2.5)
        lines.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{stroke_width}"{dash}{opacity}/>')
        if style.get("markers", True):
            for x, y in points:
                lines.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{color}"{opacity}/>')
        point_labels = style.get("point_labels", {})
        for x, y in points:
            if x in point_labels:
                lines.append(f'<text x="{sx(x) + 7:.1f}" y="{sy(y) - 7:.1f}" class="point-label">{esc(point_labels[x])}</text>')
        y0 = legend_y + idx * 20
        lines.append(f'<line x1="{legend_x}" y1="{y0 - 3}" x2="{legend_x + 16}" y2="{y0 - 3}" stroke="{color}" stroke-width="{stroke_width}"{dash}{opacity}/>')
        lines.append(f'<text x="{legend_x + 22}" y="{y0 + 2}" class="legend">{esc(name)}</text>')

    lines += svg_footer()
    write_text(path, "\n".join(lines) + "\n")


def grouped_bar_chart(path, title, xlabel, ylabel, labels, series):
    width, height = 920, 520
    left, right, top, bottom = 88, 34, 58, 112
    plot_w = width - left - right
    plot_h = height - top - bottom
    all_values = [v for _, _, values in series for v in values]
    y_ticks = nice_ticks(max(all_values) * 1.12 if all_values else 1)
    y_max = y_ticks[-1]
    lines = svg_header(width, height)
    draw_axes(lines, width, height, left, right, top, bottom, y_ticks, y_max, title, xlabel, ylabel)

    group_w = plot_w / max(len(labels), 1)
    gap = group_w * 0.2
    bar_area = group_w - gap
    bar_w = bar_area / max(len(series), 1)
    colors = [color for _, color, _ in series]
    for i, label in enumerate(labels):
        x0 = left + i * group_w + gap / 2
        for j, (_, _, values) in enumerate(series):
            value = values[i]
            h = safe_div(value, y_max) * plot_h
            x = x0 + j * bar_w
            y = top + plot_h - h
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(bar_w - 3, 1):.1f}" height="{h:.1f}" fill="{colors[j]}"/>')
        lx = left + i * group_w + group_w / 2
        lines.append(f'<text x="{lx:.1f}" y="{top + plot_h + 24}" text-anchor="middle" class="tick">{esc(label)}</text>')

    legend_x = left + 8
    legend_y = top + 16
    for idx, (name, color, _) in enumerate(series):
        y0 = legend_y + idx * 20
        lines.append(f'<rect x="{legend_x}" y="{y0 - 9}" width="12" height="12" fill="{color}"/>')
        lines.append(f'<text x="{legend_x + 18}" y="{y0 + 2}" class="legend">{esc(name)}</text>')

    lines += svg_footer()
    write_text(path, "\n".join(lines) + "\n")


def stacked_bar_chart(path, title, xlabel, ylabel, labels, stacks):
    width, height = 920, 520
    left, right, top, bottom = 88, 34, 58, 112
    plot_w = width - left - right
    plot_h = height - top - bottom
    totals = [sum(values[i] for _, _, values in stacks) for i in range(len(labels))]
    y_ticks = nice_ticks(max(totals) * 1.12 if totals else 1)
    y_max = y_ticks[-1]
    lines = svg_header(width, height)
    draw_axes(lines, width, height, left, right, top, bottom, y_ticks, y_max, title, xlabel, ylabel)

    group_w = plot_w / max(len(labels), 1)
    bar_w = group_w * 0.52
    for i, label in enumerate(labels):
        x = left + i * group_w + (group_w - bar_w) / 2
        cursor = top + plot_h
        for _, color, values in stacks:
            value = values[i]
            h = safe_div(value, y_max) * plot_h
            cursor -= h
            lines.append(f'<rect x="{x:.1f}" y="{cursor:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
        lx = left + i * group_w + group_w / 2
        lines.append(f'<text x="{lx:.1f}" y="{top + plot_h + 24}" text-anchor="middle" class="tick">{esc(label)}</text>')

    legend_x = left + 8
    legend_y = top + 16
    for idx, (name, color, _) in enumerate(stacks):
        y0 = legend_y + idx * 20
        lines.append(f'<rect x="{legend_x}" y="{y0 - 9}" width="12" height="12" fill="{color}"/>')
        lines.append(f'<text x="{legend_x + 18}" y="{y0 + 2}" class="legend">{esc(name)}</text>')

    lines += svg_footer()
    write_text(path, "\n".join(lines) + "\n")


def group_timing(rows):
    groups = defaultdict(list)
    for row in rows:
        key = (row["run_case_id"], row["implementation"], row["mpi_mode"])
        groups[key].append(row)

    summaries = []
    for key, group_rows in groups.items():
        stable = [r for r in group_rows if r["iter"] >= 1]
        if not stable:
            stable = list(group_rows)
        first = next((r for r in group_rows if r["iter"] == 0), group_rows[0])
        base = stable[0]
        total_values = [r["total_s_max"] for r in stable]
        compute_values = [r["compute_s_max"] for r in stable]
        comm_values = [r["communication_s_max"] for r in stable]
        packing_values = [r["packing_s_max"] for r in stable]
        total_mean = mean(total_values)
        row = {
            "run_case_id": key[0],
            "implementation": key[1],
            "mpi_mode": key[2],
            "study_tags": base.get("study_tags", ""),
            "study_notes": base.get("study_notes", ""),
            "case_order": base["case_order"],
            "nranks": base["nranks"],
            "n1": base["n1"],
            "n2": base["n2"],
            "n3": base["n3"],
            "nsys": base["nsys"],
            "nrow_min": base["nrow_min"],
            "nrow_max": base["nrow_max"],
            "stable_iterations": len(stable),
            "total_s_mean": total_mean,
            "total_s_std": stdev(total_values),
            "total_s_min": min(total_values),
            "total_s_max_observed": max(total_values),
            "total_s_iter0": first["total_s_max"],
            "warmup_ratio": safe_div(first["total_s_max"], total_mean),
            "compute_s_mean": mean(compute_values),
            "communication_s_mean": mean(comm_values),
            "packing_s_mean": mean(packing_values),
            "throughput_mcells_s": safe_div(base["n1"] * base["n2"] * base["n3"], total_mean) / 1.0e6,
            "cv_percent": safe_div(stdev(total_values), total_mean) * 100.0,
        }
        summaries.append(row)
    summaries.sort(key=lambda r: (r["case_order"], r["implementation"], r["mpi_mode"]))
    return summaries


def has_tag(row, tag):
    return tag in set(x.strip() for x in row.get("study_tags", "").split(";") if x.strip())


def select_summary(summaries, **criteria):
    selected = summaries
    for key, value in criteria.items():
        selected = [r for r in selected if r.get(key) == value]
    return selected


def by_key(rows, keys):
    return {tuple(row[k] for k in keys): row for row in rows}


def make_outputs():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    profile_path = latest("tdma_total_profile_merged_*.csv")
    correctness_path = latest("tdma_correctness_merged_*.csv")
    manifest_path = latest("tdma_case_manifest_merged_*.csv")

    profile_rows = load_csv(profile_path)
    correctness_rows = load_csv(correctness_path)
    manifest_rows = load_csv(manifest_path)
    timing = group_timing(profile_rows)

    write_csv(
        TABLE_DIR / "summary_by_case.csv",
        [
            "run_case_id",
            "implementation",
            "mpi_mode",
            "nranks",
            "n1",
            "n2",
            "n3",
            "nsys",
            "nrow_max",
            "stable_iterations",
            "total_s_mean",
            "total_s_std",
            "total_s_iter0",
            "warmup_ratio",
            "compute_s_mean",
            "communication_s_mean",
            "packing_s_mean",
            "throughput_mcells_s",
            "cv_percent",
            "study_tags",
        ],
        [
            [
                r["run_case_id"],
                r["implementation"],
                r["mpi_mode"],
                r["nranks"],
                r["n1"],
                r["n2"],
                r["n3"],
                r["nsys"],
                r["nrow_max"],
                r["stable_iterations"],
                f"{r['total_s_mean']:.12g}",
                f"{r['total_s_std']:.12g}",
                f"{r['total_s_iter0']:.12g}",
                f"{r['warmup_ratio']:.12g}",
                f"{r['compute_s_mean']:.12g}",
                f"{r['communication_s_mean']:.12g}",
                f"{r['packing_s_mean']:.12g}",
                f"{r['throughput_mcells_s']:.12g}",
                f"{r['cv_percent']:.12g}",
                r["study_tags"],
            ]
            for r in timing
        ],
    )

    implementations = sorted({r["implementation"] for r in profile_rows})
    profile_modes = sorted({r["mpi_mode"] for r in profile_rows})
    correctness_impl = sorted({r["implementation"] for r in correctness_rows})
    coverage_rows = [
        ["profile rows", len(profile_rows)],
        ["profile cases", len(timing)],
        ["correctness rows", len(correctness_rows)],
        ["manifest rows", len(manifest_rows)],
        ["implementations in profile", ", ".join(implementations)],
        ["implementations in correctness", ", ".join(correctness_impl)],
        ["MPI modes", ", ".join(profile_modes)],
        ["Fortran timing rows present", "yes" if "fortran-original" in implementations else "no"],
    ]
    write_text(TABLE_DIR / "0_data_coverage.md", md_table(["item", "value"], coverage_rows))

    correctness_sorted = sorted(correctness_rows, key=lambda r: (r["case_order"], r["implementation"], r["mpi_mode"]))
    max_error = max((r["max_abs_error_to_expected"] for r in correctness_sorted), default=float("nan"))
    correctness_rows_md = [
        [
            r["run_case_id"],
            r["implementation"],
            r["mpi_mode"],
            r["nranks"],
            f"{r['n1']}x{r['n2']}x{r['n3']}",
            fmt_sci(r["max_abs_error_to_expected"], 3),
        ]
        for r in correctness_sorted
    ]
    correctness_intro = f"Max absolute error to expected value across available rows: {fmt_sci(max_error, 3)}\n\n"
    write_text(
        TABLE_DIR / "1_correctness_summary.md",
        correctness_intro
        + md_table(["case", "implementation", "mode", "np", "grid", "max_abs_error"], correctness_rows_md),
    )

    central = [
        r
        for r in timing
        if r["n1"] == 128 and r["n2"] == 128 and r["n3"] == 4096 and r["nranks"] in {2, 4, 8}
    ]
    central_rows = [
        [
            r["run_case_id"],
            r["implementation"],
            r["mpi_mode"],
            r["nranks"],
            fmt_float(r["total_s_mean"] * 1000, 3),
            fmt_float(r["compute_s_mean"] * 1000, 3),
            fmt_float(r["communication_s_mean"] * 1000, 3),
            fmt_float(r["packing_s_mean"] * 1000, 3),
            fmt_float(r["throughput_mcells_s"], 1),
            fmt_float(r["warmup_ratio"], 1),
        ]
        for r in central
    ]
    write_text(
        TABLE_DIR / "2_central_case_timing.md",
        md_table(
            [
                "case",
                "implementation",
                "mode",
                "np",
                "total_ms",
                "compute_ms",
                "comm_ms",
                "packing_ms",
                "throughput_Mcells_s",
                "iter0_to_stable",
            ],
            central_rows,
        ),
    )

    strong = [r for r in timing if has_tag(r, "strong_scaling") and r["mpi_mode"] == "device"]
    strong_rows_out = []
    strong_csv_rows = []
    for size in sorted({(r["n1"], r["n2"], r["n3"]) for r in strong}):
        subset = sorted([r for r in strong if (r["n1"], r["n2"], r["n3"]) == size], key=lambda r: r["nranks"])
        base = next((r for r in subset if r["nranks"] == 2), None)
        for r in subset:
            speedup = safe_div(base["total_s_mean"], r["total_s_mean"]) if base else float("nan")
            efficiency = safe_div(speedup, r["nranks"] / 2.0)
            row = [
                f"{size[0]}x{size[1]}x{size[2]}",
                r["nranks"],
                fmt_float(r["total_s_mean"] * 1000, 3),
                fmt_float(speedup, 3),
                fmt_float(efficiency * 100, 1),
                fmt_float(r["throughput_mcells_s"], 1),
            ]
            strong_rows_out.append(row)
            strong_csv_rows.append(row)
    write_text(
        TABLE_DIR / "3_strong_scaling.md",
        md_table(["grid", "np", "total_ms", "speedup_2base", "efficiency_2base_percent", "throughput_Mcells_s"], strong_rows_out),
    )
    write_csv(TABLE_DIR / "strong_scaling.csv", ["grid", "np", "total_ms", "speedup_2base", "efficiency_2base_percent", "throughput_Mcells_s"], strong_csv_rows)

    weak = [r for r in timing if (has_tag(r, "weak_nrow_scaling") or has_tag(r, "weak_nsys_scaling")) and r["mpi_mode"] == "device"]
    weak_rows = []
    for path, tag in [("weak_nrow", "weak_nrow_scaling"), ("weak_nsys", "weak_nsys_scaling")]:
        for r in sorted([row for row in weak if has_tag(row, tag)], key=lambda x: (x["nranks"], x["n1"], x["n2"], x["n3"])):
            weak_rows.append(
                [
                    path,
                    r["nranks"],
                    f"{r['n1']}x{r['n2']}x{r['n3']}",
                    r["nsys"],
                    r["nrow_max"],
                    fmt_float(r["total_s_mean"] * 1000, 3),
                    fmt_float(r["throughput_mcells_s"], 1),
                ]
            )
    write_text(TABLE_DIR / "4_weak_scaling.md", md_table(["path", "np", "grid", "nsys", "nrow", "total_ms", "throughput_Mcells_s"], weak_rows))

    nsys_rows = [
        r
        for r in timing
        if has_tag(r, "nsys_sensitivity") and r["mpi_mode"] == "device" and r["n3"] == 4096
    ]
    nsys_rows_md = [
        [
            r["nranks"],
            r["nsys"],
            f"{r['n1']}x{r['n2']}x{r['n3']}",
            r["nrow_max"],
            fmt_float(r["total_s_mean"] * 1000, 3),
            fmt_float(r["throughput_mcells_s"], 1),
        ]
        for r in sorted(nsys_rows, key=lambda x: (x["nranks"], x["nsys"]))
    ]
    write_text(TABLE_DIR / "5_nsys_sensitivity.md", md_table(["np", "nsys", "grid", "nrow", "total_ms", "throughput_Mcells_s"], nsys_rows_md))

    nrow_rows = [
        r
        for r in timing
        if has_tag(r, "nrow_sensitivity") and r["mpi_mode"] == "device" and r["n1"] == 128 and r["n2"] == 128
    ]
    nrow_rows_md = [
        [
            r["nranks"],
            r["nrow_max"],
            f"{r['n1']}x{r['n2']}x{r['n3']}",
            fmt_float(r["total_s_mean"] * 1000, 3),
            fmt_float(safe_div(r["communication_s_mean"], r["total_s_mean"]) * 100, 1),
            fmt_float(r["throughput_mcells_s"], 1),
        ]
        for r in sorted(nrow_rows, key=lambda x: (x["nranks"], x["nrow_max"]))
    ]
    write_text(TABLE_DIR / "6_nrow_sensitivity.md", md_table(["np", "nrow", "grid", "total_ms", "comm_percent_of_total", "throughput_Mcells_s"], nrow_rows_md))

    mode_rows = [
        r
        for r in timing
        if has_tag(r, "mpi_mode_comparison") or has_tag(r, "mpi_mode_device_reference") or has_tag(r, "host_fallback")
    ]
    mode_rows = [r for r in mode_rows if r["n1"] == 128 and r["n2"] == 128 and r["n3"] == 4096 and r["nranks"] in {2, 4, 8}]
    mode_by_np = defaultdict(dict)
    for r in mode_rows:
        mode_by_np[r["nranks"]][r["mpi_mode"]] = r
    mpi_rows_md = []
    for np in sorted(mode_by_np):
        device = mode_by_np[np].get("device")
        host = mode_by_np[np].get("host")
        ratio = safe_div(host["total_s_mean"], device["total_s_mean"]) if host and device else float("nan")
        for mode in ["device", "host"]:
            r = mode_by_np[np].get(mode)
            if not r:
                continue
            mpi_rows_md.append(
                [
                    np,
                    mode,
                    fmt_float(r["total_s_mean"] * 1000, 3),
                    fmt_float(r["compute_s_mean"] * 1000, 3),
                    fmt_float(r["communication_s_mean"] * 1000, 3),
                    fmt_float(r["packing_s_mean"] * 1000, 3),
                    fmt_float(ratio, 3) if mode == "host" else "-",
                ]
            )
    write_text(TABLE_DIR / "7_mpi_mode_comparison.md", md_table(["np", "mode", "total_ms", "compute_ms", "comm_ms", "packing_ms", "host_over_device"], mpi_rows_md))

    repro = sorted(timing, key=lambda r: r["cv_percent"], reverse=True)[:10]
    repro_rows = [
        [
            r["run_case_id"],
            r["mpi_mode"],
            r["nranks"],
            fmt_float(r["total_s_mean"] * 1000, 3),
            fmt_float(r["total_s_std"] * 1000, 4),
            fmt_float(r["cv_percent"], 2),
            fmt_float(r["warmup_ratio"], 1),
        ]
        for r in repro
    ]
    write_text(TABLE_DIR / "8_reproducibility_top_cv.md", md_table(["case", "mode", "np", "total_ms_mean", "std_ms", "cv_percent", "iter0_to_stable"], repro_rows))

    warmup = sorted(timing, key=lambda r: r["warmup_ratio"], reverse=True)[:10]
    warmup_rows = [
        [
            r["run_case_id"],
            r["mpi_mode"],
            r["nranks"],
            fmt_float(r["total_s_iter0"] * 1000, 3),
            fmt_float(r["total_s_mean"] * 1000, 3),
            fmt_float(r["warmup_ratio"], 1),
        ]
        for r in warmup
    ]
    write_text(TABLE_DIR / "9_warmup_effect_top.md", md_table(["case", "mode", "np", "iter0_ms", "stable_mean_ms", "iter0_to_stable"], warmup_rows))

    # Figures
    colors = {
        "blue": "#2f80ed",
        "green": "#27ae60",
        "orange": "#f2994a",
        "red": "#eb5757",
        "purple": "#9b51e0",
        "gray": "#7b8794",
    }

    strong_series_time = []
    strong_series_throughput = []
    for color, size in zip([colors["blue"], colors["green"]], sorted({(r["n1"], r["n2"], r["n3"]) for r in strong})):
        subset = sorted([r for r in strong if (r["n1"], r["n2"], r["n3"]) == size], key=lambda r: r["nranks"])
        base = next((r for r in subset if r["nranks"] == 2), None)
        if not base:
            continue
        measured_time = [(r["nranks"], r["total_s_mean"] * 1000) for r in subset]
        ideal_time = [(r["nranks"], base["total_s_mean"] * 1000 * safe_div(2.0, r["nranks"])) for r in subset]
        measured_throughput = [(r["nranks"], r["throughput_mcells_s"]) for r in subset]
        ideal_throughput = [(r["nranks"], base["throughput_mcells_s"] * safe_div(r["nranks"], 2.0)) for r in subset]
        efficiency_labels = {}
        for r in subset:
            speedup = safe_div(base["total_s_mean"], r["total_s_mean"])
            efficiency = safe_div(speedup, r["nranks"] / 2.0)
            efficiency_labels[r["nranks"]] = f"{efficiency * 100:.1f}%"

        label = f"{size[0]}x{size[1]}x{size[2]}"
        strong_series_time.append((label, color, measured_time, {"point_labels": efficiency_labels}))
        strong_series_time.append((f"{label} ideal", color, ideal_time, {"dash": True, "markers": False, "opacity": "0.45", "stroke_width": 2}))
        strong_series_throughput.append((label, color, measured_throughput, {"point_labels": efficiency_labels}))
        strong_series_throughput.append((f"{label} ideal", color, ideal_throughput, {"dash": True, "markers": False, "opacity": "0.45", "stroke_width": 2}))
    log_line_chart(FIGURE_DIR / "1_strong_scaling_time.svg", "Strong Scaling: Stable Total Time vs Ideal", "MPI ranks / GPUs", "total_s_max mean (ms)", strong_series_time)
    log_line_chart(FIGURE_DIR / "2_strong_scaling_throughput.svg", "Strong Scaling: Throughput vs Ideal", "MPI ranks / GPUs", "million cells / s", strong_series_throughput)

    central_device = sorted(
        [
            r
            for r in central
            if r["mpi_mode"] == "device" and r["implementation"] == "cuda-cxx"
        ],
        key=lambda r: r["nranks"],
    )
    labels = [f"np={r['nranks']}" for r in central_device]
    stacked_bar_chart(
        FIGURE_DIR / "3_phase_breakdown_central_device.svg",
        "CUDA C++ Device Mode: Phase Breakdown",
        "central case 128x128x4096",
        "phase time mean (ms)",
        labels,
        [
            ("compute", colors["blue"], [r["compute_s_mean"] * 1000 for r in central_device]),
            ("communication", colors["orange"], [r["communication_s_mean"] * 1000 for r in central_device]),
            ("packing", colors["green"], [r["packing_s_mean"] * 1000 for r in central_device]),
        ],
    )

    mode_labels = [f"np={np}" for np in sorted(mode_by_np)]
    grouped_bar_chart(
        FIGURE_DIR / "4_mpi_mode_device_vs_host.svg",
        "MPI Mode Comparison: Device vs Host",
        "central case 128x128x4096",
        "total_s_max mean (ms)",
        mode_labels,
        [
            ("device", colors["blue"], [mode_by_np[np]["device"]["total_s_mean"] * 1000 for np in sorted(mode_by_np)]),
            ("host", colors["red"], [mode_by_np[np]["host"]["total_s_mean"] * 1000 for np in sorted(mode_by_np)]),
        ],
    )

    nsys_series = []
    for color, np in zip([colors["blue"], colors["green"], colors["purple"]], [2, 4, 8]):
        subset = sorted([r for r in nsys_rows if r["nranks"] == np], key=lambda r: r["nsys"])
        if subset:
            nsys_series.append((f"np={np}", color, [(r["nsys"], r["total_s_mean"] * 1000) for r in subset]))
    log_line_chart(FIGURE_DIR / "5_nsys_sensitivity.svg", "Nsys Sensitivity at n3=4096", "nsys = n1*n2", "total_s_max mean (ms)", nsys_series, lambda x: f"{int(x)}")

    nrow_series = []
    for color, np in zip([colors["blue"], colors["green"], colors["purple"]], [2, 4, 8]):
        subset = sorted([r for r in nrow_rows if r["nranks"] == np], key=lambda r: r["nrow_max"])
        if subset:
            nrow_series.append((f"np={np}", color, [(r["nrow_max"], r["total_s_mean"] * 1000) for r in subset]))
    log_line_chart(FIGURE_DIR / "6_nrow_sensitivity.svg", "Nrow Sensitivity at n1=n2=128", "local nrow = n3/np", "total_s_max mean (ms)", nrow_series, lambda x: f"{int(x)}")

    weak_nrow = sorted([r for r in weak if has_tag(r, "weak_nrow_scaling")], key=lambda r: r["nranks"])
    weak_nsys = sorted([r for r in weak if has_tag(r, "weak_nsys_scaling")], key=lambda r: r["nranks"])
    log_line_chart(
        FIGURE_DIR / "7_weak_scaling_time.svg",
        "Weak Scaling Paths",
        "MPI ranks / GPUs",
        "total_s_max mean (ms)",
        [
            ("weak nrow path", colors["blue"], [(r["nranks"], r["total_s_mean"] * 1000) for r in weak_nrow]),
            ("weak nsys path", colors["orange"], [(r["nranks"], r["total_s_mean"] * 1000) for r in weak_nsys]),
        ],
    )

    warmup_central = central_device
    grouped_bar_chart(
        FIGURE_DIR / "8_warmup_effect_central_device.svg",
        "Warm-up Effect: First Iteration vs Stable Mean",
        "central case 128x128x4096",
        "total_s_max (ms)",
        [f"np={r['nranks']}" for r in warmup_central],
        [
            ("iter 0", colors["red"], [r["total_s_iter0"] * 1000 for r in warmup_central]),
            ("iter 1-9 mean", colors["blue"], [r["total_s_mean"] * 1000 for r in warmup_central]),
        ],
    )


if __name__ == "__main__":
    make_outputs()
