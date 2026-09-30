import csv
import json
from collections import defaultdict
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import cooler

def extend_tss_clamped(input_file, extension_len, output_file=None):

    # Read features grouped by chromosome, preserving all columns
    by_chrom = defaultdict(list)
    with open(input_file, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 3 or row[0].startswith('#'):
                continue
            by_chrom[row[0]].append(row)

    out_rows = []
    for rows in by_chrom.values():
        rows.sort(key=lambda r: (int(r[1]), int(r[2])))
        starts = [int(r[1]) for r in rows]
        ends = [int(r[2]) for r in rows]
        n = len(rows)
        for i, row in enumerate(rows):
            # Left edge: floor at 0 and at the midpoint to the previous feature
            left_limit = 0 if i == 0 else (ends[i - 1] + starts[i]) // 2
            new_start = max(starts[i] - extension_len, left_limit)

            # Right edge: cap at the midpoint to the next feature
            right_limit = (ends[i] + extension_len if i == n - 1
                           else (ends[i] + starts[i + 1]) // 2)
            new_end = min(ends[i] + extension_len, right_limit)

            new_row = list(row)
            new_row[1] = str(new_start)
            new_row[2] = str(new_end)
            out_rows.append(new_row)

    if output_file is None:
        output_file = Path(input_file).with_name(
            f"{Path(input_file).stem}_extended_clamped_{extension_len}.bed"
        )
    output_file = Path(output_file)

    with open(output_file, 'w', newline='') as out:
        csv.writer(out, delimiter='\t').writerows(out_rows)

    return str(output_file)
            
def create_range_dict(input_bed, min_regions=3):

    # Collect boundary intervals per chromosome
    boundaries = defaultdict(list)
    with open(input_bed, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 3 or row[0].startswith('#'):
                continue
            boundaries[row[0]].append((int(row[1]), int(row[2])))

    range_dict = {}
    for chrom, feats in boundaries.items():
        feats.sort()  # order by start (then end)
        # Spans strictly between consecutive boundaries
        regions = [
            [prev_end, next_start]
            for (_, prev_end), (next_start, _) in zip(feats, feats[1:])
            if prev_end < next_start  # skip overlapping/adjacent features
        ]
        if len(regions) >= min_regions:
            range_dict[chrom] = regions

    return range_dict
            
def generate_df_from_cooler(mcool_path, resolution, input_range_dict,
                            sample_name, max_group=5):

    clr = cooler.Cooler(f"{mcool_path}::/resolutions/{resolution}")
    sel = clr.matrix(balance=True, as_pixels=True, join=True)
    chromnames = set(clr.chromnames)

    parts = []
    for chrom, regions in input_range_dict.items():
        if chrom not in chromnames or not regions:
            continue
        px = sel.fetch(chrom)                 # cis pixels: both anchors in `chrom`
        if px.empty:
            continue

        bal = px["balanced"].to_numpy()
        x_pos = (px["start1"].to_numpy() + px["end1"].to_numpy()) // 2
        y_pos = (px["start2"].to_numpy() + px["end2"].to_numpy()) // 2

        # Vectorised region lookup: searchsorted(starts, pos, "right") - 1, then a
        # bounds check. starts[i] <= pos is guaranteed for i >= 0, so it only re-checks
        # pos <= ends[i]. -1 marks "outside every region / in a gap".
        reg = np.asarray(regions, dtype=np.int64)
        starts, ends = reg[:, 0], reg[:, 1]

        def region_idx(pos):
            i = np.searchsorted(starts, pos, side="right") - 1
            ok = i >= 0
            ic = np.where(ok, i, 0)
            inside = ok & (pos <= ends[ic])
            return np.where(inside, ic, -1)

        idx1 = region_idx(x_pos)
        idx2 = region_idx(y_pos)
        last_idx = len(regions) - 1

        keep = (idx1 >= 0) & (idx2 >= 0) & ~np.isnan(bal)
        keep &= ~(((idx1 == 0) & (idx2 == 0)) |
                  ((idx1 == last_idx) & (idx2 == last_idx)))
        if not keep.any():
            continue

        idx1, idx2 = idx1[keep], idx2[keep]
        bal = bal[keep]
        dist = np.abs(x_pos[keep] - y_pos[keep])

        raw_diff = idx2 - idx1
        if max_group is None:
            group = raw_diff.astype(str)
        else:
            group = np.where(raw_diff > max_group, f">{max_group}",
                             raw_diff.astype(str))
        dist_rounded = (np.round(dist / resolution) * resolution).astype(np.int64)

        parts.append(pd.DataFrame({"group": group,
                                   "distance_rounded": dist_rounded,
                                   "balanced": bal}))

    cols = ["distance_rounded", "group", "n", "mean", "se", "sample", "resolution"]
    if not parts:
        df = pd.DataFrame(columns=cols)
        df["sample"] = df["sample"].astype(str)
        return df

    allpx = pd.concat(parts, ignore_index=True)
    # sem() == sqrt(var(ddof=1)/n), identical to the Welford sqrt((M2/(n-1))/n),
    # and both yield NaN at n == 1.
    df = (allpx.groupby(["group", "distance_rounded"])["balanced"]
                .agg(n="count", mean="mean", se="sem")
                .reset_index())
    df["sample"] = sample_name
    df["resolution"] = int(resolution)
    df["group"] = df["group"].astype(str)
    return df[cols]

def save_regions(range_dict, output_path):
    # Serialise a create_range_dict() result to JSON so pipeline steps can share it as a file.
    with open(output_path, "w") as f:
        json.dump(range_dict, f)
    return str(output_path)

def load_regions(input_path):
    # Load a regions JSON written by save_regions back into the {chrom: [[start, end], ...]} structure generate_df_from_cooler expects.
    with open(input_path) as f:
        return json.load(f)

def combine_replicates(rep_dfs, condition):

    cat = pd.concat(rep_dfs, ignore_index=True)
    out = (
        cat.groupby(["resolution", "group", "distance_rounded"])["mean"]
        .agg(mean="mean", se="sem", n_rep="count")
        .reset_index()
    )
    out["condition"] = condition
    return out[["condition", "resolution", "group", "distance_rounded",
                "n_rep", "mean", "se"]]

def fold_change_replicates(rep_dfs, condition, ref_group="0"):

    cat = pd.concat(rep_dfs, ignore_index=True)
    ref = (
        cat[cat["group"] == ref_group]
        [["resolution", "sample", "distance_rounded", "mean"]]
        .rename(columns={"mean": "ref_mean"})
    )
    merged = cat.merge(ref, on=["resolution", "sample", "distance_rounded"])
    merged = merged[merged["group"] != ref_group].copy()
    merged["log2fc"] = np.log2(merged["mean"] / merged["ref_mean"])

    out = (
        merged.groupby(["resolution", "group", "distance_rounded"])["log2fc"]
        .agg(log2fc="mean", se="sem", n_rep="count")
        .reset_index()
    )
    out["condition"] = condition
    return out[["condition", "resolution", "group", "distance_rounded",
                "n_rep", "log2fc", "se"]]

def _group_sort_key(g):
    # Order separation groups numerically ("2" before "10") and always place the
    # ">N" overflow bucket last, regardless of how many groups are plotted.
    g = str(g)
    if g.startswith(">"):
        return (1, int(g[1:]))
    return (0, int(g))


def _fmt_bp(x, _):
    if x >= 1e6:
        return f"{x / 1e6:g} Mb"
    if x >= 1e3:
        return f"{x / 1e3:g} kb"
    return f"{int(x)}"


def _style_log_x(ax, dense_ticks=False):
    ax.set_xscale("log")
    if dense_ticks:
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt_bp))
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs="all"))
        ax.tick_params(axis="x", labelrotation=45)
    else:
        ax.xaxis.set_major_locator(ticker.LogLocator(base=10))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt_bp))
        ax.tick_params(axis="x", labelrotation=0)
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(axis="both", labelsize=7)


def _condition_colors(conditions):
    # Deterministic condition -> colour map, shared by every figure so a given
    # condition always renders in the same colour across all plot types.
    conditions = sorted(conditions)
    palette = sns.color_palette("tab10", n_colors=len(conditions))
    return {c: palette[i] for i, c in enumerate(conditions)}


def _group_colors(groups):
    groups = sorted(groups, key=_group_sort_key)
    palette = sns.color_palette("tab10", n_colors=len(groups))
    return {g: palette[i] for i, g in enumerate(groups)}


def _fc_ylim(fc_df, Z=1.0, pad=0.05):
    # y-limits for a log2FC panel spanning the CI band (log2fc +/- Z*se) with a small linear margin.
    # Returns None when the frame has no finite fold change.
    if fc_df.empty:
        return None
    lo = (fc_df["log2fc"] - Z * fc_df["se"]).min()
    hi = (fc_df["log2fc"] + Z * fc_df["se"]).max()
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return None
    lo = min(lo, 0.0)
    hi = max(hi, 0.0)
    span = hi - lo
    margin = span * pad if span > 0 else 1.0
    return (lo - margin, hi + margin)


def _mean_ylim(combined_df, Z=1.0, pad=0.05):
    # y-limits for a balanced-mean panel. 
    # The panels are log-scaled, so the band (mean +/- Z*se) is clipped to strictly positive values and padded
    # multiplicatively (a constant factor in log space).
    if combined_df.empty:
        return None
    lower = combined_df["mean"] - Z * combined_df["se"]
    upper = combined_df["mean"] + Z * combined_df["se"]
    lower = lower[lower > 0]
    upper = upper[np.isfinite(upper)]
    if lower.empty or upper.empty:
        return None
    lo, hi = float(lower.min()), float(upper.max())
    if not (lo > 0 and hi > lo):
        return None
    factor = (hi / lo) ** pad
    return (lo / factor, hi * factor)


def compute_axis_limits(combined_df, fc_df, Z=1.0):
    # Global y-limits per figure type, computed once over *all* TSS sets so the per-TSS figures render on a common axis.
    # Each figure type keeps its own range (the fold-change figures are not forced to match one another).
    return {
        "overlay_ps":         _mean_ylim(combined_df, Z),
        "per_condition_main": _mean_ylim(combined_df, Z),
        "per_condition_fc":   _fc_ylim(fc_df, Z),
        "overlay_fc":         _fc_ylim(fc_df, Z),
    }


def _write_placeholder_figure(output_name, message):
    # Some TSS sets yield no inter-TSS regions (e.g. a chromosome dropped by min_regions), so the combined/fc tables arrive empty.
    # Emit a stub figure instead of feeding plt.subplots a zero dimension, which would raise and fail the whole workflow.
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11, wrap=True)
    # fig.savefig(f"{output_name}.png", dpi=300, bbox_inches="tight")  # PNG disabled: SVG only
    fig.savefig(f"{output_name}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_per_condition(combined_df, fc_df, output_name, Z=1.0,
                       main_ylim=None, fc_ylim=None):
    if combined_df.empty:
        _write_placeholder_figure(output_name, "No data: this TSS set produced no inter-TSS regions.")
        return
    conditions = sorted(combined_df["condition"].unique())
    resolutions = sorted(combined_df["resolution"].unique())
    groups = sorted(combined_df["group"].unique(), key=_group_sort_key)

    group_colors = _group_colors(groups)

    n_rows, n_cols = len(conditions), len(resolutions)
    fig, axes = plt.subplots(
        n_rows * 2, n_cols,
        figsize=(4.5 * n_cols, 5.5 * n_rows),
        sharex="col", squeeze=False,
    )

    # Share y within each *panel type*: 
    # every mean panel shares one y-axis and every log2FC panel shares another, so curves are comparable across conditions.
    main_axes = [axes[r * 2][c]     for r in range(n_rows) for c in range(n_cols)]
    fc_axes   = [axes[r * 2 + 1][c] for r in range(n_rows) for c in range(n_cols)]
    for group in (main_axes, fc_axes):
        for ax in group[1:]:
            ax.sharey(group[0])
    # sharey() alone keeps the inner tick labels;
    for r in range(n_rows * 2):
        for c in range(1, n_cols):
            axes[r][c].tick_params(labelleft=False)

    for r_idx, cond in enumerate(conditions):
        for c_idx, res in enumerate(resolutions):
            ax_main = axes[r_idx * 2][c_idx]
            ax_fc = axes[r_idx * 2 + 1][c_idx]

            m = combined_df[(combined_df["condition"] == cond) &
                            (combined_df["resolution"] == res)]
            f = fc_df[(fc_df["condition"] == cond) &
                      (fc_df["resolution"] == res)]

            # gray band over the genomic-distance range where groups 0 and 1
            # overlap (measured from the mean data, applied to both panels)
            x0 = m.loc[m["group"] == "0", "distance_rounded"]
            x1 = m.loc[m["group"] == "1", "distance_rounded"]
            if not x0.empty and not x1.empty:
                lo = max(x0.min(), x1.min())
                hi = min(x0.max(), x1.max())
                if lo < hi:
                    for ax in (ax_main, ax_fc):
                        ax.axvspan(lo, hi, color="0.85", alpha=0.5, linewidth=0, zorder=0)

            for grp in groups:
                g = m[m["group"] == grp].sort_values("distance_rounded")
                if g.empty:
                    continue
                x, y, ci = g["distance_rounded"].values, g["mean"].values, Z * g["se"].values
                ls = "--" if str(grp) == "0" else "-"
                ax_main.plot(x, y, color=group_colors[grp], linewidth=1.2, linestyle=ls, label=grp)
                ax_main.fill_between(x, np.maximum(y - ci, 1e-10), y + ci,
                                     alpha=0.2, color=group_colors[grp])

            _style_log_x(ax_main)
            ax_main.set_yscale("log")
            if c_idx == 0:
                ax_main.set_ylabel("balanced mean (log)", fontsize=8,
                                    color="black")
            if r_idx == 0:
                ax_main.set_title(f"resolution={res}", fontsize=9, fontweight="bold")

            for grp in groups:
                g = f[f["group"] == grp].sort_values("distance_rounded")
                if g.empty:
                    continue
                x, y, ci = g["distance_rounded"].values, g["log2fc"].values, Z * g["se"].values
                ls = "--" if str(grp) == "0" else "-"
                ax_fc.plot(x, y, color=group_colors[grp], linewidth=1.2, linestyle=ls, label=grp)
                ax_fc.fill_between(x, y - ci, y + ci, alpha=0.2, color=group_colors[grp])

            ax_fc.axhline(0, color="black", linewidth=0.8, linestyle="--")
            _style_log_x(ax_fc)
            if c_idx == 0:
                ax_fc.set_ylabel("log2 FC", fontsize=8)
            if r_idx == n_rows - 1:
                ax_fc.set_xlabel("genomic distance (bp, log)", fontsize=8, color="black")

    # Apply the shared cross-TSS limits 
    if main_ylim is not None:
        main_axes[0].set_ylim(*main_ylim)
    if fc_ylim is not None:
        fc_axes[0].set_ylim(*fc_ylim)

    fig.suptitle("P(s) and FC per condition",
                 fontsize=10)
    plt.tight_layout(rect=[0, 0, 0.92, 0.96])

    for r_idx, cond in enumerate(conditions):
        main_pos = axes[r_idx * 2][0].get_position()      # top-left (mean) panel
        fc_pos   = axes[r_idx * 2 + 1][0].get_position()  # bottom-left (log2FC) panel
        y_center = (fc_pos.y0 + main_pos.y1) / 2
        fig.text(main_pos.x0 - 0.055, y_center, cond, ha="center", va="center",
                 rotation=90, fontsize=11, fontweight="bold", color="black")

        # read handles from this condition's mean panel so group 0 (absent from
        # the log2FC panel) is included
        handles, labels = axes[r_idx * 2][0].get_legend_handles_labels()
        right_edge = axes[r_idx * 2][n_cols - 1].get_position().x1
        fig.legend(handles, labels, title="boundaries\ncrossed", loc="center left",
                   bbox_to_anchor=(right_edge + 0.01, y_center), fontsize=8,
                   title_fontsize=8, frameon=True)

    # plt.savefig(f"{output_name}.png", dpi=300, bbox_inches="tight")  # PNG disabled: SVG only
    plt.savefig(f"{output_name}.svg", bbox_inches="tight")
    plt.close(fig)


def _group_linestyles(groups):
    cycle = ["-", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 1)), (0, (1, 1)),
             (0, (3, 1, 1, 1, 1, 1))]
    styles, i = {}, 0
    for g in sorted(groups, key=_group_sort_key):
        if str(g) == "0":
            styles[g] = "--"
        else:
            styles[g] = cycle[i % len(cycle)]
            i += 1
    return styles


def plot_overlay_ps(combined_df, output_name, Z=1.0, show_se=True, ylim=None):
    # Balanced-mean P(s) decay with every condition x group overlaid on one axes per resolution
    if combined_df.empty:
        _write_placeholder_figure(output_name, "No data: this TSS set produced no inter-TSS regions.")
        return
    conditions = sorted(combined_df["condition"].unique())
    resolutions = sorted(combined_df["resolution"].unique())
    # Drop the ">N" overflow bucket: overlay only the explicit numeric groups.
    combined_df = combined_df[~combined_df["group"].astype(str).str.startswith(">")]
    if combined_df.empty:
        _write_placeholder_figure(output_name, "No data: only the >N overflow group is present.")
        return
    groups = sorted(combined_df["group"].unique(), key=_group_sort_key)

    cond_colors = _condition_colors(conditions)
    group_styles = _group_linestyles(groups)

    n_cols = len(resolutions)
    panel_h = 5.5 / 2      # one per_condition panel row
    overhead_h = 1.0       # title + x tick labels + x-axis label
    fig, axes = plt.subplots(
        1, n_cols,
        figsize=(4.5 * n_cols, panel_h + overhead_h),
        sharex=True, sharey=True, squeeze=False,
    )

    for c_idx, res in enumerate(resolutions):
        ax = axes[0][c_idx]

        res_df = combined_df[combined_df["resolution"] == res]
        x0 = res_df.loc[res_df["group"] == "0", "distance_rounded"]
        x1 = res_df.loc[res_df["group"] == "1", "distance_rounded"]
        if not x0.empty and not x1.empty:
            lo = max(x0.min(), x1.min())
            hi = min(x0.max(), x1.max())
            if lo < hi:
                ax.axvspan(lo, hi, color="0.85", alpha=0.5, linewidth=0, zorder=0)

        for cond in conditions:
            for grp in groups:
                g = combined_df[(combined_df["condition"] == cond) &
                                (combined_df["resolution"] == res) &
                                (combined_df["group"] == grp)].sort_values("distance_rounded")
                if g.empty:
                    continue
                x, y, ci = g["distance_rounded"].values, g["mean"].values, Z * g["se"].values
                ax.plot(x, y, color=cond_colors[cond], linestyle=group_styles[grp],
                        linewidth=1.2)
                if show_se:
                    ax.fill_between(x, np.maximum(y - ci, 1e-10), y + ci,
                                    alpha=0.1, color=cond_colors[cond], linewidth=0)

        _style_log_x(ax, dense_ticks=True)
        ax.set_yscale("log")
        if c_idx == 0:
            ax.set_ylabel("balanced mean (log)", fontsize=8)
        ax.set_title(f"resolution={res}", fontsize=9, fontweight="bold")
        ax.set_xlabel("genomic distance (bp, log)", fontsize=8)

    if ylim is not None:
        axes[0][0].set_ylim(*ylim)

    cond_handles = [plt.Line2D([], [], color=cond_colors[c], linewidth=1.5) for c in conditions]
    grp_handles = [plt.Line2D([], [], color="black", linestyle=group_styles[g], linewidth=1.5)
                   for g in groups]

    leg_kw = dict(fontsize=8, title_fontsize=8, frameon=True, handlelength=1.2)
    leg1 = fig.legend(cond_handles, conditions, title="condition", loc="lower left",
                      bbox_to_anchor=(1.0, 0.51), **leg_kw)
    leg2 = fig.legend(grp_handles, [str(g) for g in groups], title="boundaries\ncrossed",
                      loc="upper left", bbox_to_anchor=(1.0, 0.49), **leg_kw)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    width = 1.2 * max(leg1.get_window_extent(r).width,
                      leg2.get_window_extent(r).width) / fig.bbox.width
    leg1.remove()
    leg2.remove()
    leg1 = fig.legend(cond_handles, conditions, title="condition", loc="lower left",
                      bbox_to_anchor=(1.0, 0.51, width, 0.0), mode="expand", **leg_kw)
    fig.add_artist(leg1)
    leg2 = fig.legend(grp_handles, [str(g) for g in groups], title="boundaries\ncrossed",
                      loc="upper left", bbox_to_anchor=(1.0, 0.49, width, 0.0), mode="expand", **leg_kw)
    leg2.set_alignment("center")  # center the title within the expanded box
    fig.suptitle("P(s) stratified by TSS-boundary separation",
                 fontsize=10, y=0.99)
    plt.tight_layout(rect=[0, 0, 1.0, 0.93])
    # plt.savefig(f"{output_name}.png", dpi=300, bbox_inches="tight")  # PNG disabled: SVG only
    plt.savefig(f"{output_name}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_overlay_fc(fc_df, output_name, Z=1.0, ylim=None):
    if fc_df.empty:
        _write_placeholder_figure(output_name, "No data: this TSS set produced no inter-TSS regions.")
        return
    groups = sorted(fc_df["group"].unique(), key=_group_sort_key)
    resolutions = sorted(fc_df["resolution"].unique())
    conditions = sorted(fc_df["condition"].unique())

    cond_colors = _condition_colors(conditions)

    n_rows, n_cols = len(groups), len(resolutions)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4.5 * n_cols, 3.0 * n_rows),
        sharex=True, sharey=True, squeeze=False,
    )

    for r_idx, grp in enumerate(groups):
        for c_idx, res in enumerate(resolutions):
            ax = axes[r_idx][c_idx]
            for cond in conditions:
                g = fc_df[(fc_df["condition"] == cond) &
                          (fc_df["resolution"] == res) &
                          (fc_df["group"] == grp)].sort_values("distance_rounded")
                if g.empty:
                    continue
                x, y, ci = g["distance_rounded"].values, g["log2fc"].values, Z * g["se"].values
                ax.plot(x, y, color=cond_colors[cond], linewidth=1.2, label=cond)
                ax.fill_between(x, y - ci, y + ci, alpha=0.15, color=cond_colors[cond])

            ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
            _style_log_x(ax, dense_ticks=True)
            if c_idx == 0:
                ax.set_ylabel(f"log2 FC\nboundaries crossed = {grp}", fontsize=8)
            if r_idx == 0:
                ax.set_title(f"resolution={res}", fontsize=9, fontweight="bold")
            if r_idx == n_rows - 1:
                ax.set_xlabel("genomic distance (bp, log)", fontsize=8)

    if ylim is not None:
        axes[0][0].set_ylim(*ylim)

    cond_handles, cond_labels = axes[0][0].get_legend_handles_labels()
    leg_kw = dict(fontsize=8, title_fontsize=8, frameon=True, handlelength=1.2)
    fig.legend(cond_handles, cond_labels, title="condition", loc="center left",
               bbox_to_anchor=(1.0, 0.5), **leg_kw)
    fig.suptitle("log2 fold change relative to group crossing no boundaries",
                 fontsize=10)
    plt.tight_layout(rect=[0, 0, 1.0, 0.97])
    # plt.savefig(f"{output_name}.png", dpi=300, bbox_inches="tight")  # PNG disabled: SVG only
    plt.savefig(f"{output_name}.svg", bbox_inches="tight")
    plt.close(fig)
