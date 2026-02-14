import json
import os
import sys
import tempfile
import warnings
from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd
from pypdf import PdfReader, PdfWriter, Transformation
from upsetplot import UpSet
from upsetplot.data import from_contents

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import COLOR_CELLTYPES, COLOR_STAGES

warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

res_dir = "results/networks/analysis"
fig_dir = "figures/networks/"
os.makedirs(fig_dir, exist_ok=True)
input_json = os.path.join(res_dir, "pfc_celltype_stage_targets.json")
input_csv = os.path.join(res_dir, "pfc_celltype_stage_targets.csv")

_subplot_w = 1
_subplot_h = 1
_src_w_in = 3
_src_h_in = 2

cell_w_pt = _subplot_w * 100
cell_h_pt = _subplot_h * 80

color_celltypes = OrderedDict(COLOR_CELLTYPES)
color_stages = COLOR_STAGES

if os.path.exists(input_json):
    with open(input_json) as f:
        celltype_data = json.load(f)
else:
    csv_data = pd.read_csv(input_csv)
    celltype_data = {}
    for celltype in csv_data["CellType"].unique():
        celltype_data[celltype] = {}
        cell_data = csv_data[csv_data["CellType"] == celltype]
        for stage in cell_data["Stage"].unique():
            stage_data = cell_data[cell_data["Stage"] == stage]
            celltype_data[celltype][stage] = stage_data["Target"].unique().tolist()

celltypes = sorted(celltype_data.keys())
_all_stage_set = set(
    [s for ct in celltype_data for s in celltype_data[ct]]
    + [f"S{i}" for i in range(1, 16)]
)
all_stages = sorted(
    _all_stage_set,
    key=lambda s: int(s[1:]) if (s[:1] == "S" and s[1:].isdigit()) else 0,
    reverse=True,
)

tmp_dir = tempfile.mkdtemp()
pdf_paths = []
for celltype in celltypes:
    stage_targets = celltype_data[celltype]
    complete = {s: stage_targets.get(s, []) for s in all_stages}
    non_empty = [s for s in all_stages if len(complete[s]) > 0]
    if not non_empty:
        continue
    contents = {s: set(complete[s]) for s in non_empty}
    try:
        df = from_contents(contents)
    except Exception:
        continue
    fig = plt.figure(figsize=(_src_w_in, _src_h_in))
    upset = UpSet(
        df,
        sort_by="degree",
        sort_categories_by="input",
        max_subset_rank=15,
        show_counts="{:d}",
        element_size=15,
        totals_plot_elements=1,
    )
    for cat in non_empty:
        upset.style_categories(cat, bar_facecolor=color_stages.get(cat, "#888"))
    subplots = upset.plot(fig)
    if "totals" in subplots:
        ax_totals = subplots["totals"]
        ax_totals.tick_params(axis="x", rotation=30, labelsize=9)
        for lbl in ax_totals.get_xticklabels():
            lbl.set_ha("right")
    fig.subplots_adjust(top=0.88, wspace=0.15)
    title_col = color_celltypes.get(celltype, "#333333")
    fig.text(
        0.07,  # x position
        0.95,  # y position
        celltype,
        color=title_col,
        fontsize=11,
        ha="left",
        va="top",
        transform=fig.transFigure,
    )
    tmp_pdf = os.path.join(tmp_dir, f"plot_{len(pdf_paths)}.pdf")
    fig.savefig(tmp_pdf, format="pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    pdf_paths.append(tmp_pdf)

n_plots = len(pdf_paths)
if n_plots > 0:
    n_row = (n_plots + 1) // 2
    max_w_pt = 0.0
    max_h_pt = 0.0
    for p in pdf_paths:
        r = PdfReader(p)
        pg = r.pages[0]
        max_w_pt = max(max_w_pt, float(pg.mediabox.width))
        max_h_pt = max(max_h_pt, float(pg.mediabox.height))
    scale = min(cell_w_pt / max_w_pt, cell_h_pt / max_h_pt)
    row_h_pt = cell_h_pt * 0.85
    writer = PdfWriter()
    blank = writer.add_blank_page(width=2 * cell_w_pt, height=n_row * row_h_pt)
    scaled_h_pt = scale * max_h_pt
    for idx in range(n_plots):
        i = idx // 2
        j = idx % 2
        reader = PdfReader(pdf_paths[idx])
        src = reader.pages[0]
        tx = j * cell_w_pt
        ty = (n_row - i) * row_h_pt - scaled_h_pt
        blank.merge_transformed_page(
            src, Transformation().scale(scale, scale).translate(tx, ty)
        )
    out_path = os.path.join(fig_dir, "upset_pfc_celltypes.pdf")
    with open(out_path, "wb") as f:
        writer.write(f)
    for p in pdf_paths:
        try:
            os.remove(p)
        except OSError:
            pass
    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass
