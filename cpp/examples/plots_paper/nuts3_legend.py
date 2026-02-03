from plotting_settings import plotting_dir, set_fontsize, colors, dpi
from matplotlib.lines import Line2D
import matplotlib.pylab as plt
from matplotlib.legend_handler import HandlerTuple
import os

save_dir = os.path.join(plotting_dir, "plots")

set_fontsize()
legend_elements = [
    # Lines
    Line2D([0], [0], color=colors["Green"], lw=4, label="No NPIs"),
    Line2D([0], [0], color=colors["Teal"], lw=4, label='NPIs continued'),
    Line2D([0], [0], color=colors["Red"], lw=4, label='Strict NPIs'),
    Line2D([0], [0], color=colors["Orange"], lw=4, label='Dynamic NPIs'),

    # Thresholds
    # Line2D([0], [0], color=colors["Red"], lw=1.2, linestyle="--", label="High incidence threshold"),
]

fig, ax = plt.subplots(figsize=(1, 1.2))
ax.axis("off")

ax.legend(
    handles=legend_elements,
    loc="center",
    ncol=1,
    frameon=False
)

plt.savefig(
    os.path.join(save_dir, "nuts3_legend_strategies.png"),
    dpi=dpi,
    bbox_inches="tight"
)


set_fontsize()
legend_elements = [
    # Thresholds
    Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.8, linestyle="-", label="High incidence threshold"),
    Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.5, linestyle="-", label="Low incidence threshold"),
]

fig, ax = plt.subplots(figsize=(3, 0.8))
ax.axis("off")

ax.legend(
    handles=legend_elements,
    loc="center",
    ncol=1,
    frameon=True
)

plt.savefig(
    os.path.join(save_dir, "nuts3_legend_thresholds.png"),
    dpi=dpi,
    bbox_inches="tight",
    pad_inches=0.02
)

set_fontsize()
legend_elements = [
    # Thresholds
    (
        Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.5, linestyle="-"),
        Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.8, linestyle="-"),
    ),
]

labels = [
    "Incidence thresholds",
]
fig, ax = plt.subplots(figsize=(1, 0.4))
ax.axis("off")

ax.legend(
    handles=legend_elements,
    labels=labels,
    loc="center",
    ncol=1,
    frameon=True,
    handler_map={tuple: HandlerTuple(ndivide=None)}
)

plt.savefig(
    os.path.join(save_dir, "nuts3_legend_thresholds_single.png"),
    dpi=dpi,
    bbox_inches="tight",
    pad_inches=0.02
)


set_fontsize()
legend_elements = [
    # Thresholds
    (
        Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.5, linestyle="-"),
        Line2D([0], [0], color=colors["Blue"], lw=4, alpha = 0.8, linestyle="-"),
    ),
    Line2D([0], [0], color=colors["Orange"], lw=0, linestyle="-", marker='o', markersize=12, markeredgewidth=0, markeredgecolor=colors["Orange"], markerfacecolor=colors["Orange"]),
    Line2D([0], [0], color=colors["Orange"], lw=0, linestyle="-", marker='X', markersize=12, markeredgewidth=0, markeredgecolor=colors["Orange"], markerfacecolor=colors["Orange"]),
    Line2D([0], [0], color=colors["Orange"], lw=0, linestyle="-", marker='D', markersize=10, markeredgewidth=0, markeredgecolor=colors["Orange"], markerfacecolor=colors["Orange"]),
    # Line2D([0], [0], color=colors["Teal"], lw=4, linestyle="-"),
    # Line2D([0], [0], color=colors["Red"], lw=4, linestyle="-"),
    # Line2D([0], [0], color=colors["Green"], lw=4, linestyle="-")
]

labels = [
    "Incidence thresholds",
    "Aachen",
    "Flensburg",
    "Heilbronn"
]
fig, ax = plt.subplots(figsize=(3, 0.8))
ax.axis("off")

ax.legend(
    handles=legend_elements,
    labels=labels,
    loc="center",
    ncol=1,
    frameon=True,
    handler_map={tuple: HandlerTuple(ndivide=None)}
)

plt.savefig(
    os.path.join(save_dir, "nuts3_legend_thresholds_optimization.png"),
    dpi=dpi,
    bbox_inches="tight",
    pad_inches=0.02
)
