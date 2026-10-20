import os

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

# tinit value of the "detailed" ODE-based initialization used as reference;
# \widetilde{T} = tinit - T is the init window length shown in the label of
# each tinit subfolder.
T = 50

FILENAME = "convergence_all_compartments_max_rel.png"

# Fraction of each plot's height that contains the subplots and their tick
# labels, without the "Time step" axis title and the legend below them.
# Used to crop every row except the bottom one, so that the axis title and
# the legend are only shown once, at the bottom.
LEGEND_CROP_FRAC = 0.695

FIG_WIDTH_INCHES = 14


def _stack_plots(img_paths, labels, output_path, hspace=-0.05, label_x=-0.12, fontsize=24):
    """ Stack plots vertically into one figure with a label at the left of each row.

    The legend is cropped out of every row except the bottom one, and rows are packed tightly
    (using hspace) so that the legend of the bottom plot is the only one visible.

    @param[in] img_paths Paths of the plots to stack, from top to bottom.
    @param[in] labels Label shown at the left of each row, same order as img_paths.
    @param[in] output_path Path the combined figure is saved to.
    @param[in] hspace Vertical space between rows passed to the gridspec. Negative values make
        the plots overlap so that only the bottom row's legend remains visible.
    @param[in] label_x Horizontal position of the row labels in axes coordinates.
    @param[in] fontsize Font size of the row labels.
    """
    num_rows = len(img_paths)

    images = []
    for i, img_path in enumerate(img_paths):
        if not os.path.isfile(img_path):
            raise FileNotFoundError(f"Could not find plot at {img_path}")

        img = mpimg.imread(img_path)
        is_bottom_row = (i == num_rows - 1)
        if not is_bottom_row:
            # Crop out the legend so it is only shown for the bottom plot.
            crop_height = int(img.shape[0] * LEGEND_CROP_FRAC)
            img = img[:crop_height]
        images.append(img)

    # Height ratios based on the (cropped) images' pixel heights so that
    # each row keeps its original aspect ratio despite the rows having
    # different heights.
    height_ratios = [img.shape[0] for img in images]
    img_width = images[0].shape[1]
    fig_height_inches = FIG_WIDTH_INCHES * sum(height_ratios) / img_width

    fig = plt.figure(figsize=(FIG_WIDTH_INCHES, fig_height_inches))
    gs = fig.add_gridspec(
        num_rows, 1, height_ratios=height_ratios, hspace=hspace)

    for row, (img, label) in enumerate(zip(images, labels)):
        ax = fig.add_subplot(gs[row, 0])
        ax.imshow(img, aspect="auto")
        ax.axis("off")

        ax.text(label_x, 0.9, label,
                transform=ax.transAxes, ha="left", va="top", fontsize=fontsize)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved combined plot to {output_path}")


def combine_tinit_plots(base_dir, tinit_values, output_path=None, filename=FILENAME, hspace=-0.05):
    """ Stack the convergence plots of several tinit subfolders vertically into one figure.

    Rows are ordered from the highest to the lowest value of \\widetilde{T} = tinit - T,
    top to bottom. The legend is cropped out of every row except the bottom one, and rows
    are packed tightly (using hspace) so that the legend of the bottom plot is the only
    one visible.

    @param[in] base_dir Directory that contains one subfolder "tinit=<value>" per tinit value.
    @param[in] tinit_values List of tinit values to combine (order does not matter, they are
        sorted by \\widetilde{T} internally).
    @param[in] output_path Path the combined figure is saved to. Defaults to
        base_dir/combined_<filename>.
    @param[in] filename Name of the plot file to pick up from each tinit subfolder.
    @param[in] hspace Vertical space between rows passed to the gridspec. Negative values make
        the plots overlap so that only the bottom row's legend remains visible.
    """
    # Highest \widetilde{T} at the top, lowest at the bottom.
    sorted_tinit = sorted(
        tinit_values, key=lambda tinit: tinit - T, reverse=True)

    img_paths = [os.path.join(base_dir, f"tinit={tinit}", filename)
                 for tinit in sorted_tinit]
    labels = [rf"$\widetilde{{T}}={tinit - T}$" for tinit in sorted_tinit]

    if output_path is None:
        output_path = os.path.join(base_dir, f"combined_{filename}")
    _stack_plots(img_paths, labels, output_path,
                 hspace=hspace, label_x=-0.12, fontsize=16)


def combine_fd_order_plots(base_dir, fd_orders, output_path=None, filename=FILENAME, hspace=-0.05,
                           subfolder_template="finite_diff={fd_order}"):
    """ Stack the convergence plots of several FD orders vertically into one figure.

    Rows are ordered from the lowest to the highest FD order, top to bottom. The legend is
    cropped out of every row except the bottom one, and rows are packed tightly (using hspace)
    so that the legend of the bottom plot is the only one visible.

    @param[in] base_dir Directory that contains one subfolder per FD order, named according to
        subfolder_template.
    @param[in] fd_orders List of FD orders to combine (order does not matter, they are sorted
        internally).
    @param[in] output_path Path the combined figure is saved to. Defaults to
        base_dir/combined_<filename>.
    @param[in] filename Name of the plot file to pick up from each FD order subfolder.
    @param[in] hspace Vertical space between rows passed to the gridspec. Negative values make
        the plots overlap so that only the bottom row's legend remains visible.
    @param[in] subfolder_template Template (with a {fd_order} placeholder) used to build each
        FD order's subfolder path relative to base_dir. May contain additional path components,
        e.g. "detailed_init_exponential_t0ide=50_tmax=150_finite_diff={fd_order}/tinit=10".
    """
    sorted_fd_orders = sorted(fd_orders)

    img_paths = [os.path.join(base_dir, subfolder_template.format(fd_order=fd_order), filename)
                 for fd_order in sorted_fd_orders]
    labels = [f"FD order {fd_order}" for fd_order in sorted_fd_orders]

    if output_path is None:
        output_path = os.path.join(base_dir, f"combined_{filename}")
    _stack_plots(img_paths, labels, output_path,
                 hspace=hspace, label_x=-0.15, fontsize=22)


def combine_labeled_plots(base_dir, subfolders, labels, output_path=None, filename=FILENAME, hspace=-0.05,
                          label_x=-0.2, fontsize=22):
    """ Stack the convergence plots of several subfolders vertically into one figure, each row with a
    given label.

    The legend is cropped out of every row except the bottom one, and rows are packed tightly (using
    hspace) so that the legend of the bottom plot is the only one visible.

    @param[in] base_dir Directory that contains the subfolders.
    @param[in] subfolders Subfolder paths relative to base_dir, from top to bottom. May contain
        additional path components, e.g. "detailed_init_exponential_t0ide=50_tmax=150_finite_diff=4/tinit=50".
    @param[in] labels Label shown at the left of each row, same order as subfolders.
    @param[in] output_path Path the combined figure is saved to. Defaults to
        base_dir/combined_<filename>.
    @param[in] filename Name of the plot file to pick up from each subfolder.
    @param[in] hspace Vertical space between rows passed to the gridspec. Negative values make
        the plots overlap so that only the bottom row's legend remains visible.
    @param[in] label_x Horizontal position of the row labels in axes coordinates.
    @param[in] fontsize Font size of the row labels.
    """
    if len(subfolders) != len(labels):
        raise ValueError(
            f"Expected one label per subfolder, got {len(labels)} labels for {len(subfolders)} subfolders.")

    img_paths = [os.path.join(base_dir, subfolder, filename)
                 for subfolder in subfolders]

    if output_path is None:
        output_path = os.path.join(base_dir, f"combined_{filename}")
    _stack_plots(img_paths, labels, output_path,
                 hspace=hspace, label_x=label_x, fontsize=fontsize)


def main():

    # Name of the plot file to pick up from each tinit subfolder.
    filename = FILENAME
    # Vertical space between rows. Negative values make the plots overlap.
    hspace = -0.05

    # # Directory containing the tinit=<value> subfolders.
    # base_dir = "plots/2026-09-24/forward_ta_dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73_kahan=false_forwardfd=true/detailed_init_exponential_t0ide=50_tmax=150_finite_diff=4"
    # # tinit values to combine.
    # tinit_values = [40, 50]
    # # Path the combined figure is saved to. Set to None to default to base_dir/combined_<filename>.
    # output_path = None

    # combine_tinit_plots(base_dir, tinit_values,
    #                     output_path=output_path, filename=filename, hspace=hspace)

    # # Directory containing the FD order subfolders.
    # fd_base_dir = "plots/2026-10-11/dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73_kahan=false_forwardfd=true"
    # # FD orders to combine.
    # fd_orders = [1, 2, 3]
    # # Template mapping each FD order to its subfolder, relative to fd_base_dir.
    # subfolder_template = "detailed_init_exponential_t0ide=50_tmax=150_finite_diff={fd_order}/tinit=10"

    # combine_fd_order_plots(fd_base_dir, fd_orders, filename=filename, hspace=hspace,
    #                        subfolder_template=subfolder_template)

    # Directory containing the fine offsets subfolders.
    fineoffsets_base_dir = "plots/2026-10-16/dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73_kahan=false_forwardfd=true"
    # Subfolders relative to fineoffsets_base_dir, from top to bottom, with the corresponding labels.
    fineoffsets_subfolders = [
        "detailed_init_exponential_t0ide=50_tmax=150_finite_diff=4_fineoffsets=0-0-0-0_fineflows=true/tinit=50",
        "detailed_init_exponential_t0ide=50_tmax=150_finite_diff=4_fineoffsets=3-3-3-3_fineflows=true/tinit=50"]
    fineoffsets_labels = [r"$\Delta t_0=\Delta t$",
                          r"$\Delta t_0=10^{-3}\,\Delta t$"]

    combine_labeled_plots(fineoffsets_base_dir, fineoffsets_subfolders, fineoffsets_labels, filename=filename,
                          hspace=hspace)


if __name__ == "__main__":
    main()
