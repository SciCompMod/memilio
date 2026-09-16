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
    num_rows = len(sorted_tinit)

    images = []
    for i, tinit in enumerate(sorted_tinit):
        img_path = os.path.join(base_dir, f"tinit={tinit}", filename)
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

    for row, (img, tinit) in enumerate(zip(images, sorted_tinit)):
        ax = fig.add_subplot(gs[row, 0])
        ax.imshow(img, aspect="auto")
        ax.axis("off")

        widetilde_t = tinit - T
        ax.text(-0.12, 0.9, rf"$\widetilde{{T}}={widetilde_t}$",
                transform=ax.transAxes, ha="left", va="top", fontsize=24)

    if output_path is None:
        output_path = os.path.join(base_dir, f"combined_{filename}")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved combined plot to {output_path}")


def main():
    # Directory containing the tinit=<value> subfolders.
    base_dir = "plots/2026-09-15/S_deriv_forward_dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73_kahan=false_buffer=false/detailed_init_exponential_t0ide=50_tmax=150_finite_diff=4"
    # tinit values to combine.
    tinit_values = [40, 50]
    # Path the combined figure is saved to. Set to None to default to base_dir/combined_<filename>.
    output_path = None
    # Name of the plot file to pick up from each tinit subfolder.
    filename = FILENAME
    # Vertical space between rows. Negative values make the plots overlap.
    hspace = -0.05

    combine_tinit_plots(base_dir, tinit_values,
                        output_path=output_path, filename=filename, hspace=hspace)


if __name__ == "__main__":
    main()
