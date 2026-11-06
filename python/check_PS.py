# Copyright © 2026 GSI Helmholtzzentrum für Schwerionenforschung gmbh
# <https://gsi.de/>
#
# Everyone is permitted to copy and distribute verbatim copies of this license
# document, but changing it is not allowed.
#
# Preamble
# The GNU Affero General Public License (version 3) is a free, copyleft
# license for software and other kinds of works. The licenses for most
# software and other practical works are designed to take away your freedom
# to share and change the works. By contrast, the GNU Affero General Public
# License is intended to guarantee your freedom to share and change all
# versions of a program--to make sure it remains free software for all its
# users. We, GSI Helmholtzzentrum für Schwerionenforschung, use the GNU
# Affero General Public License for most of our software; it applies also to
# any other work released this way by its authors. You can apply it to your
# programs, too.

"""
Check the phase space.
Required original, reshoot ps and a path on where to store the plots.
"""

import sys
from pathlib import Path
import uproot
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams["figure.figsize"] = (10, 10)
import numpy as np
import os


def import_root_file(root_file_path: Path, hist_name, flow=False):
    file = uproot.open(root_file_path)
    hist = file[hist_name]
    x = hist.axis().centers(flow=flow)
    y = hist.values(flow=flow)
    return x, y


def get_data_for_plot(root_file_path, hist_name, flow=False):
    file = uproot.open(root_file_path)
    hist = file[hist_name]
    edges = hist.axis().edges(flow=flow)
    widths = hist.axis().widths(flow=flow)
    y = hist.values(flow=flow)
    return edges, y / widths


def generate_config_dict(
    position_original_index: int,
    dimension_original_index: int,
    quantity_original_index: int,
    quantity_reshoot_index: int,
):
    """
    @brief
    Generate a dictionary file used to selec the root keys
    requires the index of the selected vectors
    @param position_original_index: index of the position list [100, 1500, 2000]
    @param dimensions_PS_original_list: index of the
    """
    # define constants
    positions_PS_original_list = [
        "100",
        "1500",
        "2000",
    ]
    dimensions_PS_original_list = [
        "",
        "wide",
        "very_wide",
    ]

    quantity_PS_original_list = [
        "Theta",
        "Phi",
        "Ekin",
        "X",
        "Y",
    ]

    quantity_PS_reshoot_list = [
        "Theta",
        "Phi",
        "EkinU",
        "X",
        "Y",
    ]
    plot_units_list_y = [
        "counts [1/rad]",
        "counts [1/rad]",
        "counts [u/MeV]",
        "counts [1/mm]",
        "counts [1/mm]",
    ]

    plot_units_list_x = [
        "Theta [rad]",
        "Phi [rad]",
        "Ek [MeV/u]",
        "x pos [mm]",
        "y pos [mm]",
    ]
    config_dic = {
        "original_config": {
            "position": positions_PS_original_list[position_original_index],
            "dimension": dimensions_PS_original_list[dimension_original_index],
            "quantity": quantity_PS_original_list[quantity_original_index],
        },
        "reshoot_config": {
            "quantity": quantity_PS_reshoot_list[quantity_reshoot_index]
        },
        "units_x": plot_units_list_x[quantity_original_index],
        "units_y": plot_units_list_y[quantity_original_index],
    }

    return config_dic


def calculate_chi2(Y1, Y2):
    delta2 = (Y1 - Y2) ** 2
    den = Y1 + Y2  # quadratic sum of errors:
    return delta2 / den


def _generate_key_original_ps(config_dic, z):
    pos = config_dic["original_config"]["position"]
    dim = config_dic["original_config"]["dimension"]
    quan = config_dic["original_config"]["quantity"]
    if dim == "":
        original_PS_key = f"{quan}_{pos}_{z}"
    else:
        original_PS_key = f"{quan}_{dim}_{pos}_{z}"
    return original_PS_key


def generate_key_reshoot_ps(config_dic, z):
    quan = config_dic["reshoot_config"]["quantity"]
    reshoot_PS_key = f"{quan}_{1}_{z}"
    return reshoot_PS_key


def calculate_normalization(original_root, reshoot_root, config_dic, z):
    """
    Calcualte the normalization based on a given histogram.
    Such normalization is requried to scale reshooted PS to the original PS.
    """
    original_root_key = _generate_key_original_ps(config_dic, z)
    reshoot_root_key = generate_key_reshoot_ps(config_dic, z)
    # x_original, y_original = import_root_file(
    #    original_root, original_root_key, flow=True
    # )
    # x_reshoot, y_reshoot = import_root_file(reshoot_root, reshoot_root_key, flow=True)

    counts_original = calculate_total_counts(
        original_root, original_root_key, flow=True
    )
    counts_reshoot = calculate_total_counts(reshoot_root, reshoot_root_key, flow=True)
    # counts_reshoot = np.trapezoid(y_reshoot[1:-1], x_reshoot)
    # counts_original = np.trapezoid(y_original[1:-1], x_original)
    # overflow_orig = y_original[-1]
    # underflow_orig = y_original[0]
    # overflow_res = y_reshoot[-1]
    # underflow_res = y_reshoot[0]

    return counts_original / counts_reshoot


def calculate_total_counts(root_histogram, key, flow=False):
    x, y = import_root_file(root_histogram, key, flow=flow)

    integral = np.trapezoid(y[1:-1], x[1:-1])
    return integral + y[0] + y[-1]
    # return np.sum(y)


def define_particle_list(z_max):
    particle_list = [str(a) for a in range(0, z_max)]
    extra_particles = ["e-", "e+", "gamma", "mu-", "mu+", "pi-", "pi+"]
    for p in extra_particles:
        particle_list.append(p)
    return particle_list


def process_data(
    original_root, reshoot_root, config_dic, z_max, normalization, output_pdf, pbar
):
    particle_list = define_particle_list(z_max)
    for p in particle_list:
        original_root_key = _generate_key_original_ps(config_dic, p)
        reshoot_root_key = generate_key_reshoot_ps(config_dic, p)
        x_original_plot, y_original_plot = get_data_for_plot(
            original_root, original_root_key
        )
        x_reshoot_plot, y_reshoot_plot = get_data_for_plot(
            reshoot_root, reshoot_root_key
        )

        fig, ax = plt.subplots()
        ax.stairs(
            y_reshoot_plot * normalization,
            x_reshoot_plot,
            label="output PS * normalization",
            color="orange",
        )

        ax.stairs(
            y_reshoot_plot,
            x_reshoot_plot,
            color="crimson",
            label="output PS",
        )
        ax.stairs(y_original_plot, x_original_plot, label="input PS", color="navy")

        x_original, y_original = import_root_file(original_root, original_root_key)
        x_reshoot, y_reshoot = import_root_file(reshoot_root, reshoot_root_key)
        original_counts = calculate_total_counts(original_root, original_root_key)
        reshoot_counts = calculate_total_counts(reshoot_root, reshoot_root_key)
        str_to_print = f"original counts = {original_counts:.1f}, reshoot_counts = {reshoot_counts:.1f}"
        ax.text(0.5, 0.8, str_to_print, transform=plt.gca().transAxes)
        str_to_print = f"normalization = {normalization:.3f}, reshoot_counts*normalization = {reshoot_counts * normalization:.1f}"
        ax.text(0.5, 0.75, str_to_print, transform=plt.gca().transAxes)
        plot_units_x = config_dic["units_x"]
        plot_units_y = config_dic["units_y"]

        # means and averages are calculated without overflow, on the data plotted
        norm_reshoot = np.trapezoid(y_reshoot, x_reshoot)
        if norm_reshoot != 0:
            x_mean_reshoot = (
                np.trapezoid(y_reshoot * x_reshoot, x_reshoot) / norm_reshoot
            )
        else:
            x_mean_reshoot = np.nan
        norm_original = np.trapezoid(y_original, x_original)
        if norm_original != 0:
            x_mean_original = (
                np.trapezoid(y_original * x_original, x_original) / norm_original
            )
        else:
            x_mean_original = np.nan

        str_to_print = f"original average = {x_mean_original:.2f} {plot_units_x}, reshoot average = {x_mean_reshoot:.2f} {plot_units_x}"

        ax.text(0.5, 0.7, str_to_print, transform=plt.gca().transAxes)
        ax.axvline(x_mean_original, label="original average", color="green")
        ax.axvline(x_mean_reshoot, label="reshoot average", color="lightgreen")
        ax.grid()
        ax.legend()
        if norm_reshoot != 0:
            ax.set_yscale("log")
        ax.set_xlabel(plot_units_x)
        ax.set_ylabel(plot_units_y)
        title_config = (
            f"original_key={original_root_key}___reshoot_key={reshoot_root_key}"
        )
        ax.set_title(title_config)
        output_pdf.savefig(fig)
        plt.close(fig)
        pbar.update(1)


def main():
    if len(sys.argv) != 6:
        exp_msg = "Error: please provide the following arguments:\n"
        exp_msg += "\t<path_to_input.root>: path to the input root file to be converted using the G4_PS converter\n"
        exp_msg += "\t<path_to_output.root>: path to the output root file generated after the simulations using the converted PS via the macros\n"
        exp_msg += "\t<path_to_output_folder>: path to an output folder where the plots are saved\n"
        exp_msg += "\t<user_configuration>: index [0,2] of the configuration selected: 0 = 100, 1 = 1500, 2 = 2000\n"
        exp_msg += "\t<user_dimensions>: index [0,2] of the dimensions selected: 0 = normal, 1 = wide, 2 = very_wide\n"
        exp_msg += "Usage example: \n"
        exp_msg += "python3 check_PS.py /home/questo/data/sim_d/GCR_2025/G4_GCR_PS/PSmod_1000_34789751/PSmod_1000_34789751.root /home/questo/code/geant/G4_phase_space_tester/running/PSmod_1000_34789751/PSmod_1000_34789751.root /tmp/francisbaldo/ 0 0 \n"
        raise AttributeError(exp_msg)

    z_max = 27
    original_root = Path(sys.argv[1])
    if not original_root.exists():
        exp_msg = f"error: cannot find original file {original_root}"
        raise IOError(exp_msg)
    reshoot_root = Path(sys.argv[2])
    if not reshoot_root.exists():
        exp_msg = f"error: cannot find reshoot file {reshoot_root}"
        raise IOError(exp_msg)

    output_path = Path(sys.argv[3])
    if not output_path.exists():
        exp_msg = f"error: output directory do not exists ({output_path})"
        raise IOError(exp_msg)

    user_configuration = int(sys.argv[4])
    user_dimensions = int(sys.argv[5])

    config_dic = generate_config_dict(user_configuration, user_dimensions, 2, 2)
    normalization = calculate_normalization(original_root, reshoot_root, config_dic, 1)
    print(f"normalization = {normalization}")
    n_of_keys = 5
    total_intructions = n_of_keys * (z_max + 7)

    with PdfPages(f"{output_path}/validation.pdf") as output_pdf:
        with tqdm(total=total_intructions, desc="Processing", position=0) as pbar:
            for i in range(n_of_keys):  ## loop over the quantityes
                config_dic = generate_config_dict(
                    user_configuration, user_dimensions, i, i
                )
                # print(f"Processing key = {config_dic['original_config']['quantity']}")
                # print(config_dic)
                process_data(
                    original_root,
                    reshoot_root,
                    config_dic,
                    z_max,
                    normalization,
                    output_pdf,
                    pbar,
                )


if __name__ == "__main__":
    main()
