##Copyright © 2024 GSI Helmholtzzentrum für Schwerionenforschung gmbh
#<https://gsi.de/>
# 
#Everyone is permitted to copy and distribute verbatim copies of this license
#document, but changing it is not allowed.
# 
#Preamble
#The GNU General Public License is a free, copyleft license for software and
#other kinds of works. The licenses for most software and other practical
#works are designed to take away your freedom to share and change the works.
#By contrast, the GNU General Public License is intended to guarantee your
#freedom to share and change all versions of a program--to make sure it
#remains free software for all its users. We, GSI Helmholtzzentrum für
#Schwerionenforschung, use the GNU General Public License for most of our
#software; it applies also to any other work released this way by its authors.
#You can apply it to your programs, too.

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


def import_root_file(root_file_path: Path, hist_name):
    file = uproot.open(root_file_path)
    hist = file[hist_name]
    x = hist.axis().centers()
    y = hist.values()
    return x, y


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

    plot_units_list = [
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
        "units": plot_units_list[quantity_original_index],
    }

    return config_dic


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
    x_original, y_original = import_root_file(original_root, original_root_key)
    x_reshoot, y_reshoot = import_root_file(reshoot_root, reshoot_root_key)

    counts_reshoot = np.sum(y_reshoot)
    counts_original = np.sum(y_original)

    return counts_original / counts_reshoot


def process_data(
    original_root, reshoot_root, config_dic, z_max, normalization, output_pdf, pbar
):
        for z in range(0, z_max):
            original_root_key = _generate_key_original_ps(config_dic, z)
            reshoot_root_key = generate_key_reshoot_ps(config_dic, z)
            x_original, y_original = import_root_file(original_root, original_root_key)
            x_reshoot, y_reshoot = import_root_file(reshoot_root, reshoot_root_key)

            fig, ax = plt.subplots()
            ax.plot(
                x_reshoot,
                y_reshoot * normalization,
                label="output PS * normalization",
                marker="d",
                color="orange",
            )

            ax.plot(
                x_reshoot,
                y_reshoot,
                marker="p",
                color="crimson",
                label="output PS",
            )
            ax.plot(x_original, y_original, label="input PS", marker="x", color="navy")
            str_to_print = f"original counts = {np.sum(y_original):.1f}, reshoot_counts = {np.sum(y_reshoot):.1f}"
            ax.text(0.5, 0.8, str_to_print, transform=plt.gca().transAxes)
            str_to_print = f"normalization = {normalization:.3f}, reshoot_counts*normalization = {np.sum(y_reshoot)*normalization:.1f}"
            ax.text(0.5, 0.75, str_to_print, transform=plt.gca().transAxes)
            plot_units = config_dic["units"]

            norm_reshoot = np.trapezoid(y_reshoot, x_reshoot)
            x_mean_reshoot = np.trapezoid(y_reshoot * x_reshoot, x_reshoot) / norm_reshoot
            norm_original = np.trapezoid(y_original, x_original)
            x_mean_original = np.trapezoid(y_original * x_original, x_original) / norm_original
            str_to_print = f"original average = {x_mean_original:.2f} {plot_units}, reshoot average = {x_mean_reshoot:.2f} {plot_units}"

            ax.text(0.5, 0.7, str_to_print, transform=plt.gca().transAxes)
            ax.axvline(x_mean_original, label="original average", color="green")
            ax.axvline(x_mean_reshoot, label="reshoot average", color="lightgreen")
            ax.grid()
            ax.legend()
            ax.set_yscale("log")
            ax.set_xlabel(plot_units)
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
    total_intructions = n_of_keys*z_max

    with PdfPages(f'{output_path}/validation.pdf') as output_pdf:
        with tqdm(total=total_intructions, desc="Processing", position = 0) as pbar:
            for i in range(n_of_keys): ## loop over the quantityes
                config_dic = generate_config_dict(user_configuration, user_dimensions, i, i)
                #print(f"Processing key = {config_dic['original_config']['quantity']}")
                #print(config_dic)
                process_data(
                    original_root, reshoot_root, config_dic, z_max, normalization, output_pdf, pbar
                )


if __name__ == "__main__":
    main()
