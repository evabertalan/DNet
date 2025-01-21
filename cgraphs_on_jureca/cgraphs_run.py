import os
import argparse
import glob
from proteingraphanalyser import ProteinGraphAnalyser


def MD_conserved_graph(psf, dcds, name, target_folder):
    MAX_WATER = 3
    OCCUPANCY = 0.4
    c_dcd = ProteinGraphAnalyser(
        target_folder=target_folder,
        psf_file=psf,
        dcd_files=dcds,
        sim_name=name,
        plot_parameters={"graph_color": "#666666", "formats": ["png", "eps"]},
    )
    c_dcd.calculate_graphs(
        max_water=MAX_WATER,
        check_angle=True,
        selection="protein or resname LYR or resname HSE",
        additional_donors=["NH"],
        additional_acceptors=["NH"],
        residuewise=True,
        wrap_dcd=True,
    )
    c_dcd.plot_graphs(
        label_nodes=True,
        xlabel="PCA projected membrane plane",
        ylabel="Membrane normal (Å)",
        occupancy=OCCUPANCY,
        color_data=True,
        node_color_selection="protein or resname LYR or resname HSE",
        node_color_map="coolwarm_r",
        res_id_label_shift=19,
    )
    # c_dcd.plot_graphs(label_nodes=True, xlabel='PCA projected membrane plane', ylabel='Membrane normal ($\AA$)')


def main():
    parser = argparse.ArgumentParser(
        description="Process pKa from molecular dynamics trajectories."
    )
    parser.add_argument("psf", help="Path to the PSF file")
    parser.add_argument(
        "dcd",
        nargs="+",
        help="Path to the DCD files. The path can contain regex to select multiple files by a name pattern.",
    )
    parser.add_argument(
        "workfolder",
    )

    args = parser.parse_args()

    base = os.path.basename(args.psf)
    base_name, ext = os.path.splitext(base)

    dcd_files = []
    for dcd_file in args.dcd:
        dcd_files += glob.glob(dcd_file)
    dcd_files.sort()

    MD_conserved_graph(args.psf, dcd_files, base_name, args.workfolder)


if __name__ == "__main__":
    main()
