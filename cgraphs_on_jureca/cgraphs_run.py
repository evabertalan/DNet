import os
import argparse
import glob
from cgraphs import ConservedGraph


def MD_conserved_graph(psfs, dcds, names, target_folder):
    MAX_WATER = 3
    OCCUPANCY = 0.1
    c_dcd = ConservedGraph(
        type_option="dcd",
        dcd_files=dcds,
        psf_files=psfs,
        sim_names=names,
        target_folder=target_folder,
        pdb_root_folder=target_folder,
        plot_parameters={"graph_color": "#666666", "formats": ["png", "eps"]},
    )
    c_dcd.calculate_graphs(
        graph_type="water_wire",
        max_water=MAX_WATER,
        check_angle=True,
        selection="protein or resname LYR or resname HSE",
        additional_donors=["NH"],
        additional_acceptors=["NH"],
    )
    c_dcd.get_conserved_graph(conservation_threshold=0.8, occupancy=OCCUPANCY)
    c_dcd.plot_graphs(
        label_nodes=True,
        xlabel="PCA projected membrane plane",
        ylabel="Membrane normal ($\AA$)",
        color_data=True,
        node_color_selection="protein or resname LYR or resname HSE",
        node_color_map="coolwarm_r",
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
    print(base_name)

    dcd_files = []
    for dcd_file in args.dcd:
        dcd_files += glob.glob(dcd_file)
    dcd_files.sort()

    MD_conserved_graph([args.psf], [dcd_files], [base_name], args.workfolder)


if __name__ == "__main__":
    main()
