import os
import argparse
import glob
from dnet_graphs import DNetGraphs


def MD_conserved_graph(
    psf, dcds, name, target_folder, max_water, occupancy, res_id_label_shift
):
    dnet_graphs = DNetGraphs(
        target_folder=target_folder,
        psf_file=psf,
        dcd_files=dcds,
        sim_name=name,
        plot_parameters={"graph_color": "#666666", "formats": ["png", "eps"]},
    )
    dnet_graphs.calculate_graphs(
        max_water=max_water,
        check_angle=True,
        selection="protein or resname LYR or resname HSE",
        additional_donors=["NH"],
        additional_acceptors=["NH"],
        residuewise=True,
        wrap_dcd=True,
    )
    dnet_graphs.plot_graphs(
        label_nodes=True,
        xlabel="PCA projected membrane plane (Å)",
        ylabel="Membrane normal (Å)",
        occupancy=occupancy,
        color_data=True,
        node_color_selection="protein or resname LYR or resname HSE",
        node_color_map="coolwarm_r",
        res_id_label_shift=19,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Process pKa from molecular dynamics trajectories."
    )
    parser.add_argument("--psf", help="Path to the PSF file")
    parser.add_argument(
        "--dcd",
        nargs="+",
        help="Path to the DCD files. The path can contain regex to select multiple files by a name pattern.",
    )
    parser.add_argument(
        "--workfolder",
    )

    parser.add_argument(
        "--max_water",
    )

    parser.add_argument(
        "--occupancy",
    )

    args = parser.parse_args()

    base = os.path.basename(args.psf)
    base_name, ext = os.path.splitext(base)

    dcd_files = []
    for dcd_file in args.dcd:
        dcd_files += glob.glob(dcd_file)
    dcd_files.sort()

    MD_conserved_graph(
        args.psf,
        dcd_files,
        base_name,
        args.workfolder,
        int(args.max_water),
        float(args.occupancy),
        int(args.res_id_label_shift),
    )


if __name__ == "__main__":
    main()
