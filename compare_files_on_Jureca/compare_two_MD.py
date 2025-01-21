import os
import argparse
import glob
from cgraphs import CompareTwo
import sys


def MD_compare_two(psf1, psf2, dcd1, dcd2, name1, name2, color1, color2, target_folder):
    MAX_WATER = 3
    OCCUPANCY = 0.5
    comp = CompareTwo(
        "dcd",
        psf1=psf1,
        psf2=psf2,
        dcd1=dcd1,
        dcd2=dcd2,
        target_folder=target_folder,
        name1=name1,
        name2=name2,
    )
    comp.calculate_graphs(
        graph_type="water_wire", max_water=MAX_WATER, distance=3.5, cut_angle=60
    )
    comp.construct_comparison_objects(occupancy=OCCUPANCY)
    comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=True)
    comp.plot_graph_comparison(color1=color1, color2=color2, label_nodes=False)


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("psf1", help="")
    parser.add_argument("psf2", help="")
    parser.add_argument(
        "workfolder",
    )

    parser.add_argument(
        "--dcd1",
        nargs="+",
        help="Path to the DCD files. The path can contain regex to select multiple files by a name pattern.",
    )
    parser.add_argument(
        "--dcd2",
        nargs="+",
        help="Path to the DCD files. The path can contain regex to select multiple files by a name pattern.",
    )
    parser.add_argument("--color1", help="")
    parser.add_argument("--color2", help="")

    args = parser.parse_args()

    base_1 = os.path.basename(args.psf1)
    base_2 = os.path.basename(args.psf2)
    base_name_1, ext = os.path.splitext(base_1)
    base_name_2, ext = os.path.splitext(base_2)

    dcd_files_1 = []
    for dcd_file in args.dcd1:
        dcd_files_1 += glob.glob(dcd_file)
    dcd_files_1.sort()

    dcd_files_2 = []
    for dcd_file in args.dcd2:
        dcd_files_2 += glob.glob(dcd_file)
    dcd_files_2.sort()

    MD_compare_two(
        args.psf1,
        args.psf2,
        dcd_files_1,
        dcd_files_2,
        base_name_1,
        base_name_2,
        args.color1,
        args.color2,
        args.workfolder,
    )


if __name__ == "__main__":
    main()
