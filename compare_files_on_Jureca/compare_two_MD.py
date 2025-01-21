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
    dcd_1 = sys.argv[1].split()

    # The second positional argument is the second list of files (list2)
    dcd_2 = sys.argv[2].split()

    parser.add_argument("psf_1", help="")
    parser.add_argument("psf_2", help="")
    parser.add_argument("color1", help="")
    parser.add_argument("color2", help="")
    parser.add_argument(
        "workfolder",
    )

    args = parser.parse_args()

    base_1 = os.path.basename(args.psf_1)
    base_2 = os.path.basename(args.psf_2)
    base_name_1, ext = os.path.splitext(base_1)
    base_name_2, ext = os.path.splitext(base_2)

    dcd_files_1 = []
    for dcd_file in dcd_1:
        dcd_files_1 += glob.glob(dcd_file)
    dcd_files_1.sort()

    dcd_files_2 = []
    for dcd_file in dcd_2:
        dcd_files_2 += glob.glob(dcd_file)
    dcd_files_2.sort()

    MD_compare_two(
        args.psf_1,
        args.psf_2,
        dcd_files_1,
        dcd_files_2,
        base_name_1,
        base_name_2,
        color1,
        color2,
        args.workfolder,
    )


if __name__ == "__main__":
    main()
