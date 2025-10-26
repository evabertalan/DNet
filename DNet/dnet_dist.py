import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis.*")

import os
import ast
import argparse
import glob
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import apply_PBC
from MDAnalysis.transformations import wrap
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree


class DNetDist:
    def __init__(self, psf, dcd, output_folder, wrap_dcd=True):

        _, ext = os.path.splitext(psf)
        self.psf = psf
        self.dcd = dcd
        self.u = mda.Universe(self.psf, self.dcd)

        if wrap_dcd:
            self.u.trajectory.add_transformations(wrap(self.u.atoms))

        base = os.path.basename(psf)
        self.base_name, ext = os.path.splitext(base)
        self.output_folder = output_folder

    def _parse_cgraphs_edges(self, graphs_input):
        if os.path.exists(graphs_input):
            with open(graphs_input) as f:
                for line in f.readlines():
                    if line.startswith("List of edges:"):
                        edges = ast.literal_eval(line.split(": ")[-1])
                return edges

        else:
            raise FileNotFoundError(f"The {graphs_input} file does not exist.")

    def _parse_cgraphs_nodes(self, graphs_input):
        try:
            with open(graphs_input) as f:
                for line in f.readlines():
                    if line.startswith("List of nodes:"):
                        nodes = ast.literal_eval(line.split(": ")[-1])
            return nodes
        except (FileNotFoundError, ValueError, SyntaxError) as e:
            print(f"Error reading cgraphs file: {e}")
            return None

    def calculate_distances(
        self,
        graphs_input,
        max_water_distance=3.5,
        step_size=1,
        start=None,
        stop=None,
        selection="",
    ):
        self.step_size = step_size

        self.edges = self._parse_cgraphs_edges(graphs_input)
        self.nodes = self._parse_cgraphs_nodes(graphs_input)
        self.max_water_distance = max_water_distance

        if selection:
            self.u.select_atoms(selection)

        bonding_groups = (
            "OH2",
            "OW",
            "NE",
            "NH1",
            "NH2",
            "ND2",
            "SG",
            "NE2",
            "ND1",
            "NZ",
            "OG",
            "OG1",
            "NE1",
            "OH",
            "OE1",
            "OE2",
            "N16",
            "OD1",
            "OD2",
            "SD",
        )
        selection_bonding_group = " or name ".join(bonding_groups)
        water_definition = "resname TIP3 and name OH2"

        box = self.u.dimensions

        if len(self.edges[0][0].split("-")) == 3:
            # if Bridge analysis was done with residuewise=True, for the distance calculation we take the CA atom
            # not recommended for distance calumniation, using residuewise=False gives the distance we are interested in
            e1 = [
                self.u.select_atoms(
                    f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name CA"
                )
                for edge in self.edges
            ]
            e2 = [
                self.u.select_atoms(
                    f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name CA"
                )
                for edge in self.edges
            ]

            n = [
                self.u.select_atoms(
                    f"segid {node.split('-')[0]} and resname {node.split('-')[1]} and resid {node.split('-')[2]} and name CA"
                )
                for node in self.nodes
            ]

        else:  # for residuewise=False, distance can be accurately calculated between the actual H-bonding atom pairs
            e1 = [
                self.u.select_atoms(
                    f"segid {edge[0].split('-')[0]} and resname {edge[0].split('-')[1]} and resid {edge[0].split('-')[2]} and name {edge[0].split('-')[3]}"
                )
                for edge in self.edges
            ]
            e2 = [
                self.u.select_atoms(
                    f"segid {edge[1].split('-')[0]} and resname {edge[1].split('-')[1]} and resid {edge[1].split('-')[2]} and name {edge[1].split('-')[3]}"
                )
                for edge in self.edges
            ]

            n = [
                self.u.select_atoms(
                    f"segid {node.split('-')[0]} and resname {node.split('-')[1]} and resid {node.split('-')[2]} and name {node.split('-')[3]}"
                )
                for node in self.nodes
            ]

            sidechain_atoms = [
                self.u.select_atoms(
                    f"(segid {node.split('-')[0]} and resname {node.split('-')[1]} and resid {node.split('-')[2]}) and (name {selection_bonding_group})"
                )
                for node in self.nodes
            ]

            water_atoms = self.u.select_atoms(water_definition)

            self.distances_over_time = []
            self.water_around_group_over_time = []
            self.water_around_total_res_over_time = []
            self.frames = []

            for ts in self.u.trajectory[slice(start, stop, self.step_size)]:
                self.frames.append(ts.frame)
                box = self.u.dimensions

                water_pos = water_atoms.positions.copy()

                res_positions = np.concatenate([res.positions for res in n])
                dists = distance_array(res_positions, water_pos, box=box)
                water_around_group = (
                    (dists < self.max_water_distance).sum(axis=1).tolist()
                )

                tree = cKDTree(water_pos, boxsize=box[:3])
                water_around_total_res = []

                for group in sidechain_atoms:
                    group_pos = group.positions
                    neighbor_lists = tree.query_ball_point(
                        group_pos, r=self.max_water_distance
                    )
                    unique_waters = np.unique(np.concatenate(neighbor_lists))
                    water_around_total_res.append(len(unique_waters))

                frame_distances = []
                for res1, res2 in zip(e1, e2):
                    box = self.u.dimensions
                    dist = distance_array(res1.positions, res2.positions, box=box)[0][0]
                    frame_distances.append(dist)

                self.distances_over_time.append(frame_distances)
                self.water_around_group_over_time.append(water_around_group)
                self.water_around_total_res_over_time.append(water_around_total_res)

    def write_results_to_df(self):
        self.distances_df = pd.DataFrame(
            np.array(self.distances_over_time),
            columns=[f"{e[0]} - {e[1]}" for e in self.edges],
            index=self.frames,
        )
        # add frame as index column
        self.distances_df.to_csv(
            f"{self.output_folder}/{self.base_name}_pair_distances.csv", index=True
        )

        self.waters_df = pd.DataFrame(
            np.array(self.water_around_group_over_time),
            columns=self.nodes,
            index=self.frames,
        )
        dist = str(self.max_water_distance).replace(".", "_")
        self.waters_df.to_csv(
            f"{self.output_folder}/{self.base_name}_waters_within_{dist}_of_group.csv",
            index=True,
        )

        self.total_water_df = pd.DataFrame(
            np.array(self.water_around_total_res_over_time),
            columns=self.nodes,
            index=self.frames,
        )
        self.total_water_df.to_csv(
            f"{self.output_folder}/{self.base_name}_total_waters_within_{dist}_of_res.csv",
            index=True,
        )

    def plot_results(self):
        self.distances_df = pd.read_csv(
            f"{self.output_folder}/{self.base_name}_pair_distances.csv"
        )
        fig, axes = plt.subplots(
            nrows=len(self.distances_df.columns), ncols=1, figsize=(20, 80), sharex=True
        )

        for i, column in enumerate(self.distances_df.columns):
            axes[i].plot(self.distances_df.index, self.distances_df[column])
            axes[i].spines["right"].set_visible(False)
            axes[i].spines["top"].set_visible(False)
            axes[i].set_title(f"{column}")
            axes[i].set_ylabel("Distance")

        # Set the x-axis label only on the last subplot
        axes[-1].set_xlabel("Frame")

        # Adjust layout to prevent overlap
        plt.tight_layout()
        fig.savefig("./dist_time_series_per_res.png")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze distances of H-bond graph edges and water molecules around graph nodes calculated with C-Graphs from molecular dynamics trajectories."
    )
    parser.add_argument(
        "psf",
        help="Path to the PSF file (Protein Structure File) required to load the molecular system.",
    )
    parser.add_argument(
        "dcd",
        nargs="+",
        help="Path(s) to the DCD trajectory file(s). You can use wildcard patterns (e.g., '*.dcd') to select multiple files.",
    )
    parser.add_argument(
        "graphs_input",
        help="Path to the _info.txt C-Grpahs file containing graph edges and nodes for distance calculations.",
    )
    parser.add_argument(
        "--output_folder",
        help="Path to the folder where output files (CSV and plots) will be saved.",
    )

    parser.add_argument(
        "--selection",
        help="MDAnalysis selection string to restrict distance calculations to a subset of atoms from the cgraphs input.",
    )

    parser.add_argument(
        "--start",
        type=int,
        help="Starting frame index for trajectory analysis. If not provided, starts from the first frame.",
    )
    parser.add_argument(
        "--stop",
        type=int,
        help="Stopping frame index for trajectory analysis. If not provided, processes until the last frame.",
    )
    parser.add_argument(
        "--step",
        type=int,
        help="Step size for iterating through the trajectory frames. For example, '--step 10' "
        "processes every 10th frame to reduce computation time.",
    )

    parser.add_argument(
        "--max_water_distance",
        default=3.5,
        help="Maximum distance (in Å) within which water molecules are considered for analysis. Default is 3.5 Å.",
    )

    parser.add_argument(
        "--wrap_dcd",
        help="Apply wrapping transformations to keep molecules within the simulation box. Default is True",
    )

    args = parser.parse_args()

    dcd_files = []
    for dcd_file in args.dcd:
        dcd_files += glob.glob(dcd_file)
    dcd_files.sort()
    if not dcd_files:
        raise FileNotFoundError("No valid DCD files found.")

    output_folder = (
        args.output_folder if args.output_folder else os.path.dirname(args.psf)
    )
    os.makedirs(output_folder, exist_ok=True)

    if args.wrap_dcd:
        wrap_dcd = args.wrap_dcd.lower() == "true"
    else:
        wrap_dcd = True
    dist_traj = DNetDist(args.psf, dcd_files, output_folder, wrap_dcd=wrap_dcd)
    dist_traj.calculate_distances(
        graphs_input=args.graphs_input,
        max_water_distance=float(args.max_water_distance),
        step_size=args.step,
        start=args.start,
        stop=args.stop,
        selection=args.selection,
    )
    dist_traj.write_results_to_df()
    # dist_traj.plot_results()


if __name__ == "__main__":
    main()
