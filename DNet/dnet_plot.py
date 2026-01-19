import ast
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os


class DNetPlot:
    def __init__(
        self,
        graphs_info_txt,
        pKas_for_frame_csv,
        pair_distances_csv,
        water_within_csv,
        total_water_within_csv,
        plot_folder,
        sim_name,
        step=1,
    ):
        self.plot_folder = plot_folder
        self.sim_name = sim_name
        self.graph_nodes = self._get_H_bond_nodes(graphs_info_txt)
        pKa_nodes = [("-").join(node.split("-")[0:3]) for node in self.graph_nodes]

        if pKas_for_frame_csv and Path(pKas_for_frame_csv).exists():
            pKas = pd.read_csv(pKas_for_frame_csv, index_col="frame")
            pKas = self._handle_HSE(self.graph_nodes, pKas)
            self.pKas = pKas.loc[:, pKas.columns.isin(pKa_nodes)][::step]
            self.pKas.to_csv(
                Path(self.plot_folder, f"pKa_nodes_per_freame_{sim_name}.csv")
            )
        else:
            print(
                "No pKa time series file was found. Module continues without plotting pKa values."
            )
            self.pKas = None

        self.distances = pd.read_csv(pair_distances_csv, index_col=0)
        self.distances.to_csv(
            Path(self.plot_folder, f"edge_distances_per_freame_{sim_name}.csv")
        )

        self.water_numbers = pd.read_csv(water_within_csv, index_col=0)
        self.water_numbers.to_csv(
            Path(self.plot_folder, f"water_aroun_atom_per_freame_{sim_name}.csv")
        )

        total_water_around_res = pd.read_csv(total_water_within_csv, index_col=0)
        total_water_around_res.columns = (
            total_water_around_res.columns.str.split("-").str[:3].str.join("-")
        )
        self.total_water_around_res = total_water_around_res.loc[
            :, ~total_water_around_res.columns.duplicated()
        ]
        self.total_water_around_res.to_csv(
            Path(
                self.plot_folder,
                f"total_water_around_res_per_freame_{sim_name}.csv",
            )
        )

    def _handle_HSE(self, graph_nodes, pKa_df):
        for node in graph_nodes:
            seg, res_name, res_id, *rest = node.split("-")

            his_col = f"{seg}-HIS-{res_id}"
            if res_name == "HSE" and his_col in pKa_df.columns:
                pKa_df.rename(
                    columns={
                        col: f"{seg}-HSE-{res_id}" if col == his_col else col
                        for col in pKa_df.columns
                    },
                    inplace=True,
                )
            if res_name == "HSD" and his_col in pKa_df.columns:
                pKa_df.rename(
                    columns={
                        col: f"{seg}-HSD-{res_id}" if col == his_col else col
                        for col in pKa_df.columns
                    },
                    inplace=True,
                )
        return pKa_df

    def _get_H_bond_nodes(self, file):
        nodes = None
        try:
            with open(file, "r") as f:
                for line in f:
                    if line.startswith("List of nodes:"):
                        nodes = ast.literal_eval(line.split(": ")[-1])
                        break
        except (FileNotFoundError, ValueError) as e:
            print(f"Error processing file {file}: {e}")
        return nodes

    def _ax_util(self, ax, title=None, xlabel=None, ylabel=None, only_integers=False):
        title_fs = 36
        label_fs = 34
        tick_fs = 32
        ax.set_title(title, fontsize=title_fs)
        ax.set_xlabel(xlabel, fontsize=label_fs)
        ax.set_ylabel(ylabel, fontsize=label_fs)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.tick_params(axis="x", labelsize=tick_fs)
        ax.tick_params(axis="y", labelsize=tick_fs)
        if only_integers:
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        return ax

    def _shift_resid_index(self, node_names, shift=0):
        node_names = node_names.split(" - ")
        updated_node_names = []
        for n in node_names:
            parts = n.split("-")
            parts[2] = str(int(parts[2]) + shift)
            n_joined = ("-").join(parts)
            updated_node_names.append(n_joined)
        return (" - ").join(updated_node_names)

    def calculate_PMF(self, distances, T=300, num_bins=100, min_th=0.001):
        # Boltzmann constant in kcal/(mol⋅K) (kB = 1.380649e-23 J/K)
        # https://en.wikipedia.org/wiki/Boltzmann_constant
        kB_kcalmol = 1.987204259e-3
        self.T = T  # Temperature in Kelvin
        kT_kcalmol = kB_kcalmol * self.T

        # Create histogram (frequency counts) density=True, means it is normalized
        hist, bin_edges = np.histogram(distances, bins=num_bins, density=True)

        # The midpoints of the bins represent the reaction coordinate
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # set minimum probability threshold to avoid rare states on the PMF
        threshold = min_th * np.sum(
            hist
        )  # TODO: maybe define it in a more intuitive way
        # maybe use
        # Calculate the probability distribution (P(x))
        # Hist is already normalized by the 'density=True' --> hist = P(x)
        P_x = hist[np.where(hist > threshold)[0]]
        bin_centers = bin_centers[np.where(hist > threshold)[0]]

        # Avoid log of zero (bins with zero probability)
        # replacing zeros with a small number
        P_x[P_x == 0] = 1e-10

        # Compute the PMF W(x) = -kT ln P(x) (in units of kcal/mol)
        PMF = -kT_kcalmol * np.log(P_x)

        # Normalize the PMF to set the minimum to 0
        PMF -= np.min(PMF)

        return bin_centers, PMF

    def count_substates(self, distances, num_bins=100):
        hist, bin_edges = np.histogram(distances, bins=100, density=True)
        if hist[0] > 0.1:
            # this is to pad values, if the histogram starts with a sharp high peak
            # the smoothing algorithm takes in consideration the neighboring bins
            # if the first bin has a high peak, will be missed
            hist = np.append(np.linspace(0, 0.01, 6), hist)
            bin_edges = np.append(np.zeros(6), bin_edges)
        # Simple moving average
        smoothed_hist = np.convolve(hist, np.ones(9) / 9, mode="same")

        peaks, properties = find_peaks(smoothed_hist, prominence=0.03)
        return len(peaks)

    def create_combined_plots(
        self,
        frame_to_time=100,
        pmf_last_nth_frames=20000,
        plot_formats=["png"],
        res_id_label_shift=0,
    ):

        pKa_color = "#227351"
        total_water_color = "#8f5038"
        water_color = "#cc7f27"
        dist_color = "#335080"
        text_fs = 32

        num_bins = 100

        substates_per_node = []
        unique_hbond_per_res = []
        substates_per_edge = []

        for graph_node in self.graph_nodes:
            graph_node = ("-").join(graph_node.split("-")[0:3])
            total_number_of_states = 0

            if self.pKas and graph_node in self.pKas.columns:
                pKa_column = graph_node
            else:
                pKa_column = None

            dist_columns = [
                col
                for col in self.distances.columns
                if any(
                    ("-").join(part.split("-")[0:3]) == graph_node
                    for part in col.split(" - ")
                )
            ]

            # for plotting we need to sort the columns consistently:
            rearranged_dist_columns = [
                (
                    edge
                    if graph_node == ("-").join(edge.split(" - ")[0].split("-")[0:3])
                    else " - ".join(edge.split(" - ")[::-1])
                )
                for edge in dist_columns
            ]
            distcance_columns = sorted(
                rearranged_dist_columns,
                key=lambda pair: int(pair.split(" - ")[1].split("-")[2]),
            )

            water_columns = sorted(
                [
                    col
                    for col in self.water_numbers.columns
                    if any(
                        ("-").join(part.split("-")[0:3]) == graph_node
                        for part in col.split(" - ")
                    )
                ]
            )
            pka_rows = 1 if pKa_column else 0
            row = 0

            fig, ax = plt.subplots(
                nrows=pka_rows + 1 + len(distcance_columns) + len(water_columns),
                ncols=4,
                figsize=(40, (2 + len(distcance_columns) + len(water_columns)) * 8),
                gridspec_kw={"width_ratios": [2, 1, 1, 1]},
            )
            if pKa_column:
                ax[row, 0].plot(
                    self.pKas.index / frame_to_time,
                    self.pKas[pKa_column],
                    color=pKa_color,
                )
                ax[row, 0] = self._ax_util(
                    ax[0, 0],
                    title=self._shift_resid_index(pKa_column, res_id_label_shift),
                    xlabel="Time (ns)",
                    ylabel="pKa",
                )

                ax[row, 1].hist(
                    self.pKas[pKa_column],
                    bins=num_bins,
                    edgecolor=pKa_color,
                    color=pKa_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[row, 1] = self._ax_util(ax[0, 1], xlabel="Frequency", ylabel="pKa")

                last_x_pKa = self.pKas[
                    self.pKas.index > max(self.pKas.index) - pmf_last_nth_frames
                ]
                ax[row, 2].hist(
                    last_x_pKa[pKa_column],
                    bins=num_bins,
                    edgecolor=pKa_color,
                    color=pKa_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[row, 2] = self._ax_util(ax[0, 2], xlabel="Frequency", ylabel="pKa")
                ax[row, 2].text(
                    0.95,
                    0.95,
                    f"last {pmf_last_nth_frames/frame_to_time:.0f} ns",
                    horizontalalignment="right",  # Align text to the right
                    verticalalignment="top",  # Align text to the top
                    transform=ax[0, 2].transAxes,  # Use normalized coordinates (0 to 1)
                    fontsize=text_fs,
                )
                ax[row, 3].axis("off")
                row += 1

            ax[row, 0].plot(
                self.total_water_around_res.index / frame_to_time,
                self.total_water_around_res[graph_node],
                color=total_water_color,
            )
            ax[row, 0] = self._ax_util(
                ax[row, 0],
                title=f"total # of waters within 3.5 Å of {self._shift_resid_index(graph_node, res_id_label_shift)}",
                xlabel="Time (ns)",
                ylabel="#waters",
                only_integers=True,
            )

            ax[row, 1].hist(
                self.total_water_around_res[graph_node],
                bins=num_bins,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                rasterized=True,
                orientation="horizontal",
            )
            ax[row, 1] = self._ax_util(
                ax[row, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
            )

            last_x_total_water = self.total_water_around_res[
                self.total_water_around_res.index
                > max(self.total_water_around_res.index) - pmf_last_nth_frames
            ]
            ax[row, 2].hist(
                last_x_total_water[graph_node],
                bins=num_bins,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                rasterized=True,
                orientation="horizontal",
            )
            ax[row, 2] = self._ax_util(
                ax[row, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
            )
            ax[row, 2].text(
                0.95,
                0.95,
                f"last {pmf_last_nth_frames/frame_to_time:.0f} ns",
                horizontalalignment="right",  # Align text to the right
                verticalalignment="top",  # Align text to the top
                transform=ax[1, 2].transAxes,  # Use normalized coordinates (0 to 1)
                fontsize=text_fs,
            )
            ax[row, 3].axis("off")

            for k, wat_col in enumerate(water_columns):
                x = pka_rows + 1 + k

                ax[x, 0].plot(
                    self.water_numbers.index / frame_to_time,
                    self.water_numbers[wat_col],
                    color=water_color,
                )
                ax[x, 0] = self._ax_util(
                    ax[x, 0],
                    title=f"# of waters within 3.5 Å of {self._shift_resid_index(wat_col, res_id_label_shift)}",
                    xlabel="Time (ns)",
                    ylabel="#waters",
                    only_integers=True,
                )

                ax[x, 1].hist(
                    self.water_numbers[wat_col],
                    bins=num_bins,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[x, 1] = self._ax_util(
                    ax[x, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
                )

                last_x_water = self.water_numbers[
                    self.water_numbers.index
                    > max(self.water_numbers.index) - pmf_last_nth_frames
                ]
                ax[x, 2].hist(
                    last_x_water[wat_col],
                    bins=num_bins,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[x, 2] = self._ax_util(
                    ax[x, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
                )
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {pmf_last_nth_frames/frame_to_time:.0f} ns",
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=ax[x, 2].transAxes,
                    fontsize=text_fs,
                )
                ax[x, 3].axis("off")

            total_substates = 0
            for j, dist_col in enumerate(distcance_columns):

                og_column_name = (
                    dist_col
                    if dist_col in self.distances.columns
                    else " - ".join(dist_col.split(" - ")[::-1])
                )

                x = 1 + j + pka_rows + len(water_columns)
                ax[x, 0].plot(
                    self.distances.index / frame_to_time,
                    self.distances[og_column_name],
                    color=dist_color,
                )
                ax[x, 0] = self._ax_util(
                    ax[x, 0],
                    title=self._shift_resid_index(dist_col, res_id_label_shift),
                    xlabel="Time (ns)",
                    ylabel="Distance (Å)",
                )

                # in case .txt format is requested
                # H_bond_df = pd.DataFrame(data={'time': self.distances.index / frame_to_time, 'distance': self.distances[og_column_name]})
                # H_bond_df.to_csv(Path(self.plot_folder, f"H_bond_distance_time_series_{self.sim_name}_{self._shift_resid_index(dist_col.replace(' - ', '__'), res_id_label_shift)}.txt"), sep='\t', index=False)

                ax[x, 1].hist(
                    self.distances[og_column_name],
                    bins=num_bins,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[x, 1] = self._ax_util(
                    ax[x, 1], xlabel="Frequency", ylabel="Distance (Å)"
                )

                last_x_dist = self.distances[
                    self.distances.index
                    > max(self.distances.index) - pmf_last_nth_frames
                ]

                num_substates = self.count_substates(
                    last_x_dist[og_column_name], num_bins=num_bins
                )
                substates_per_edge.append(
                    [dist_col.split(" - ")[0], dist_col.split(" - ")[1], num_substates]
                )
                total_number_of_states += num_substates

                ax[x, 2].hist(
                    last_x_dist[og_column_name],
                    bins=num_bins,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    rasterized=True,
                    orientation="horizontal",
                )
                ax[x, 2] = self._ax_util(
                    ax[x, 2], xlabel="Frequency", ylabel="Distance (Å)"
                )
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {pmf_last_nth_frames/frame_to_time:.0f} ns\nPN = {num_substates}",
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=ax[x, 2].transAxes,
                    fontsize=text_fs,
                )

                bin_centers, PMF = self.calculate_PMF(
                    last_x_dist[og_column_name], num_bins=num_bins
                )

                pmfs = pd.DataFrame(data={"distance": bin_centers, "PMF": PMF})
                pmfs.to_csv(
                    Path(
                        self.plot_folder,
                        f"PMF_{self.sim_name}_{self._shift_resid_index(dist_col.replace(' - ', '__'), res_id_label_shift)}.csv",
                    )
                )

                # in case .txt format is requested
                # pmfs.to_csv(
                #     Path(
                #         self.plot_folder,
                #         f"PMF_{self.sim_name}_{self._shift_resid_index(dist_col.replace(' - ', '__'), res_id_label_shift)}.txt",), sep='\t', index=False
                # )

                ax[x, 3].plot(PMF, bin_centers, color=dist_color, linewidth=2)
                ax[x, 3] = self._ax_util(
                    ax[x, 3], xlabel="PMF (kcal/mol)", ylabel="Distance (Å)"
                )
                total_substates += num_substates

            substates_per_node.append(
                [
                    graph_node,
                    len(distcance_columns),
                    total_substates,
                    round(total_substates / len(distcance_columns), 2),
                ]
            )

            def split_node_name(n):
                return ("-").join(n.split("-")[:-1])

            edges_per_res = [
                f'{ split_node_name(e.split(" - ")[0])} - {split_node_name(e.split(" - ")[1])}'
                for e in distcance_columns
            ]
            unique_res, count = np.unique(edges_per_res, return_counts=True)
            unique_hbond_per_res.append([graph_node, len(unique_res)])

            # print("total_number_of_states", total_number_of_states)
            fig.suptitle(
                self._shift_resid_index(graph_node, res_id_label_shift),
                fontsize=text_fs * 1.2,
            )
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            # fig.tight_layout(h_pad=4.0)
            for img_format in plot_formats:
                fig.savefig(
                    Path(
                        self.plot_folder,
                        f"{self.sim_name}_{self._shift_resid_index(graph_node, res_id_label_shift)}_dist_combined.{img_format}",
                    ),
                    format=img_format,
                )
            plt.close()

        df = (
            pd.DataFrame(
                substates_per_node,
                columns=["amino_acid", "num_Hbonds", "total_substates", "ratio"],
            )
            .drop_duplicates()
            .sort_values(by="amino_acid")
        )
        df.to_csv(
            Path(
                self.plot_folder,
                f"{self.sim_name}_substates_per_amino_acid.csv",
            ),
            index=False,
        )

        with open(
            Path(self.plot_folder, f"{self.sim_name}_PN_per_edges.txt"), "w"
        ) as f:
            for row in substates_per_edge:
                f.write(" ".join(map(str, row)) + "\n")

        pd.DataFrame(
            unique_hbond_per_res,
            columns=["amino_acid", "num_Hbonds_per_res"],
        ).drop_duplicates().sort_values(by="amino_acid").to_csv(
            Path(
                self.plot_folder,
                f"{self.sim_name}_Hbonds_per_amino_acid_residuewise.csv",
            ),
            index=False,
        )


def main():

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--plot_folder",
        required=True,
        help="",
    )
    parser.add_argument(
        "--graphs_info_txt",
        required=True,
        help="",
    )

    parser.add_argument(
        "--pKas_for_frame_csv",
        required=False,
        help="",
    )

    parser.add_argument(
        "--pair_distances_csv",
        required=True,
        help="",
    )

    parser.add_argument(
        "--water_within_csv",
        required=True,
        help="",
    )

    parser.add_argument(
        "--total_water_within_csv",
        required=True,
        help="",
    )

    parser.add_argument(
        "--sim_name",
        help="",
    )
    parser.add_argument(
        "--step",
        default=1,
        help="",
        type=int,
    )
    parser.add_argument(
        "--frame_to_time",
        help="",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--pmf_last_nth_frames",
        help="",
        default=20000,
        type=int,
    )
    parser.add_argument(
        "--plot_formats",
        help="",
        default="['png']",
    )

    parser.add_argument(
        "--res_id_label_shift",
        default=0,
        help="",
        type=int,
    )
    args = parser.parse_args()

    for file in [
        args.graphs_info_txt,
        args.pair_distances_csv,
        args.water_within_csv,
        args.total_water_within_csv,
    ]:
        if not os.path.isfile(file):
            raise FileNotFoundError(
                f"The file {file} does not exist. Please provide all required files to plot the data."
            )

    if not os.path.isdir(args.plot_folder):
        os.makedirs(args.plot_folder)

    base = os.path.basename(args.pair_distances_csv)
    base_name, ext = os.path.splitext(base)

    sim_name = args.sim_name if args.sim_name else base_name.split("_pair_distances")[0]

    plotter = DNetPlot(
        graphs_info_txt=args.graphs_info_txt,
        pKas_for_frame_csv=args.pKas_for_frame_csv,
        pair_distances_csv=args.pair_distances_csv,
        water_within_csv=args.water_within_csv,
        total_water_within_csv=args.total_water_within_csv,
        plot_folder=args.plot_folder,
        sim_name=sim_name,
        step=args.step,
    )

    plotter.create_combined_plots(
        frame_to_time=args.frame_to_time,
        pmf_last_nth_frames=args.pmf_last_nth_frames,
        plot_formats=ast.literal_eval(args.plot_formats),
        res_id_label_shift=args.res_id_label_shift,
    )


if __name__ == "__main__":
    main()
