import ast
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


class MDProfilePlotter:
    def __init__(
        self,
        cgraphs_info_file,
        pKa_info_file,
        distance_csv,
        water_number_csv,
        total_water_around_res_csv,
        path_to_save_output,
        sim_name,
    ):
        self.path_to_save_output = path_to_save_output
        self.sim_name = sim_name
        graph_nodes = self._get_H_bond_nodes(cgraphs_info_file)
        pKa_nodes = [("-").join(node.split("-")[0:3]) for node in graph_nodes]

        pKas = pd.read_csv(pKa_info_file, index_col="frame")
        self.pKas = pKas.loc[:, self.pKas.columns.isin(pKa_nodes)]
        self.pKas.to_csv(
            Path(self.path_to_save_output, f"pKa_nodes_per_freame_{sim_name}.csv")
        )

        self.distances = pd.read_csv(distance_csv, index_col=0)
        self.water_numbers = pd.read_csv(water_number_csv, index_col=0)

        total_water_around_res = pd.read_csv(total_water_around_res_csv, index_col=0)
        total_water_around_res.columns = (
            total_water_around_res.columns.str.split("-").str[:3].str.join("-")
        )
        self.total_water_around_res = total_water_around_res.loc[
            :, ~total_water_around_res.columns.duplicated()
        ]

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

    def _ax_util(ax, title=None, xlabel=None, ylabel=None, only_integers=False):
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

    def _shift_resid_index(node_names, shift=0):
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
        k_B_kcalmol = 1.987204259e-3
        self.T = T  # Temperature in Kelvin
        kT_kcalmol = k_B_kcalmol * self.T

        # Create histogram (frequency counts) density=True, means it is normalized
        hist, bin_edges = np.histogram(distances, bins=num_bins, density=True)

        # The midpoints of the bins represent the reaction coordinate
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # set minimum probability threshold to avoid rare states on the PMF
        threshold = min_th * np.sum(
            hist
        )  # TODO: maybe define it in a more intuitive way
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

    def create_combined_plots(self):

        pKa_color = "#227351"
        total_water_color = "#8f5038"
        water_color = "#cc7f27"
        dist_color = "#335080"
        text_fs = 32

        frame_to_time = 100
        end_frame_pmf = 20000
        shift = 19

        for i, pKa_column in enumerate(self.pKas.columns):
            distcance_columns = [
                col
                for col in self.distances.columns
                if any(
                    ("-").join(part.split("-")[0:3]) == pKa_column
                    for part in col.split(" - ")
                )
            ]
            water_columns = [
                col
                for col in self.water_numbers.columns
                if any(
                    ("-").join(part.split("-")[0:3]) == pKa_column
                    for part in col.split(" - ")
                )
            ]

            fig, ax = plt.subplots(
                nrows=2 + len(distcance_columns) + len(water_columns),
                ncols=4,
                figsize=(40, (2 + len(distcance_columns) + len(water_columns)) * 8),
                gridspec_kw={"width_ratios": [2, 1, 1, 1]},
            )
            ax[0, 0].plot(
                self.pKas.index / frame_to_time, self.pKas[pKa_column], color=pKa_color
            )
            ax[0, 0] = self._ax_util(
                ax[0, 0],
                title=self._shift_resid_index(pKa_column, shift),
                xlabel="Time (ns)",
                ylabel="pKa",
            )

            ax[0, 1].hist(
                self.pKas[pKa_column],
                bins=100,
                edgecolor=pKa_color,
                color=pKa_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[0, 1] = self._ax_util(ax[0, 1], xlabel="Frequency", ylabel="pKa")

            last_x_pKa = self.pKas[
                self.pKas.index > max(self.pKas.index) - end_frame_pmf
            ]
            ax[0, 2].hist(
                last_x_pKa[pKa_column],
                bins=100,
                edgecolor=pKa_color,
                color=pKa_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[0, 2] = self._ax_util(ax[0, 2], xlabel="Frequency", ylabel="pKa")
            ax[0, 2].text(
                0.95,
                0.95,
                f"last {end_frame_pmf/frame_to_time:.0f} ns",
                horizontalalignment="right",  # Align text to the right
                verticalalignment="top",  # Align text to the top
                transform=ax[0, 2].transAxes,  # Use normalized coordinates (0 to 1)
                fontsize=text_fs,
            )
            ax[0, 3].axis("off")

            ax[1, 0].plot(
                self.total_water_around_res.index / frame_to_time,
                self.total_water_around_res[pKa_column],
                color=total_water_color,
            )
            ax[1, 0] = self._ax_util(
                ax[1, 0],
                title=self._shift_resid_index(pKa_column, shift),
                xlabel="Time (ns)",
                ylabel="#waters",
                only_integers=True,
            )

            ax[1, 1].hist(
                self.total_water_around_res[pKa_column],
                bins=100,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[1, 1] = self._ax_util(
                ax[1, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
            )

            last_x_total_water = self.total_water_around_res[
                self.total_water_around_res.index
                > max(self.total_water_around_res.index) - end_frame_pmf
            ]
            ax[1, 2].hist(
                last_x_total_water[pKa_column],
                bins=100,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[1, 2] = self._ax_util(
                ax[1, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
            )
            ax[1, 2].text(
                0.95,
                0.95,
                f"last {end_frame_pmf/frame_to_time:.0f} ns",
                horizontalalignment="right",  # Align text to the right
                verticalalignment="top",  # Align text to the top
                transform=ax[1, 2].transAxes,  # Use normalized coordinates (0 to 1)
                fontsize=text_fs,
            )
            ax[1, 3].axis("off")

            for k, wat_col in enumerate(water_columns):
                x = 2 + k

                ax[x, 0].plot(
                    self.water_numbers.index / frame_to_time,
                    self.water_numbers[wat_col],
                    color=water_color,
                )
                ax[x, 0] = self._ax_util(
                    ax[x, 0],
                    title=f"# of waters within 3.5 Å of {self._shift_resid_index(wat_col, shift)}",
                    xlabel="Time (ns)",
                    ylabel="#waters",
                    only_integers=True,
                )

                ax[x, 1].hist(
                    self.water_numbers[wat_col],
                    bins=100,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 1] = self._ax_util(
                    ax[x, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
                )

                last_x_water = self.water_numbers[
                    self.water_numbers.index
                    > max(self.water_numbers.index) - end_frame_pmf
                ]
                ax[x, 2].hist(
                    last_x_water[wat_col],
                    bins=100,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 2] = self._ax_util(
                    ax[x, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
                )
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {end_frame_pmf/frame_to_time:.0f} ns",
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=ax[x, 2].transAxes,
                    fontsize=text_fs,
                )
                ax[x, 3].axis("off")

            for j, dist_col in enumerate(distcance_columns):
                x = j + 2 + len(water_columns)
                ax[x, 0].plot(
                    self.distances.index / frame_to_time,
                    self.distances[dist_col],
                    color=dist_color,
                )
                ax[x, 0] = self._ax_util(
                    ax[x, 0],
                    title=self._shift_resid_index(dist_col, shift),
                    xlabel="Time (ns)",
                    ylabel="Distance (Å)",
                )

                ax[x, 1].hist(
                    self.distances[dist_col],
                    bins=100,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 1] = self._ax_util(
                    ax[x, 1], xlabel="Frequency", ylabel="Distance (Å)"
                )

                last_x_dist = self.distances[
                    self.distances.index > max(self.distances.index) - end_frame_pmf
                ]
                ax[x, 2].hist(
                    last_x_dist[dist_col],
                    bins=100,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 2] = self._ax_util(
                    ax[x, 2], xlabel="Frequency", ylabel="Distance (Å)"
                )
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {end_frame_pmf/frame_to_time:.0f} ns",
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=ax[x, 2].transAxes,
                    fontsize=text_fs,
                )

                bin_centers, PMF = self.calculate_PMF(last_x_dist[dist_col])

                pmfs = pd.DataFrame(cols=bin_centers, data=PMF)
                pmfs.to_csv(
                    Path(
                        self.path_to_save_output, f"PMF_{self.sim_name,}_{dist_col}.csv"
                    )
                )

                ax[x, 3].plot(PMF, bin_centers, color=dist_color, linewidth=2)
                ax[x, 3] = self._ax_util(
                    ax[x, 3], xlabel="PMF (kcal/mol)", ylabel="Distance (Å)"
                )

            fig.tight_layout(h_pad=4.0)
            fig.savefig(
                Path(
                    self.path_to_save_output,
                    f"{self.sim_name}_{self._shift_resid_index(pKa_column, shift)}_pKa_dist_combined.png",
                )
            )
