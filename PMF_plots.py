import os
import pandas as pd
import ast
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from pathlib import Path


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

        graph_nodes = self._get_H_bond_nodes(cgraphs_info_file)
        pKa_nodes = [("-").join(node.split("-")[0:3]) for node in graph_nodes]

        pKas = pd.read_csv(pKa_info_file, index_col="frame")
        self.pKas = pKas.loc[:, pKas.columns.isin(pKa_nodes)]
        self.pKas.to_csv(
            Path(path_to_save_output, f"pKa_nodes_per_freame_{sim_name}.csv")
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

    def calculate_PMF(self, distances, T=300, num_bins=100, min_th=0.001):
        #   Constants
        k_B = 1.380649e-23  # Boltzmann constant in J/K
        self.T = T  # Temperature in Kelvin
        kT = k_B * T  # Thermal energy (Joules)

        # Convert kT to units of kcal/mol (1 kT in kcal/mol at 300K is approx. 0.596)
        kT_kcalmol = 0.596

        # Create the histogram (frequency counts) density=True, means it is normalized
        hist, bin_edges = np.histogram(distances, bins=num_bins, density=True)

        # The midpoints of the bins represent the reaction coordinate
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Step 3: Calculate the probability distribution (P(x))
        # Since hist is already normalized by the 'density=True' argument,
        # hist can be treated as P(x)
        threshold = min_th * np.sum(hist)
        P_x = hist[np.where(hist > threshold)[0]]
        bin_centers = bin_centers[np.where(hist > threshold)[0]]

        # Step 4: Calculate the PMF using W(x) = -kT ln P(x)
        # Avoid log of zero by replacing zeros with a small number
        # (or ignoring bins with zero probability)
        P_x[P_x == 0] = 1e-10  # Handling log(0) by a very small number

        # Compute the PMF (in units of kcal/mol)
        PMF = -kT_kcalmol * np.log(P_x)

        # Step 5: Normalize the PMF to set the minimum to 0
        PMF -= np.min(PMF)

        return bin_centers, PMF

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

    def _shift_resid_in_index(node_names, shift=0):
        node_names = node_names.split(" - ")
        updated_node_names = []
        for n in node_names:
            parts = n.split("-")
            parts[2] = str(int(parts[2]) + shift)
            n_joined = ("-").join(parts)
            updated_node_names.append(n_joined)
        return (" - ").join(updated_node_names)

    def create_combined_plots():
        # distances = distances.where(distances <= 20, np.nan)

        pKa_color = "#227351"
        total_water_color = "#8f5038"
        water_color = "#cc7f27"
        dist_color = "#335080"
        test_fs = 32

        frame_to_time = 100
        end_frame_pmf = 20000
        shift = 19

        for i, pKa_column in enumerate(pKas.columns):
            distcance_columns = [
                col
                for col in distances.columns
                if any(
                    ("-").join(part.split("-")[0:3]) == pKa_column
                    for part in col.split(" - ")
                )
            ]
            water_columns = [
                col
                for col in water_numbers.columns
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
            ax[0, 0].plot(pKas.index / frame_to_time, pKas[pKa_column], color=pKa_color)
            ax[0, 0] = ax_util(
                ax[0, 0],
                title=shift_resid_in_index(pKa_column, shift),
                xlabel="Time (ns)",
                ylabel="pKa",
            )

            ax[0, 1].hist(
                pKas[pKa_column],
                bins=100,
                edgecolor=pKa_color,
                color=pKa_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[0, 1] = ax_util(ax[0, 1], xlabel="Frequency", ylabel="pKa")

            last_x_pKa = pKas[pKas.index > max(pKas.index) - end_frame_pmf]
            ax[0, 2].hist(
                last_x_pKa[pKa_column],
                bins=100,
                edgecolor=pKa_color,
                color=pKa_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[0, 2] = ax_util(ax[0, 2], xlabel="Frequency", ylabel="pKa")
            ax[0, 2].text(
                0.95,
                0.95,
                f"last {end_frame_pmf/frame_to_time:.0f} ns",
                horizontalalignment="right",  # Align text to the right
                verticalalignment="top",  # Align text to the top
                transform=ax[0, 2].transAxes,  # Use normalized coordinates (0 to 1)
                fontsize=test_fs,
            )
            ax[0, 3].axis("off")

            ax[1, 0].plot(
                total_water_around_res.index / frame_to_time,
                total_water_around_res[pKa_column],
                color=total_water_color,
            )
            ax[1, 0] = ax_util(
                ax[1, 0],
                title=shift_resid_in_index(pKa_column, shift),
                xlabel="Time (ns)",
                ylabel="#waters",
                only_integers=True,
            )

            ax[1, 1].hist(
                total_water_around_res[pKa_column],
                bins=100,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[1, 1] = ax_util(
                ax[1, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
            )

            last_x_total_water = total_water_around_res[
                total_water_around_res.index
                > max(total_water_around_res.index) - end_frame_pmf
            ]
            ax[1, 2].hist(
                last_x_total_water[pKa_column],
                bins=100,
                edgecolor=total_water_color,
                color=total_water_color,
                alpha=0.4,
                orientation="horizontal",
            )
            ax[1, 2] = ax_util(
                ax[1, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
            )
            ax[1, 2].text(
                0.95,
                0.95,
                f"last {end_frame_pmf/frame_to_time:.0f} ns",
                horizontalalignment="right",  # Align text to the right
                verticalalignment="top",  # Align text to the top
                transform=ax[1, 2].transAxes,  # Use normalized coordinates (0 to 1)
                fontsize=test_fs,
            )
            ax[1, 3].axis("off")

            for k, wat_col in enumerate(water_columns):
                x = 2 + k

                ax[x, 0].plot(
                    water_numbers.index / frame_to_time,
                    water_numbers[wat_col],
                    color=water_color,
                )
                ax[x, 0] = ax_util(
                    ax[x, 0],
                    title=f"# of waters within 3.5 Å of {shift_resid_in_index(wat_col, shift)}",
                    xlabel="Time (ns)",
                    ylabel="#waters",
                    only_integers=True,
                )

                ax[x, 1].hist(
                    water_numbers[wat_col],
                    bins=100,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 1] = ax_util(
                    ax[x, 1], xlabel="Frequency", ylabel="#waters", only_integers=True
                )

                last_x_water = water_numbers[
                    water_numbers.index > max(water_numbers.index) - end_frame_pmf
                ]
                ax[x, 2].hist(
                    last_x_water[wat_col],
                    bins=100,
                    edgecolor=water_color,
                    color=water_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 2] = ax_util(
                    ax[x, 2], xlabel="Frequency", ylabel="#waters", only_integers=True
                )
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {end_frame_pmf/frame_to_time:.0f} ns",
                    horizontalalignment="right",  # Align text to the right
                    verticalalignment="top",  # Align text to the top
                    transform=ax[x, 2].transAxes,  # Use normalized coordinates (0 to 1)
                    fontsize=test_fs,
                )
                ax[x, 3].axis("off")

            for j, dist_col in enumerate(distcance_columns):
                x = j + 2 + len(water_columns)
                ax[x, 0].plot(
                    distances.index / frame_to_time,
                    distances[dist_col],
                    color=dist_color,
                )
                ax[x, 0] = ax_util(
                    ax[x, 0],
                    title=shift_resid_in_index(dist_col, shift),
                    xlabel="Time (ns)",
                    ylabel="Distance (Å)",
                )

                ax[x, 1].hist(
                    distances[dist_col],
                    bins=100,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 1] = ax_util(ax[x, 1], xlabel="Frequency", ylabel="Distance (Å)")

                last_x_dist = distances[
                    distances.index > max(distances.index) - end_frame_pmf
                ]
                ax[x, 2].hist(
                    last_x_dist[dist_col],
                    bins=100,
                    edgecolor=dist_color,
                    color=dist_color,
                    alpha=0.4,
                    orientation="horizontal",
                )
                ax[x, 2] = ax_util(ax[x, 2], xlabel="Frequency", ylabel="Distance (Å)")
                ax[x, 2].text(
                    0.95,
                    0.95,
                    f"last {end_frame_pmf/frame_to_time:.0f} ns",
                    horizontalalignment="right",  # Align text to the right
                    verticalalignment="top",  # Align text to the top
                    transform=ax[x, 2].transAxes,  # Use normalized coordinates (0 to 1)
                    fontsize=test_fs,
                )

                # save the pmf as welll as the data files
                bin_centers, PMF = calc_PMF(last_x_dist[dist_col], ax[j + 1, 3])

                ax[x, 3].plot(PMF, bin_centers, color=dist_color, linewidth=2)
                ax[x, 3] = ax_util(
                    ax[x, 3], xlabel="PMF (kcal/mol)", ylabel="Distance (Å)"
                )

            fig.tight_layout(h_pad=4.0)
            fig.savefig(
                Path(
                    plot_folder,
                    display_name,
                    f"{display_name}_{shift_resid_in_index(pKa_column, shift)}_pKa_dist_combined.png",
                )
            )
