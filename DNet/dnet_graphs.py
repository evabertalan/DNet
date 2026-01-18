import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis.*")


import helperfunctions as _hf
import copy
import numpy as np
import pandas as pd
import MDAnalysis as _mda
import mdhbond as mdh
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import os
import argparse
import glob
import ast


class DNetGraphs:
    def __init__(
        self,
        target_folder,
        psf_file,
        dcd_files,
        sim_name=None,
        plot_parameters={},
        dont_save_graph_objects=False,
    ):

        self.plot_parameters = _hf.get_plot_parameters(plot_parameters)
        self.target_folder = target_folder
        self.workfolder = _hf.create_directory(Path(target_folder, "workfolder"))

        if not dont_save_graph_objects:
            self.graph_object_folder = _hf.create_directory(
                Path(self.workfolder, "graph_objects")
            )

        self.helper_files_folder = _hf.create_directory(
            Path(self.workfolder, ".helper_files")
        )

        self.logger = _hf.create_logger(self.helper_files_folder)

        self.psf_file = psf_file
        self.dcd_files = dcd_files
        self.multi_segments = None

        if not sim_name:
            base = os.path.basename(psf_file)
            self.sim_name, ext = os.path.splitext(base)
        else:
            self.sim_name = sim_name

        self.graph_coord_object = {
            "sim_name": self.sim_name,
            "psf": self.psf_file,
            "dcd": self.dcd_files,
        }
        self.graph_type = "water_wire"

    def _add_node_positions_from_structure(self, selected_atoms, graph, residuewise):
        node_positions = {}
        for i, residue in enumerate(selected_atoms):
            chain, res_name, res_id = (
                residue.segid,
                residue.resname,
                residue.resid,
            )
            if residuewise:
                res = f"{chain}-{res_name}-{res_id}"
            else:
                atom_name = residue.name
                res = f"{chain}-{res_name}-{res_id}-{atom_name}"

            if res in graph.nodes:
                node_positions.update({res: residue.position})
        return node_positions

    def calculate_graphs(
        self,
        selection="protein",
        max_water=3,
        # exclude_backbone_backbone=True,
        include_backbone_sidechain=False,
        distance=3.5,
        cut_angle=60.0,
        check_angle=True,
        additional_donors=[],
        additional_acceptors=[],
        step=1,
        start=None,
        stop=None,
        residuewise=True,
        wrap_dcd=False,
        connected_component_root=None,
        occupancy=None,
        dont_save_graph_objects=False,
        collect_angles=False,
    ):
        self.distance = distance
        self.connected_component_root = connected_component_root
        self.logger.info(f"H-bond criteria cut off distance: {self.distance} A")

        self.include_backbone_sidechain = include_backbone_sidechain

        if include_backbone_sidechain:
            self.logger.info("Including sidechain-backbone interactions")
            additional_donors.append("N")
            additional_acceptors.append("O")

        self.selection = selection
        self.logger.info(f"Atom selection string: {self.selection}")

        self.max_water = max_water
        self.logger.info(
            f"Maximum number of water in water bridges is set to : {self.max_water}"
        )

        if additional_donors or additional_acceptors:
            self.logger.info(
                f"""List of additional donors: {additional_donors}
                List of additional acceptors: {additional_acceptors}
                """
            )

        if check_angle:
            self.logger.info(f"H-bond criteria cut off angle: {cut_angle} degree")

        self.residuewise = residuewise
        self.logger.info(
            f"Calculating H-bonds {'residuewise' if self.residuewise else 'atomwise'}"
        )

        self.logger.info(
            f"""Loading {len(self.dcd_files)} trajectory files for {self.sim_name}.
            From frame {1 if start is None else start} until frame {'last' if stop is None else stop} with a step size of {step}."""
        )
        self.logger.info("This step takes some time...")
        wba = mdh.WireAnalysis(
            self.selection,
            self.psf_file,
            self.dcd_files,
            residuewise=self.residuewise,
            check_angle=check_angle,
            add_donors_without_hydrogen=not check_angle,
            additional_donors=additional_donors,
            additional_acceptors=additional_acceptors,
            distance=distance,
            cut_angle=cut_angle,
            wrap_dcd=wrap_dcd,
            step=step,
            start=start,
            stop=stop,
        )

        angles_per_frame = wba.set_water_wires(
            water_in_convex_hull=max_water, max_water=max_water
        )
        wba.compute_average_water_per_wire()
        if connected_component_root:
            if occupancy:
                wba.filter_occupancy(occupancy)
            wba.filter_connected_component(connected_component_root)

        self.graph_coord_object.update({"wba": wba})

        self.graph = wba.filtered_graph
        self.graph_coord_object.update({"graph": self.graph})

        u = _mda.Universe(self.psf_file, self.dcd_files)

        selected_atoms = u.select_atoms(self.selection)

        if len([seg.segid for seg in u.select_atoms("protein").segments]) > 1:
            self.plot_parameters["show_chain_label"] = True
            self.multi_segments = [
                seg.segid for seg in u.select_atoms("protein").segments
            ]
        self.graph_coord_object.update({"selected_atoms": selected_atoms})

        if self.connected_component_root:
            plot_folder = _hf.create_directory(
                Path(
                    self.workfolder,
                    f"{self.max_water}_water_wires_connected_components",
                    self.connected_component_root,
                    self.sim_name,
                )
            )
        else:
            plot_folder = _hf.create_directory(
                Path(self.workfolder, f"{self.max_water}_water_wires", self.sim_name)
            )

        if collect_angles:
            angles_per_frame.to_csv(Path(plot_folder, f"{self.sim_name}_angles.csv")),

        df = pd.DataFrame.from_dict(
            _hf.edge_info(wba, self.graph.edges), orient="index"
        ).reset_index()
        df.columns = ["edge", "water", "occupancy"]
        df["edge"] = df["edge"].str.replace(":", "_")

        waters = f"_max_{self.max_water}_water_bridges"
        df.to_csv(
            Path(
                plot_folder,
                f"{self.sim_name}{waters}_water_occupancy_all_edge_info.txt",
            ),
            sep="\t",
            index=False,
        )

        if occupancy:
            df[df["occupancy"] >= occupancy].to_csv(
                Path(
                    plot_folder,
                    f"{self.sim_name}{waters}_water_{occupancy}_occupancy_all_edge_info.txt",
                ),
                sep="\t",
                index=False,
            )

        if not dont_save_graph_objects:
            if connected_component_root:
                root = f"_{self.connected_component_root}_"
                self.water_graphs_folder = _hf.create_directory(
                    Path(
                        self.graph_object_folder,
                        f"{self.max_water}_water_wires_connected_components",
                    )
                )
            else:
                root = ""
                self.water_graphs_folder = _hf.create_directory(
                    Path(self.graph_object_folder, f"{self.max_water}_water_wires")
                )
            wba.dump_to_file(
                Path(
                    self.water_graphs_folder,
                    f"{self.sim_name}{root}{self.max_water}_water_wires_graph.pickle",
                )
            )
            _hf.pickle_write_file(
                Path(
                    self.helper_files_folder,
                    f"{self.sim_name}{root}{self.max_water}_water_nx_graphs.pickle",
                ),
                self.graph,
            )

            _hf.json_write_file(
                Path(
                    self.helper_files_folder,
                    f"{self.sim_name}{root}{self.max_water}_water_graph_edge_info.json",
                ),
                _hf.edge_info(wba, self.graph.edges),
            )

            graph_coord_object_loc = Path(
                self.helper_files_folder,
                f"{self.sim_name}{root}{self.max_water}_water_wires_coord_objects.pickle",
            )
            _hf.pickle_write_file(
                graph_coord_object_loc,
                self.graph_coord_object,
            )
            self.logger.info(f"Graph object is saved as: {graph_coord_object_loc}")

        self.node_positions = self._add_node_positions_from_structure(
            selected_atoms, self.graph, self.residuewise
        )
        self.pca_positions = _hf.calculate_pca_positions(self.node_positions)

    def _get_node_positions(self, pca=True):
        node_pos = {}
        for node in self.graph.nodes:
            n = _hf.get_node_name(node)
            if (
                n not in self.node_positions.keys()
                or n.split("-")[1] in _hf.water_types
            ):
                chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                coords = (
                    self.graph_coord_object["selected_atoms"]
                    .select_atoms("resid " + res_id)
                    .positions[0]
                )
                if coords is not None:
                    node_pos.update({n: list(coords)})
            else:
                node_pos.update({n: self.node_positions[n]})
        if pca:
            return _hf.calculate_pca_positions(node_pos)
        else:
            return node_pos

    def plot_graphs(
        self,
        label_nodes=True,
        label_edges=True,
        xlabel="PCA projected xy plane",
        ylabel="Z coordinates (Å)",
        occupancy=None,
        color_propka=False,
        color_data=False,
        node_color_selection=None,
        node_color_map="viridis",
        color_edges_by=None,
        res_id_label_shift=0,
        color_edge_by_occupnacy=False,
    ):
        wba = self.graph_coord_object["wba"]
        if occupancy:
            wba.filter_occupancy(occupancy)
            graph = wba.filtered_graph
        else:
            graph = self.graph_coord_object["graph"]

        self.logger.debug(f"Creating water wire graph for {self.sim_name}")
        fig, ax = _hf.create_plot(
            title=f"""Water wire graph of structure {self.sim_name}
            Selection:{self.selection[1:-16]}""",
            xlabel=xlabel,
            ylabel=ylabel,
            plot_parameters=self.plot_parameters,
        )
        node_pca_pos = self._get_node_positions()
        node_pca_pos = _hf.check_projection_sign(node_pca_pos, self.pca_positions)

        edge_value_dict = {}
        if color_edges_by:
            edge_value_dict = _hf.read_edge_color_data(color_edges_by)
            edge_value_dict = {
                k: v for k, v in edge_value_dict.items() if k in graph.edges
            }
            cmap, norm, edge_colors = _hf.get_edge_color_map(
                list(edge_value_dict.values())
            )
            color_bar_label = "Edge value"

        elif color_edge_by_occupnacy:
            values = np.linspace(0.1, 1, 10)
            cmap, norm, occupany_colors = _hf.get_edge_color_map(values)
            color_bar_label = "H-bond occupancy"

        waters, occ_per_wire, _ = _hf.get_edge_params(wba, graph.edges)
        for e in graph.edges:
            e0 = _hf.get_node_name(e[0])
            e1 = _hf.get_node_name(e[1])
            if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                x = [edge_line[0][0], edge_line[1][0]]
                y = [edge_line[0][1], edge_line[1][1]]

                if e in edge_value_dict.keys():
                    val = edge_value_dict[e]
                    color = edge_colors.to_rgba(val)

                elif color_edge_by_occupnacy:
                    occ = occ_per_wire[list(graph.edges).index(e)]
                    color = occupany_colors.to_rgba(occ)

                else:
                    color = self.plot_parameters["graph_color"]

                ax.plot(
                    x,
                    y,
                    color=color,
                    marker="o",
                    linewidth=self.plot_parameters["edge_width"],
                    markersize=self.plot_parameters["node_size"] * 0.01,
                    markerfacecolor=self.plot_parameters["graph_color"],
                    markeredgecolor=self.plot_parameters["graph_color"],
                )

                if label_edges:
                    ax.annotate(
                        np.round(waters[list(graph.edges).index(e)], 1),
                        (x[0] + (x[1] - x[0]) / 2, y[0] + (y[1] - y[0]) / 2),
                        color="indianred",
                        fontsize=self.plot_parameters["edge_label_size"],
                        weight="bold",
                    )
                    if occupancy:
                        ax.annotate(
                            int(occ_per_wire[list(graph.edges).index(e)] * 100),
                            (
                                x[0] + (x[1] - x[0]) / 2,
                                y[0] + (y[1] - 1.0 - y[0]) / 2,
                            ),
                            color="green",
                            fontsize=self.plot_parameters["edge_label_size"],
                        )

        color_info = {}
        lab = " with labels" if label_nodes else ""
        if color_propka and color_data:
            self.logger.warning(
                f"Can not color plot by propka and external data values at the same time. Please select just one coloring option!"
            )
        else:
            struct_object = self.graph_coord_object["selected_atoms"]
            selected_nodes = struct_object.select_atoms(str(node_color_selection))
            if color_propka:
                propka_file = Path(self.target_folder, f"{self.sim_name}.propka")
                if propka_file.exists():
                    color_info = _hf.read_propka_file(propka_file, selected_nodes)
                else:
                    self.logger.warning(
                        f"{self.sim_name}.propka not found. To color residues by pKa values, place the propka file in the PDB folder, next to the PDB file."
                    )
                if len(color_info):
                    value_colors, cmap, norm = _hf.get_color_map(
                        color_info, color_map=node_color_map
                    )
                    self.logger.info(f"Color {self.sim_name} by pKa values{lab}.")
                    color_bar_label = "pKa value"
                else:
                    self.logger.info(
                        f"{self.sim_name}.propka does not contain the selected residues. Please update the Residues to color!"
                    )

            elif color_data:
                color_info = _hf.read_color_data_file(
                    self.sim_name, self.target_folder, selected_nodes
                )
                if color_info is None:
                    self.logger.error(
                        f"No {self.sim_name}_data.txt file was found in {self.target_folder}."
                    )
                elif len(color_info):
                    value_colors, cmap, norm = _hf.get_color_map(
                        color_info, color_map=node_color_map
                    )
                    color_bar_label = "Amino acid data value"
                    self.logger.info(
                        f"Color {self.sim_name} by values from external data file{lab}."
                    )
                else:
                    self.logger.error(
                        f"The content of {self.sim_name}_data.txt is invalid or no residues found in the file matching the selection for node coloring"
                    )

        markers = ["o", ",", "v", "p", "D", "*", "h", "H", "X"]
        for n, values in node_pca_pos.items():
            if n in graph.nodes:
                if self.multi_segments:
                    marker_shape = markers[self.multi_segments.index(n.split("-")[0])]
                else:
                    marker_shape = "o"

                if n.split("-")[1] in _hf.water_types:
                    ax.scatter(
                        values[0],
                        values[1],
                        marker=marker_shape,
                        color=self.plot_parameters["water_node_color"],
                        s=self.plot_parameters["node_size"] * 0.7,
                        zorder=5,
                    )
                elif n.split("-")[1] in _hf.amino_d.keys():

                    color = (
                        value_colors[n]
                        if n in color_info.keys()
                        else self.plot_parameters["graph_color"]
                    )

                    ax.scatter(
                        values[0],
                        values[1],
                        color=color,
                        marker=marker_shape,
                        s=self.plot_parameters["node_size"],
                        zorder=5,
                        edgecolors=self.plot_parameters["graph_color"],
                    )
                else:
                    color = (
                        value_colors[n]
                        if n in color_info.keys()
                        else self.plot_parameters["non_prot_color"]
                    )
                    ax.scatter(
                        values[0],
                        values[1],
                        color=color,
                        marker=marker_shape,
                        s=self.plot_parameters["node_size"],
                        zorder=5,
                        edgecolors=self.plot_parameters["graph_color"],
                    )

        if label_nodes:
            self.logger.info(f"Shifting resid labels with {res_id_label_shift}")
            for n in graph.nodes:
                n = _hf.get_node_name(n)
                if n in node_pca_pos.keys():
                    values = node_pca_pos[n]
                    chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                    if res_name in _hf.water_types:
                        pass  # temporary turn off water labels
                        # ax.annotate(
                        #     f"W{res_id}",
                        #     (values[0] + 0.2, values[1] - 0.25),
                        #     fontsize=self.plot_parameters["node_label_size"],
                        # )
                    elif res_name in _hf.amino_d.keys():
                        res_label = (
                            f"{chain_id}-{_hf.amino_d[res_name]}{res_id}"
                            if self.plot_parameters["show_chain_label"]
                            else f"{_hf.amino_d[res_name]}{int(res_id)+res_id_label_shift}"
                        )

                        ax.annotate(
                            res_label,
                            (values[0] + 0.2, values[1] - 0.26),
                            fontsize=self.plot_parameters["node_label_size"],
                        )
                    else:
                        res_label = (
                            f"{chain_id}-{res_name}{res_id}"
                            if self.plot_parameters["show_chain_label"]
                            else f"{res_name}{int(res_id)+res_id_label_shift}"
                        )
                        ax.annotate(
                            res_label,
                            (values[0] + 0.2, values[1] - 0.25),
                            fontsize=self.plot_parameters["node_label_size"],
                            color=self.plot_parameters["non_prot_color"],
                        )

        if color_info or color_edge_by_occupnacy or color_edges_by:
            cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            if all(isinstance(x, np.integer) for x in cbar.get_ticks()):
                tick_centers = (norm.boundaries[:-1] + norm.boundaries[1:]) / 2
                cbar.ax.set_yticks(tick_centers)
                cbar.ax.set_yticklabels(norm.boundaries[:-1])

            cbar.ax.tick_params(
                labelsize=self.plot_parameters["plot_tick_fontsize"] - 3
            )
            cbar.set_label(
                label=color_bar_label,
                size=self.plot_parameters["plot_label_fontsize"] - 3,
            )

        plt.tight_layout()
        is_label = "_labeled" if label_nodes else ""
        is_propka = "_pKa_color" if color_info and color_propka else ""
        is_conservation = "_data_color" if color_info and color_data else ""
        is_backbone = (
            "_backbone"
            if hasattr(self, "include_backbone_sidechain")
            and self.include_backbone_sidechain
            else ""
        )

        if self.connected_component_root:
            plot_folder = _hf.create_directory(
                Path(
                    self.workfolder,
                    f"{self.max_water}_water_wires_connected_components",
                    self.connected_component_root,
                    self.sim_name,
                )
            )
        else:
            plot_folder = _hf.create_directory(
                Path(self.workfolder, f"{self.max_water}_water_wires", self.sim_name)
            )

        waters = f"_max_{self.max_water}_water_bridges"
        occ = f"_min_occupancy_{occupancy}" if occupancy else ""
        root = (
            f"_{self.connected_component_root}" if self.connected_component_root else ""
        )
        for form in self.plot_parameters["formats"]:
            plt.savefig(
                Path(
                    plot_folder,
                    f"{self.sim_name}{root}{waters}{occ}_graph{is_propka}{is_conservation}{is_backbone}{is_label}.{form}",
                ),
                format=form,
                dpi=self.plot_parameters["plot_resolution"],
            )
        if is_label:
            _hf.write_text_file(
                Path(
                    plot_folder,
                    f"{self.sim_name}{root}{waters}{occ}_water_wire_graph_info.txt",
                ),
                [
                    "Water wire graph of " + self.sim_name,
                    "\nSelection string: " + str(self.selection[0:-15]),
                    "\nNumber of maximum water molecules allowed in the bridge: "
                    + str(self.max_water),
                    (
                        "\nMinimum H-bond occupancy: " + str(occupancy)
                        if occupancy
                        else ""
                    ),
                    "\n",
                    "\nNumber of nodes in "
                    + self.sim_name
                    + ": "
                    + str(len(graph.nodes)),
                    "\nNumber of edges in "
                    + self.sim_name
                    + ": "
                    + str(len(graph.edges)),
                    "\n",
                    "\nList of nodes: " + str(graph.nodes),
                    "\n",
                    "\nList of edges: " + str(graph.edges),
                ],
            )
        plt.close()

    def get_linear_length(self, objects, graph):
        connected_components = _hf.get_connected_components(graph)
        struct_object = (
            objects["structure"]
            if self.type_option == "pdb"
            else objects["selected_atoms"]
        )
        return _hf.calculate_connected_compontents_coordinates(
            connected_components, struct_object, option=self.type_option
        )

    def plot_linear_lengths(self, occupancy=None, label_nodes=True):
        self.logger.info("Plotting linear lengths for continuous network components")
        self.logger.debug("Creating linear length plot for " + self.sim_name)
        if "graph" in self.graph_coord_object.keys():

            wba = self.graph_coord_object["wba"]
            if occupancy:
                wba.filter_occupancy(occupancy)
                graph = wba.filtered_graph
            else:
                graph = self.graph_coord_object["graph"]

            self.logger.debug(
                "Creating "
                + self.graph_type
                + " linear length plot for: "
                + self.sim_name
            )
            connected_components_coordinates = self.get_linear_length(
                self.graph_coord_object, graph
            )
            plot_parameters = copy.deepcopy(self.plot_parameters)
            plot_parameters["figsize"] = (
                min(1 + int(len(connected_components_coordinates)), 100),
                plot_parameters["figsize"][1],
            )

            fig, ax = _hf.create_plot(
                title=f"Linear length of continuous water wire subnetworks \nalong the Z-axis in structure {self.sim_name}\nSelection: {self.selection[1:-16]}",
                xlabel="# of nodes",
                ylabel="Z-axis coordinates (Å)",
                plot_parameters=plot_parameters,
            )

            for i, g in enumerate(connected_components_coordinates):
                for j in range(len(g)):
                    node_name = connected_components_coordinates[i][j][0]
                    chain_id, res_name, res_id = (
                        node_name.split("-")[0],
                        node_name.split("-")[1],
                        node_name.split("-")[2],
                    )
                    if res_name in _hf.water_types:
                        color = self.plot_parameters["water_node_color"]
                    elif res_name in _hf.amino_d.keys():
                        color = self.plot_parameters["graph_color"]
                    else:
                        color = self.plot_parameters["non_prot_color"]
                    z_coords = connected_components_coordinates[i][j][1][2]
                    ax.scatter(
                        i,
                        z_coords,
                        color=color,
                        s=self.plot_parameters["node_size"] * 0.8,
                    )
                    if label_nodes:
                        if res_name in _hf.water_types:
                            pass
                        elif res_name in _hf.amino_d.keys():
                            res_label = (
                                f"{chain_id}-{_hf.amino_d[res_name]}{res_id}"
                                if self.plot_parameters["show_chain_label"]
                                else f"{_hf.amino_d[res_name]}{res_id}"
                            )
                            ax.annotate(
                                res_label,
                                (i, z_coords),
                                fontsize=self.plot_parameters["node_label_size"],
                                zorder=6,
                            )
                        else:
                            res_label = (
                                f"{chain_id}-{res_name}{res_id}"
                                if self.plot_parameters["show_chain_label"]
                                else f"{res_name}{res_id}"
                            )
                            ax.annotate(
                                res_label,
                                (i, z_coords),
                                fontsize=self.plot_parameters["node_label_size"],
                                zorder=6,
                            )

            ax.set_xticks(np.arange(len(connected_components_coordinates)))
            ax.set_xticklabels([len(c) for c in connected_components_coordinates])

            plt.tight_layout()
            is_label = "_labeled" if label_nodes else ""
            is_backbone = (
                "_backbone"
                if hasattr(self, "include_backbone_sidechain")
                and self.include_backbone_sidechain
                else ""
            )

            if self.connected_component_root:
                plot_folder = _hf.create_directory(
                    Path(
                        self.workfolder,
                        f"{self.max_water}_water_wires_connected_components",
                        self.connected_component_root,
                        self.sim_name,
                    )
                )
            else:
                plot_folder = _hf.create_directory(
                    Path(
                        self.workfolder, f"{self.max_water}_water_wires", self.sim_name
                    )
                )

            waters = (
                f"_max_{self.max_water}_water_bridges" if self.max_water > 0 else ""
            )
            occ = f"_min_occupancy_{occupancy}" if occupancy else ""
            root = (
                f"_{self.connected_component_root}"
                if self.connected_component_root
                else ""
            )
            for form in self.plot_parameters["formats"]:
                plt.savefig(
                    f"{plot_folder}{self.sim_name}{root}{waters}{occ}_linear_length{is_backbone}{is_label}.{form}",
                    format=form,
                    dpi=self.plot_parameters["plot_resolution"],
                )
            plt.close()

        else:
            self.logger.warning(
                f"{self.sim_name} has no {self.graph_type} graph. Linear length can not be calculated for this structure."
            )


def main():
    parser = argparse.ArgumentParser(
        description="Analyze molecular dynamics trajectories and generate water wire graphs."
    )
    parser.add_argument(
        "psf",
        help="Path to the PSF (Protein Structure File) used for molecular dynamics simulations.",
    )
    parser.add_argument(
        "dcd",
        nargs="+",
        help="Path(s) to the DCD (trajectory) file(s). Supports wildcard patterns (e.g., '*.dcd').",
    )
    parser.add_argument(
        "--output_folder",
        type=str,
        help="Directory where output files and plots will be saved. If not provided, defaults to the location of the PSF file.",
    )

    parser.add_argument(
        "--max_water",
        type=int,
        default=3,
        help="Maximum number of water molecules allowed in water wire connections (default: 3).",
    )
    parser.add_argument(
        "--occupancy",
        type=float,
        default=0.1,
        help="Minimum hydrogen bond occupancy required to include an edge in the graph (default: 0.1, which means 10% occupancy).",
    )
    parser.add_argument(
        "--distance",
        type=float,
        default=3.5,
        help="The distance criterion for the  H-bond search, measured between the heavy atoms. The default value is 3.5Å.",
    )
    parser.add_argument(
        "--cut_angle",
        type=float,
        default=60,
        help="Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. The default value is 60°.",
    )
    parser.add_argument(
        "--selection",
        type=str,
        default="protein",
        help="Atom selection string for defining the region of interest in the molecular system for graph calculation (default: 'protein').",
    )

    parser.add_argument(
        "--additional_donors",
        type=str,
        default="[]",
        help="""List of additional hydrogen bond donor atoms, formatted as a Python list (e.g., "['N', 'S']").""",
    )
    parser.add_argument(
        "--additional_acceptors",
        type=str,
        default="[]",
        help="""List of additional hydrogen bond acceptor atoms, formatted as a Python list (e.g., "['O', 'F']").""",
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
        default=1,
        help="Step size for iterating through the trajectory frames. For example, '--step 10' "
        "processes every 10th frame to reduce computation time. (Default is 1).",
    )

    parser.add_argument(
        "--residuewise",
        default=True,
        action="store_true",
        help="Calculate hydrogen bonds at the residue level instead of the atomic level (default: True).",
    )

    parser.add_argument(
        "--atomwise",
        dest="residuewise",
        action="store_false",
        help="Calculate hydrogen bonds at the atomic level instead of the residue level (overrides --residuewise).",
    )

    parser.add_argument(
        "--wrap_dcd",
        type=str,
        choices=["true", "false"],
        default="true",
        help="Apply periodic boundary condition wrapping to keep molecules inside the simulation box. Use 'true' or 'false' (default: true).",
    )

    parser.add_argument(
        "--root",
        type=str,
        help="Root node for connected component search. In a form of SEGNAME-RESANEM-RESID, e.g: A-ASP-213",
    )
    parser.add_argument(
        "--res_id_label_shift",
        default=0,
        type=int,
        help="Shift residue ID labels by a specified amount in plots (default: 0).",
    )

    parser.add_argument(
        "--color_data",
        action="store_true",
        help="Color nodes in the graph based on external data values.",
    )

    parser.add_argument(
        "--node_color_selection",
        default="protein",
        help="Selection criteria for which nodes to color in the graph (default: 'protein').",
    )

    parser.add_argument(
        "--node_color_map",
        default="coolwarm_r",
        help="Colormap used for node coloring (default: 'coolwarm_r').",
    )

    parser.add_argument(
        "--plot_parameters",
        default="{}",
        help="""Dictionary of plot parameters formatted as a string : {'graph_color': '#666666', 'formats': ['png', 'eps']}" """,
    )

    parser.add_argument(
        "--include_backbone",
        action="store_true",
        help="Include interactions between backbone and sidechain atoms in the analysis.",
    )

    parser.add_argument(
        "--no_label_plots",
        action="store_true",
        help="Creates all the graph plots without the nodes and edges labels as well.",
    )

    parser.add_argument(
        "--dont_save_graph_objects",
        action="store_true",
        help="Don't save the metadata and full graph objects of the calculations. Use this flag if there is not enough space for the calculation results or when the graph objects are not needed for further calculations or analysis.",
    )

    parser.add_argument(
        "--color_edge_by_occupnacy",
        action="store_true",
        help="Color graph edges according to a gray color scale matching the occupancy values.",
    )

    parser.add_argument(
        "--color_edges_by_file",
        type=str,
        help="Path to the .txt file, which contains values according to edges need to be colored. Each row of the .txt file needs to be: edge1 edge2 value ",
    )

    parser.add_argument(
        "--collect_angles",
        default=False,
        action="store_true",
        help="Create a csv file with the angles of all donor-acceptor pairs that are within the set H-bond distance criterion in each frame (default: False).",
    )
    args = parser.parse_args()

    base = os.path.basename(args.psf)
    base_name, ext = os.path.splitext(base)

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

    if args.color_edge_by_occupnacy and args.color_edges_by_file:
        raise ValueError(
            "color_edge_by_occupnacy and color_edges_by_file can't be used at the same time. Edges can be colored only by either of the values."
        )

    dnet_graphs = DNetGraphs(
        target_folder=output_folder,
        psf_file=args.psf,
        dcd_files=dcd_files,
        sim_name=base_name,
        plot_parameters=ast.literal_eval(args.plot_parameters),
        dont_save_graph_objects=args.dont_save_graph_objects,
    )
    dnet_graphs.calculate_graphs(
        max_water=int(args.max_water),
        check_angle=True,
        selection=args.selection,
        additional_donors=ast.literal_eval(args.additional_donors),
        additional_acceptors=ast.literal_eval(args.additional_acceptors),
        residuewise=args.residuewise,
        distance=args.distance,
        cut_angle=args.cut_angle,
        wrap_dcd=wrap_dcd,
        step=args.step,
        start=args.start,
        stop=args.stop,
        include_backbone_sidechain=args.include_backbone,
        occupancy=float(args.occupancy),
        connected_component_root=args.root,
        dont_save_graph_objects=args.dont_save_graph_objects,
        collect_angles=args.collect_angles,
    )

    dnet_graphs.plot_graphs(
        label_nodes=True,
        xlabel="PCA projected membrane plane (Å)",
        ylabel="Membrane normal (Å)",
        occupancy=float(args.occupancy),
        color_data=args.color_data,
        node_color_selection=args.node_color_selection,
        node_color_map=args.node_color_map,
        res_id_label_shift=int(args.res_id_label_shift),
        color_edge_by_occupnacy=args.color_edge_by_occupnacy,
        color_edges_by=args.color_edges_by_file,
    )

    if args.no_label_plots:
        dnet_graphs.plot_graphs(
            label_nodes=False,
            label_edges=False,
            xlabel="PCA projected membrane plane (Å)",
            ylabel="Membrane normal (Å)",
            occupancy=float(args.occupancy),
            color_data=args.color_data,
            node_color_selection=args.node_color_selection,
            node_color_map=args.node_color_map,
            res_id_label_shift=int(args.res_id_label_shift),
            color_edge_by_occupnacy=args.color_edge_by_occupnacy,
            color_edges_by=args.color_edges_by_file,
        )


if __name__ == "__main__":
    main()
