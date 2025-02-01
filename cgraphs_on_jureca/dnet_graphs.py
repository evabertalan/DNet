import helperfunctions as _hf
import copy
import numpy as np
import MDAnalysis as _mda
import mdhbond as mdh
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import os


class DNetGraphs:
    def __init__(
        self,
        target_folder,
        psf_file,
        dcd_files,
        sim_name=None,
        plot_parameters={},
    ):

        self.plot_parameters = _hf.get_plot_parameters(plot_parameters)
        self.target_folder = target_folder
        self.workfolder = _hf.create_directory(Path(target_folder, "workfolder"))

        self.graph_object_folder = _hf.create_directory(
            Path(self.workfolder, "graph_objects")
        )

        self.helper_files_folder = _hf.create_directory(
            Path(self.workfolder, ".helper_files")
        )

        self.logger = _hf.create_logger(self.helper_files_folder)

        self.psf_file = psf_file
        self.dcd_files = dcd_files

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
        for i, resisdue in enumerate(selected_atoms):
            chain, res_name, res_id = (
                resisdue.segid,
                resisdue.resname,
                resisdue.resid,
            )
            if residuewise:
                res = f"{chain}-{res_name}-{res_id}"
            else:
                atom_name = resisdue.name
                res = f"{chain}-{res_name}-{res_id}-{atom_name}"

            if res in graph.nodes:
                node_positions.update({res: resisdue.position})
        return node_positions

    def calculate_graphs(
        self,
        selection="protein",
        max_water=3,
        exclude_backbone_backbone=True,
        include_backbone_sidechain=False,
        distance=3.5,
        cut_angle=60.0,
        check_angle=True,
        additional_donors=[],
        additional_acceptors=[],
        calcualte_distance=True,
        step=1,
        residuewise=True,
        wrap_dcd=False,
    ):
        self.distance = distance
        self.logger.info(f"H-bond criteria cut off distance: {self.distance} A")

        self.include_backbone_sidechain = include_backbone_sidechain

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

        self.water_graphs_folder = _hf.create_directory(
            Path(self.graph_object_folder, f"{self.max_water}_water_wires")
        )

        if check_angle:
            self.logger.info(f"H-bond criteria cut off angle: {cut_angle} degree")

        self.residuewise = residuewise
        self.logger.info(
            f"Calculating H-bonds {'residuewise' if self.residuewise else 'atomwise'}"
        )

        self.logger.info(
            f"Loading {len(self.dcd_files)} trajectory files for {self.sim_name}"
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
        )
        wba.set_water_wires(water_in_convex_hull=max_water, max_water=max_water)
        wba.compute_average_water_per_wire()
        self.graph_coord_object.update({"wba": wba})

        wba.dump_to_file(
            Path(
                self.water_graphs_folder,
                f"{self.sim_name}_{self.max_water}_water_wires_graph.pickle",
            )
        )

        self.graph = wba.filtered_graph
        self.graph_coord_object.update({"graph": self.graph})

        u = _mda.Universe(self.psf_file, self.dcd_files)
        selected_atoms = u.select_atoms(self.selection)
        self.graph_coord_object.update({"selected_atoms": selected_atoms})

        _hf.pickle_write_file(
            Path(
                self.helper_files_folder,
                f"{self.sim_name}_{self.max_water}_water_nx_graphs.pickle",
            ),
            self.graph,
        )

        _hf.json_write_file(
            Path(
                self.helper_files_folder,
                f"{self.sim_name}_{self.max_water}_water_graph_edge_info.json",
            ),
            _hf.edge_info(wba, self.graph.edges),
        )

        graph_coord_object_loc = Path(
            self.helper_files_folder,
            f"{self.sim_name}_{self.max_water}_water_wires_coord_objects.pickle",
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
        color_edges=False,
        res_id_label_shift=0,
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

        edge_colors_data = {}
        if color_edges:
            edge_value_dict = _hf.read_edge_color_data(self.sim_name, self.psf_file)
            edge_colors_data, cmap, norm = _hf.get_color_map(edge_value_dict)

        for e in graph.edges:
            e0 = _hf.get_node_name(e[0])
            e1 = _hf.get_node_name(e[1])
            if e0 in node_pca_pos.keys() and e1 in node_pca_pos.keys():
                edge_line = [node_pca_pos[e0], node_pca_pos[e1]]
                x = [edge_line[0][0], edge_line[1][0]]
                y = [edge_line[0][1], edge_line[1][1]]

                color = (
                    edge_colors_data[e]
                    if e in edge_colors_data.keys()
                    else self.plot_parameters["graph_color"]
                )
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
                    waters, occ_per_wire, _ = _hf.get_edge_params(wba, graph.edges)
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
            self.logger.info(
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
                    self.logger.info(
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

        for n, values in node_pca_pos.items():
            if n in graph.nodes:
                if n.split("-")[1] in _hf.water_types:
                    ax.scatter(
                        values[0],
                        values[1],
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
                        s=self.plot_parameters["node_size"],
                        zorder=5,
                        edgecolors=self.plot_parameters["graph_color"],
                    )

        if label_nodes:
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

        if color_info:
            cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            cbar.ax.tick_params(labelsize=self.plot_parameters["plot_tick_fontsize"])
            cbar.set_label(
                label=color_bar_label,
                size=self.plot_parameters["plot_label_fontsize"],
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

        plot_folder = _hf.create_directory(
            Path(self.workfolder, f"{self.max_water}_water_wires", self.sim_name)
        )

        waters = f"_max_{self.max_water}_water_bridges" if self.max_water > 0 else ""
        occ = f"_min_occupancy_{occupancy}" if occupancy else ""
        for form in self.plot_parameters["formats"]:
            plt.savefig(
                Path(
                    plot_folder,
                    f"{self.sim_name}{waters}{occ}_graph{is_propka}{is_conservation}{is_backbone}{is_label}.{form}",
                ),
                format=form,
                dpi=self.plot_parameters["plot_resolution"],
            )
        if is_label:
            _hf.write_text_file(
                Path(
                    plot_folder,
                    f"{self.sim_name}{waters}{occ}_water_wire_graph_info.txt",
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

    def get_linear_lenght(self, objects, graph):
        connected_components = _hf.get_connected_components(graph)
        struct_object = (
            objects["structure"]
            if self.type_option == "pdb"
            else objects["selected_atoms"]
        )
        return _hf.calculate_connected_compontents_coordinates(
            connected_components, struct_object, option=self.type_option
        )

    def plot_linear_lenghts(self, occupancy=None, label_nodes=True):
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
            connected_components_coordinates = self.get_linear_lenght(
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

            plot_folder = _hf.create_directory(
                Path(self.workfolder, f"{self.max_water}_water_wires", {self.sim_name})
            )

            waters = (
                f"_max_{self.max_water}_water_bridges" if self.max_water > 0 else ""
            )
            occ = f"_min_occupancy_{occupancy}" if occupancy else ""
            for form in self.plot_parameters["formats"]:
                plt.savefig(
                    f"{plot_folder}{self.sim_name}{waters}{occ}_linear_length{is_backbone}{is_label}.{form}",
                    format=form,
                    dpi=self.plot_parameters["plot_resolution"],
                )
            plt.close()
        else:
            self.logger.warning(
                f"{self.sim_name} has no {self.graph_type} graph. Linear length can not be calculated for this structure."
            )
