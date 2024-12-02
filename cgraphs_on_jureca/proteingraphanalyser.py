from . import helperfunctions as _hf
import copy
import networkx as nx
import numpy as np
import MDAnalysis as _mda
from MDAnalysis.analysis import distances
from . import mdhbond as mdh
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import os


class ProteinGraphAnalyser:
    def __init__(
        self,
        target_folder,
        psf_file,
        dcd_files,
        sim_name=None,
        plot_parameters={},
    ):

        self.plot_parameters = _hf.get_plot_parameters(plot_parameters)
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

    def add_reference_from_structure(self, s, g):
        residuewise = len(list(g.nodes)[0].split("-")) == 3
        for i, resisdue in enumerate(s):
            if residuewise:
                chain, res_name, res_id = (
                    resisdue.segid,
                    resisdue.resname,
                    resisdue.resid,
                )
                res = f"{chain}-{res_name}-{res_id}"
            else:
                chain, res_name, res_id, atom_name = (
                    resisdue.segid,
                    resisdue.resname,
                    resisdue.resid,
                    resisdue.name,
                )
                res = f"{chain}-{res_name}-{res_id}-{atom_name}"
            if res not in self.reference_coordinates.keys() and res in g.nodes:
                self.reference_coordinates.update({res: resisdue.position})

    def calculate_graphs(
        self,
        graph_type="water_wire",
        selection="protein",
        max_water=3,
        exclude_backbone_backbone=True,
        include_backbone_sidechain=False,
        include_waters=True,
        distance=3.5,
        cut_angle=60.0,
        check_angle=True,
        additional_donors=[],
        additional_acceptors=[],
        calcualte_distance=True,
        step=1,
        residuewise=True,
    ):
        self.distance = distance
        self.cut_angle = cut_angle
        self.include_waters = include_waters
        self.include_backbone_sidechain = include_backbone_sidechain
        assert type(additional_donors) is list and type(additional_acceptors) is list
        if additional_donors or additional_acceptors:
            self.logger.info(
                f"List of additional donors: {additional_donors}\nList of additional acceptors: {additional_acceptors}"
            )
        if selection != "protein":
            self.logger.info(f"Atom selection string: {selection}")
        self.graph_type = graph_type
        self.logger.info(f"Calculating graphs for {self.graph_type} analysis.")
        if self.type_option == "pdb":
            self.selection = f"({selection}) or resname BWX"
            # self.get_reference_coordinates(self.reference_pdb)
            try:
                self.graph_type in ["water_wire", "hbond"]
            except ValueError:
                raise ValueError('Given graph_type has to be "water_wire" or "hbond"')
            self.logger.info(f"H-bond criteria cut off distance: {distance} A")
            if check_angle:
                self.logger.info(f"H-bond criteria cut off angle: {cut_angle} degree")
            if self.graph_type == "water_wire":
                self.water_graphs_folder = _hf.create_directory(
                    f"{self.graph_object_folder}{self.max_water}_water_wires/"
                )
                self.logger.info(
                    f"Maximum number of water in water bridges is set to: {max_water}"
                )
                self.max_water = max_water
                for name, objects in self.graph_coord_objects.items():
                    pdb_file, pdb_code = objects["file"], name
                    if len(_hf.water_in_pdb(pdb_file)) == 0:
                        self.logger.warning(
                            f"There are no water molecules in {pdb_code}. Water wire can not be calculated."
                        )
                    else:
                        self.logger.info(
                            f"Calculating {self.graph_type} graph for: {pdb_code}"
                        )
                        if len(_hf.water_in_pdb(pdb_file)) == 0:
                            self.logger.warning(
                                f"There are no water molecules in {pdb_code}. Water wire can not be calculated. Please use the H-bond network option."
                            )
                        else:
                            try:
                                wba = mdh.WireAnalysis(
                                    self.selection,
                                    pdb_file,
                                    residuewise=True,
                                    check_angle=check_angle,
                                    add_donors_without_hydrogen=not check_angle,
                                    additional_donors=additional_donors,
                                    additional_acceptors=additional_acceptors,
                                    distance=distance,
                                    cut_angle=cut_angle,
                                )
                                wba.set_water_wires(max_water=max_water)
                                wba.compute_average_water_per_wire()
                                g = wba.filtered_graph
                                nx.write_gpickle(
                                    g,
                                    self.water_graphs_folder
                                    + pdb_code
                                    + "_"
                                    + self.graph_type
                                    + "_graphs.pickle",
                                )
                                self.graph_coord_objects[pdb_code].update({"graph": g})
                                tmp = copy.copy(wba)
                                tmp._universe = None
                                self.graph_coord_objects[pdb_code].update({"wba": tmp})
                                edge_info = _hf.edge_info(wba, g.edges)
                                _hf.json_write_file(
                                    self.helper_files_folder
                                    + pdb_code
                                    + "_"
                                    + self.graph_type
                                    + "_graph_edge_info.json",
                                    edge_info,
                                )
                                s = objects["structure"].select_atoms(self.selection)
                                self.add_reference_from_structure(s, g)
                            except:
                                self.logger.warning(
                                    f"Graph can't be calculated for {pdb_file}"
                                )

            elif self.graph_type == "hbond":
                if include_backbone_sidechain:
                    self.logger.info("Including sidechain-backbone interactions")
                    additional_donors.append("N")
                    additional_acceptors.append("O")
                if not exclude_backbone_backbone:
                    self.logger.info("Including backbone-backbone interactions")
                for name, objects in self.graph_coord_objects.items():
                    pdb_file, pdb_code = objects["file"], name
                    self.logger.info(
                        f"Calculating {self.graph_type} graph for: {pdb_code}"
                    )
                    try:
                        hba = mdh.HbondAnalysis(
                            self.selection,
                            pdb_file,
                            residuewise=True,
                            check_angle=check_angle,
                            add_donors_without_hydrogen=not check_angle,
                            additional_donors=additional_donors,
                            additional_acceptors=additional_acceptors,
                            distance=distance,
                            cut_angle=cut_angle,
                        )
                        hba.set_hbonds_in_selection(
                            exclude_backbone_backbone=exclude_backbone_backbone
                        )
                        if len(_hf.water_in_pdb(pdb_file)) > 0 and include_waters:
                            hba.set_hbonds_in_selection_and_water_around(distance)
                        g = hba.filtered_graph
                        nx.write_gpickle(
                            g,
                            self.graph_object_folder
                            + pdb_code
                            + "_"
                            + self.graph_type
                            + "_graphs.pickle",
                        )
                        self.graph_coord_objects[pdb_code].update({"graph": g})
                        s = objects["structure"].select_atoms(
                            f"resname BWX or {self.selection} or {_hf.water_def}"
                        )
                        self.add_reference_from_structure(s, g)
                    except:
                        self.logger.warning(f"Graph can't be calculated for {pdb_file}")

                    if calcualte_distance:
                        try:
                            dist_hba = mdh.HbondAnalysis(
                                self.selection,
                                pdb_file,
                                residuewise=False,
                                check_angle=check_angle,
                                add_donors_without_hydrogen=not check_angle,
                                additional_donors=additional_donors,
                                additional_acceptors=additional_acceptors,
                                distance=distance,
                                cut_angle=cut_angle,
                            )
                            dist_hba.set_hbonds_in_selection(
                                exclude_backbone_backbone=exclude_backbone_backbone
                            )
                            if len(_hf.water_in_pdb(pdb_file)) > 0 and include_waters:
                                dist_hba.set_hbonds_in_selection_and_water_around(
                                    distance
                                )
                            dist_g = dist_hba.filtered_graph
                            nx.write_gpickle(
                                dist_g,
                                self.graph_object_folder
                                + pdb_code
                                + "_"
                                + self.graph_type
                                + "_graphs.pickle",
                            )
                            self.graph_coord_objects[pdb_code].update(
                                {"dist_graph": dist_g}
                            )
                        except:
                            self.logger.warning(
                                f"Graph can't be calculated for {pdb_file}"
                            )

        elif self.type_option == "dcd" and self.graph_type == "water_wire":
            self.selection = selection
            self.max_water = max_water
            self.water_graphs_folder = _hf.create_directory(
                self.graph_object_folder + str(self.max_water) + "_water_wires/"
            )
            self.logger.info(
                f"Maximum number of water in water bridges is set to : {max_water}"
            )
            self.logger.info(f"H-bond criteria cut off distance: {distance} A")
            if check_angle:
                self.logger.info(f"H-bond criteria cut off angle: {cut_angle} degree")
            for i, (name, files) in enumerate(self.graph_coord_objects.items()):
                self.logger.info(
                    f"Loading "
                    + str(len(files["dcd"]))
                    + " trajectory files for "
                    + name
                )
                self.logger.info("This step takes some time.")
                wba = mdh.WireAnalysis(
                    self.selection,
                    files["psf"],
                    files["dcd"],
                    # residuewise=False,
                    residuewise=True,
                    check_angle=check_angle,
                    add_donors_without_hydrogen=not check_angle,
                    additional_donors=additional_donors,
                    additional_acceptors=additional_acceptors,
                    distance=distance,
                    cut_angle=cut_angle,
                )
                wba.set_water_wires(water_in_convex_hull=max_water, max_water=max_water)
                wba.compute_average_water_per_wire()
                self.graph_coord_objects[name].update({"selection": self.selection})
                _hf.pickle_write_file(
                    self.helper_files_folder
                    + name
                    + "_"
                    + str(self.max_water)
                    + "_water_wires_coord_objects.pickle",
                    {name: self.graph_coord_objects[name]},
                )
                g = wba.filtered_graph
                self.graph_coord_objects[name].update({"graph": g})
                tmp = copy.copy(wba)
                tmp._universe = None
                self.graph_coord_objects[name].update({"wba": tmp})
                u = _mda.Universe(files["psf"], files["dcd"])
                mda = u.select_atoms(self.selection)
                self.graph_coord_objects[name].update({"mda": mda})
                wba_loc = (
                    self.water_graphs_folder
                    + name
                    + "_"
                    + str(self.max_water)
                    + "_water_wires_graph.pickle"
                )
                wba.dump_to_file(wba_loc)
                nx.write_gpickle(
                    g,
                    self.helper_files_folder
                    + name
                    + "_"
                    + self.graph_type
                    + "_"
                    + str(max_water)
                    + "_water_nx_graphs.pickle",
                )
                edge_info = _hf.edge_info(wba, g.edges)
                _hf.json_write_file(
                    self.helper_files_folder
                    + name
                    + "_"
                    + self.graph_type
                    + "_"
                    + str(max_water)
                    + "_water_graph_edge_info.json",
                    edge_info,
                )
                self.add_reference_from_structure(mda, g)
                # if i == 0:
                #     self.get_reference_coordinates(mda)
                #     self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)
                self.logger.info("Graph object is saved as: " + wba_loc)

        else:
            raise ValueError(
                'For dcd analysis only graph_type="water_wire" is supported.'
            )
        self.pca_positions = _hf.calculate_pca_positions(self.reference_coordinates)

    def _get_node_positions(self, objects, pca=True):
        node_pos = {}
        for node in objects["graph"].nodes:
            n = _hf.get_node_name(node)
            if (
                n not in self.reference_coordinates.keys()
                or n.split("-")[1] in _hf.water_types
            ):
                chain_id, res_name, res_id = _hf.get_node_name_pats(n)
                if self.type_option == "pdb":
                    chain = objects["structure"].select_atoms("segid " + chain_id)
                    if res_name in _hf.water_types:
                        coords = _hf.get_water_coordinates(chain, res_id)
                    elif res_name in _hf.amino_d.keys():
                        coords = chain.select_atoms(
                            f"resname BWX or protein and name CA and resid {res_id}"
                        ).positions[0]
                    else:
                        coords = chain.select_atoms(
                            f"resname {res_name} and resid {res_id}"
                        ).positions[0]
                elif self.type_option == "dcd":
                    coords = objects["mda"].select_atoms("resid " + res_id).positions[0]
                if coords is not None:
                    node_pos.update({n: list(coords)})
            else:
                node_pos.update({n: self.reference_coordinates[n]})
        if pca:
            return _hf.calculate_pca_positions(node_pos)
        else:
            return node_pos

    def _get_edge_distance(self, objects, name, folder, denumber_waters=False):
        bond_distances = {}

        for edge in objects["dist_graph"].edges:
            e1_chain_id, e1_res_name, e1_res_id, e1_group_name = _hf.get_node_name_pats(
                edge[0], with_group=True
            )
            e2_chain_id, e2_res_name, e2_res_id, e2_group_name = _hf.get_node_name_pats(
                edge[1], with_group=True
            )

            chain_1 = objects["structure"].select_atoms("segid " + e1_chain_id)
            chain_2 = objects["structure"].select_atoms("segid " + e2_chain_id)

            n1 = chain_1.select_atoms(
                f"resname {e1_res_name} and resid {e1_res_id} and name {e1_group_name}"
            ).positions
            n2 = chain_2.select_atoms(
                f"resname {e2_res_name} and resid {e2_res_id} and name {e2_group_name}"
            ).positions

            # assert len(n1) == 1
            # assert len(n2) == 1
            if len(n1) and len(n2):
                dist_arr = distances.distance_array(n1[0], n2[0])

                if denumber_waters:
                    e0 = f"{e1_chain_id}-{e1_res_name}-{e1_res_id}"
                    e1 = f"{e2_chain_id}-{e2_res_name}-{e2_res_id}"
                    if edge[0].split("-")[1] in _hf.water_types:
                        e0 = "HOH"
                    if edge[1].split("-")[1] in _hf.water_types:
                        e1 = "HOH"

                else:
                    # e0 = edge[0]
                    # e1 = edge[1]
                    e0 = f"{e1_chain_id}-{e1_res_name}-{e1_res_id}"
                    e1 = f"{e2_chain_id}-{e2_res_name}-{e2_res_id}"

                if f"{e0}_{e1}" in bond_distances.keys():
                    d = np.mean([bond_distances[f"{e0}_{e1}"], dist_arr[0][0]])
                    bond_distances.update({f"{e0}_{e1}": d})
                elif f"{e1}_{e0}" in bond_distances.keys():
                    d = np.mean([bond_distances[f"{e1}_{e0}"], dist_arr[0][0]])
                    bond_distances.update({f"{e1}_{e0}": d})
                else:
                    bond_distances.update({f"{e0}_{e1}": dist_arr[0][0]})

            else:
                print(
                    f"Skip distance for {e1_res_name}-{e1_res_id} and {e2_res_name}-{e2_res_id} as one of the residues not found in the PDB file."
                )

        _hf.write_text_file(
            folder + name + "_H-bond_graph_distances.txt",
            [f"{edge} {distance}\n" for edge, distance in bond_distances.items()],
        )
        return bond_distances

    def plot_graphs(
        self,
        label_nodes=True,
        label_edges=True,
        xlabel="PCA projected xy plane",
        ylabel="Z coordinates ($\AA$)",
        occupancy=None,
        color_propka=False,
        color_data=False,
        node_color_selection=None,
        node_color_map="viridis",
        calcualte_distances=False,
        color_bfactor=False,
        color_edges=False,
    ):
        if occupancy is not None or hasattr(self, "occupancy"):
            if occupancy is not None:
                occupancy = occupancy
            else:
                occupancy = self.occupancy
        for name, objects in self.graph_coord_objects.items():
            if "graph" in objects.keys():
                if occupancy:
                    wba = objects["wba"]
                    # wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(occupancy)
                    graph = wba.filtered_graph
                else:
                    graph = objects["graph"]
                self.logger.debug("Creating " + self.graph_type + " plot for: " + name)
                plot_name = "H-bond" if self.graph_type == "hbond" else "Water wire"
                fig, ax = _hf.create_plot(
                    title=f"{plot_name} graph of structure {name}\nSelection:{self.selection[1:-16]}",
                    xlabel=xlabel,
                    ylabel=ylabel,
                    plot_parameters=self.plot_parameters,
                )
                node_pca_pos = self._get_node_positions(objects)
                node_pca_pos = _hf.check_projection_sign(
                    node_pca_pos, self.pca_positions
                )

                f = _hf.create_directory(
                    self.workfolder + "/H-bond_graphs/" + name + "/"
                )
                if calcualte_distances and self.graph_type == "hbond":
                    bond_distances = self._get_edge_distance(objects, name, f)

                edge_colors_data = {}
                if color_edges:  # currently only implemented for dcd
                    edge_value_dict = _hf.read_edge_color_data(name, objects["psf"])
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

                        if label_edges and self.graph_type == "water_wire":
                            waters, occ_per_wire, _ = _hf.get_edge_params(
                                objects["wba"], graph.edges
                            )
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

                        elif (
                            calcualte_distances
                            and label_edges
                            and self.graph_type == "hbond"
                        ):
                            if f"{e[0]}_{e[1]}" in bond_distances:
                                dist = bond_distances[f"{e[0]}_{e[1]}"]
                            elif f"{e[1]}_{e[0]}" in bond_distances:
                                dist = bond_distances[f"{e[1]}_{e[0]}"]

                            ax.annotate(
                                round(dist, 1),
                                (
                                    x[0] + (x[1] - x[0]) / 2,
                                    y[0] + (y[1] - 1.0 - y[0]) / 2,
                                ),
                                color="blue",
                                fontsize=self.plot_parameters["edge_label_size"],
                            )
                color_info = {}
                lab = " with labels" if label_nodes else ""
                if color_propka and color_data:
                    self.logger.info(
                        f"Can not color plot by propka and external data values at the same time. Please select just one coloring option!"
                    )
                elif color_propka:
                    struct_object = (
                        objects["structure"]
                        if self.type_option == "pdb"
                        else objects["mda"]
                    )
                    selected_nodes = struct_object.select_atoms(
                        str(node_color_selection)
                    )
                    try:
                        color_info = _hf.read_propka_file(
                            f"{self.pdb_root_folder}/{name}.propka", selected_nodes
                        )
                    except:
                        self.logger.info(
                            f"{name}.propka not found. To color residues by pKa values, place the propka file in the PDB folder, next to the PDB file."
                        )
                    try:
                        len(color_info)
                        value_colors, cmap, norm = _hf.get_color_map(
                            color_info, color_map=node_color_map
                        )
                        self.logger.info(f"Color {name} by pKa values{lab}.")
                        color_bar_label = "pKa value"
                    except:
                        self.logger.info(
                            f"{name}.propka does not contain the selected residues. Please update the Residues to color!"
                        )

                elif color_data:
                    struct_object = (
                        objects["structure"]
                        if self.type_option == "pdb"
                        else objects["mda"]
                    )
                    selected_nodes = struct_object.select_atoms(
                        str(node_color_selection)
                    )
                    color_info = _hf.read_color_data_file(
                        name, self.pdb_root_folder, selected_nodes
                    )
                    # breakpoint()
                    try:
                        len(color_info)
                        value_colors, cmap, norm = _hf.get_color_map(
                            color_info, color_map=node_color_map
                        )
                        color_bar_label = "Amino acid data value"
                        self.logger.info(
                            f"Color {name} by values from external data file{lab}."
                        )
                    except:
                        self.logger.error(
                            f"No {name}_data.txt file was found in {self.pdb_root_folder} or content is invalid. To enable coloring by data values please add a corresponding file."
                        )

                elif color_bfactor:
                    try:
                        self.type_option == "pdb"
                    except:
                        self.logger.error(
                            "B-factor coloring only supported for PDB files"
                        )

                    selected_nodes = objects["structure"].select_atoms(
                        str(node_color_selection)
                    )
                    try:
                        bfactors = selected_nodes.tempfactors
                        color_bar_label = "B-factor"
                        for i, resisdue in enumerate(selected_nodes):
                            chain, res_name, res_id = (
                                resisdue.segid,
                                resisdue.resname,
                                resisdue.resid,
                            )
                            res = chain + "-" + res_name + "-" + str(res_id)
                            color_info.update({res: bfactors[i]})

                        value_colors, cmap, norm = _hf.get_color_map(
                            color_info, color_map=node_color_map
                        )
                    except:
                        self.logger.error("No B-factors were available")

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
                            # color = self.plot_parameters['graph_color']
                            ax.scatter(
                                values[0],
                                values[1],
                                color=color,
                                s=self.plot_parameters["node_size"],
                                zorder=5,
                                edgecolors=self.plot_parameters["graph_color"],
                            )
                            # ax.scatter(values[0],values[1], s=self.plot_parameters['node_size'], zorder=5,  c=df.z, cmap='Greens_r')
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
                                pass  # CHECK: temporary turn off water labels
                                # ax.annotate(f'W{res_id}', (values[0]+0.2, values[1]-0.25), fontsize=self.plot_parameters['node_label_size'])
                            elif res_name in _hf.amino_d.keys():
                                l = (
                                    f"{chain_id}-{_hf.amino_d[res_name]}{res_id}"
                                    if self.plot_parameters["show_chain_label"]
                                    else f"{_hf.amino_d[res_name]}{int(res_id)+19}"
                                )
                                bbox = dict(
                                    facecolor="white", alpha=0.3, edgecolor="white"
                                )  # CHECK: add white transparent box behind labels for visibility
                                ax.annotate(
                                    l,
                                    (values[0] + 0.2, values[1] - 0.26),
                                    fontsize=self.plot_parameters["node_label_size"],
                                )  # , bbox=bbox)
                            else:
                                l = (
                                    f"{chain_id}-{res_name}{res_id}"
                                    if self.plot_parameters["show_chain_label"]
                                    else f"{res_name}{int(res_id)+19}"
                                )
                                ax.annotate(
                                    l,
                                    (values[0] + 0.2, values[1] - 0.25),
                                    fontsize=self.plot_parameters["node_label_size"],
                                    color=self.plot_parameters["non_prot_color"],
                                )

                if color_info:
                    cbar = fig.colorbar(
                        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax
                    )
                    cbar.ax.tick_params(
                        labelsize=self.plot_parameters["plot_tick_fontsize"]
                    )
                    cbar.set_label(
                        label=color_bar_label,
                        size=self.plot_parameters["plot_label_fontsize"],
                    )

                plt.tight_layout()
                is_label = "_labeled" if label_nodes else ""
                is_propka = "_pKa_color" if color_info and color_propka else ""
                is_bfactor = "_bfactor_color" if color_info and color_bfactor else ""
                is_conservation = "_data_color" if color_info and color_data else ""
                is_distance = "_distance" if calcualte_distances else ""
                is_backbone = (
                    "_backbone"
                    if hasattr(self, "include_backbone_sidechain")
                    and self.include_backbone_sidechain
                    else ""
                )
                is_water = (
                    "_no_water"
                    if hasattr(self, "include_waters") and not self.include_waters
                    else ""
                )
                if self.graph_type == "hbond":
                    plot_folder = _hf.create_directory(
                        self.workfolder + "/H-bond_graphs/" + name + "/"
                    )
                    if color_bfactor:
                        _hf.write_text_file(
                            plot_folder + name + "_H-bond_graph_bfactors.txt",
                            [
                                f"{node} {round(color_info[node], 2)}\n"
                                for node in graph.nodes
                                if node in color_info.keys()
                            ],
                        )

                    for form in self.plot_parameters["formats"]:
                        plt.savefig(
                            f"{plot_folder}{name}_H-bond_graph{is_propka}{is_bfactor}{is_conservation}{is_distance}{is_backbone}{is_water}{is_label}.{form}",
                            format=form,
                            dpi=self.plot_parameters["plot_resolution"],
                        )
                    if is_label:
                        _hf.write_text_file(
                            plot_folder
                            + name
                            + is_backbone
                            + is_water
                            + "_H-bond_graph_info.txt",
                            [
                                "H-bond graph of " + name,
                                "\nSelection string: " + str(self.selection[0:-15]),
                                "\n",
                                "\nNumber of nodes in "
                                + name
                                + ": "
                                + str(len(graph.nodes)),
                                "\nNumber of edges in "
                                + name
                                + ": "
                                + str(len(graph.edges)),
                                "\n",
                                "\nList of nodes: " + str(graph.nodes),
                                "\n",
                                "\nList of edges: " + str(graph.edges),
                            ],
                        )
                elif self.graph_type == "water_wire":
                    plot_folder = _hf.create_directory(
                        self.workfolder
                        + "/"
                        + str(self.max_water)
                        + "_water_wires/"
                        + name
                        + "/"
                    )

                    waters = (
                        "_max_" + str(self.max_water) + "_water_bridges"
                        if self.max_water > 0
                        else ""
                    )
                    occ = "_min_occupancy_" + str(occupancy) if occupancy else ""
                    for form in self.plot_parameters["formats"]:
                        plt.savefig(
                            f"{plot_folder}{name}{waters}{occ}_graph{is_propka}{is_conservation}{is_backbone}{is_label}.{form}",
                            format=form,
                            dpi=self.plot_parameters["plot_resolution"],
                        )
                    if is_label:
                        _hf.write_text_file(
                            plot_folder
                            + name
                            + waters
                            + occ
                            + "_water_wire_graph_info.txt",
                            [
                                "Water wire graph of " + name,
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
                                + name
                                + ": "
                                + str(len(graph.nodes)),
                                "\nNumber of edges in "
                                + name
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
            objects["structure"] if self.type_option == "pdb" else objects["mda"]
        )
        return _hf.calculate_connected_compontents_coordinates(
            connected_components, struct_object, option=self.type_option
        )

    def plot_linear_lenghts(self, occupancy=None, label_nodes=True):
        self.logger.info("Plotting linear lengths for continuous network components")
        for name, objects in self.graph_coord_objects.items():
            self.logger.debug("Creating linear length plot for " + name)
            if "graph" in objects.keys():
                if occupancy:
                    wba = objects["wba"]
                    # wba = copy.deepcopy(objects['wba'])
                    wba.filter_occupancy(occupancy)
                    graph = wba.filtered_graph
                else:
                    graph = objects["graph"]
                self.logger.debug(
                    "Creating " + self.graph_type + " linear length plot for: " + name
                )
                connected_components_coordinates = self.get_linear_lenght(
                    objects, graph
                )
                plot_parameters = copy.deepcopy(self.plot_parameters)
                plot_parameters["figsize"] = (
                    min(1 + int(len(connected_components_coordinates)), 100),
                    plot_parameters["figsize"][1],
                )

                plot_name = "H-bond" if self.graph_type == "hbond" else "water wire"
                fig, ax = _hf.create_plot(
                    title=f"Linear length of continuous {plot_name} subnetworks \nalong the Z-axis in structure {name}\nSelection: {self.selection[1:-16]}",
                    xlabel="# of nodes",
                    ylabel="Z-axis coordinates ($\AA$)",
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
                                l = (
                                    f"{chain_id}-{_hf.amino_d[res_name]}{res_id}"
                                    if self.plot_parameters["show_chain_label"]
                                    else f"{_hf.amino_d[res_name]}{res_id}"
                                )
                                ax.annotate(
                                    l,
                                    (i, z_coords),
                                    fontsize=self.plot_parameters["node_label_size"],
                                    zorder=6,
                                )
                            else:
                                l = (
                                    f"{chain_id}-{res_name}{res_id}"
                                    if self.plot_parameters["show_chain_label"]
                                    else f"{res_name}{res_id}"
                                )
                                ax.annotate(
                                    l,
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
                is_water = (
                    "_no_water"
                    if hasattr(self, "include_waters") and not self.include_waters
                    else ""
                )

                if self.graph_type == "hbond":
                    plot_folder = _hf.create_directory(
                        self.workfolder + "/H-bond_graphs/" + name + "/"
                    )
                    for form in self.plot_parameters["formats"]:
                        plt.savefig(
                            f"{plot_folder}{name}_H-bond_linear_length{is_backbone}{is_water}{is_label}.{form}",
                            format=form,
                            dpi=self.plot_parameters["plot_resolution"],
                        )
                elif self.graph_type == "water_wire":
                    plot_folder = _hf.create_directory(
                        self.workfolder
                        + "/"
                        + str(self.max_water)
                        + "_water_wires/"
                        + name
                        + "/"
                    )
                    waters = (
                        "_" + str(self.max_water) + "_water_bridges"
                        if self.max_water > 0
                        else ""
                    )
                    occ = "_min_occupancy_" + str(occupancy) if occupancy else ""
                    for form in self.plot_parameters["formats"]:
                        plt.savefig(
                            f"{plot_folder}{name}{waters}{occ}_linear_length{is_backbone}{is_water}{is_label}.{form}",
                            format=form,
                            dpi=self.plot_parameters["plot_resolution"],
                        )
                plt.close()
            else:
                self.logger.warning(
                    f"{name} has no {self.graph_type} graph. Linear length can not be calculated for this structure."
                )
