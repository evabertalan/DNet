import streamlit as st
import time
import re
import json


def generate_bash(env, glob, pKa, dist, plot, graphs):
    # header and set up module if script is run on e.g: JURECA vs. locally
    bash_script = f"""#!/bin/bash
if command -v module &>/dev/null; then
    echo "Module system detected. Loading Python..."
    module load Python
else
    echo ""
fi

source {env['venv_path']}

# Setup Logging
OUTPUT_FOLDER="{env['output_folder']}"
mkdir -p "$OUTPUT_FOLDER"
LOGFILE="$OUTPUT_FOLDER/dnet_{env['dnet_run_name']}.log"

echo "===============================================================" >> "$LOGFILE"
echo "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet Workflow ..." >> "$LOGFILE"

# Global Files
PSF_FILE="{env['psf']}"
DCD_FILES=$(ls {env['dcd']})
WRAP_DCD={env['wrap_dcd']}
psf_name_no_ext=$(basename "$PSF_FILE" | cut -d. -f1)
cd {env['dnet_dir']} || exit

#Global Parameters
SELECTION="{glob["selection"]}"
DISTANCE_CUT_OFF={glob["distance_cut_off"]}
ANGLE_CUT_OFF={glob["angle_cut_off"]}
OCCUPANCY={glob["occupancy"]}
MAX_WATER={glob["max_water"]}
DONORS="{glob["donors"]}"
ACCEPTORS="{glob["acceptors"]}"
STEP={glob["step"]}
PLOT_PARAMETERS="{glob["plot_parameters"]}"

"""

    if pKa["run_pKa"]:
        bash_script += "\n# --- Module: DNet-pKa ---"
        bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet-pKa calculation..." >> "$LOGFILE"\n"""
        bash_script += "\nPKA_FOLDER=${OUTPUT_FOLDER}/pKa"
        bash_script += "\nmkdir -p $PKA_FOLDER"
        bash_script += f'''\npython3 -m dnet_pka "$PSF_FILE" $DCD_FILES \\
                        --output_folder "$PKA_FOLDER" \\
                        --selection "{pKa["pKa_selection"]}"'''
        if pKa["pKa_start"]:
            bash_script += f' --start {pKa["pKa_start"]}'
        if pKa["pKa_stop"]:
            bash_script += f' --stop {pKa["pKa_stop"]}'
        if pKa["pKa_step"]:
            bash_script += f' --step {pKa["pKa_step"]}'
        bash_script += ' >> "$LOGFILE" 2>&1 &\n'

    bash_script += "\n# --- Module: DNet-Graphs atomwise ---"
    bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet-Graphs atomwise calculation..." >> "$LOGFILE"\n"""
    bash_script += "\nATOMWISE_FOLDER=${OUTPUT_FOLDER}/atomwise_graphs"
    bash_script += "\nmkdir -p $ATOMWISE_FOLDER"
    bash_script += f'''\npython3 -m dnet_graphs "$PSF_FILE" $DCD_FILES \\
                        --atomwise \\
                        --output_folder "$ATOMWISE_FOLDER" \\
                        --cut_angle "$ANGLE_CUT_OFF" \\
                        --distance "$DISTANCE_CUT_OFF" \\
                        --max_water "$MAX_WATER" \\
                        --occupancy "$OCCUPANCY" \\
                        --step {glob["step"]} \\
                        --selection "$SELECTION" \\
                        --plot_parameters "$PLOT_PARAMETERS" \\
                        --additional_donors "$DONORS" --additional_acceptors "$ACCEPTORS"'''
    if glob["start"]:
        bash_script += f' --start {glob["start"]}'
    if glob["stop"]:
        bash_script += f' --stop {glob["stop"]}'
    if glob["no_label_plots"]:
        bash_script += f" --no_label_plots"
    if glob["dont_save_graph_objects"]:
        bash_script += f" --dont_save_graph_objects"
    if glob["collect_angles"]:
        bash_script += f" --collect_angles"
    if glob["res_id_label_shift"]:
        bash_script += f""" --res_id_label_shift '{glob["res_id_label_shift"]}'"""

    bash_script += ' >> "$LOGFILE" 2>&1 &\n'
    bash_script += "\nwait\n"

    if dist["run_dist"]:
        bash_script += "\n# --- Module: DNet-Dist ---"
        bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet-Dist calculation..." >> "$LOGFILE"\n"""
        bash_script += "\nDISTANCE_FOLDER=${OUTPUT_FOLDER}/distances"
        bash_script += "\nmkdir -p $DISTANCE_FOLDER"
        bash_script += f'\nMAX_WATER_DISTANCE={dist["dist_max_water_distance"]}'
        bash_script += "\nGRAPHS_INPUT=${ATOMWISE_FOLDER}/workfolder/${MAX_WATER}_water_wires/${psf_name_no_ext}/${psf_name_no_ext}_max_${MAX_WATER}_water_bridges_min_occupancy_${OCCUPANCY}_water_wire_graph_info.txt"
        bash_script += f"""\npython3 -m dnet_dist "$PSF_FILE" $DCD_FILES "$GRAPHS_INPUT"\\
                        --output_folder "$DISTANCE_FOLDER" \\
                        --step {dist["dist_step"]} \\
                        --max_water_distance {dist["dist_max_water_distance"]}"""
        if dist["dist_start"]:
            bash_script += f' --start {dist["dist_start"]}'
        if dist["dist_stop"]:
            bash_script += f' --stop {dist["dist_stop"]}'

        bash_script += ' >> "$LOGFILE" 2>&1\n'

    if plot["plot_run"]:
        bash_script += "\n# --- Module: DNet-Plot ---"
        bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet-Plot calculation..." >> "$LOGFILE"\n"""
        bash_script += "\nPLOT_FOLDER=${OUTPUT_FOLDER}/plots"
        bash_script += "\nmkdir -p $PLOT_FOLDER"
        bash_script += "\nGRAPHS_INFO_TXT=${ATOMWISE_FOLDER}/workfolder/${MAX_WATER}_water_wires/${psf_name_no_ext}/${psf_name_no_ext}_max_${MAX_WATER}_water_bridges_min_occupancy_${OCCUPANCY}_water_wire_graph_info.txt"
        bash_script += (
            "\nPKAS_FOR_FRAME_CSV=${PKA_FOLDER}/pkas_for_frames_${psf_name_no_ext}.csv"
        )
        bash_script += "\nPAIR_DISTANCES_CSV=${DISTANCE_FOLDER}/${psf_name_no_ext}_pair_distances.csv"
        bash_script += "\nWATER_WITHIN_CSV=${DISTANCE_FOLDER}/${psf_name_no_ext}_waters_within_${MAX_WATER_DISTANCE//./_}_of_group.csv"
        bash_script += "\nTOTAL_WATER_WITHIN_CSV=${DISTANCE_FOLDER}/${psf_name_no_ext}_total_waters_within_${MAX_WATER_DISTANCE//./_}_of_res.csv"
        bash_script += f"""\npython3 -m dnet_plot --plot_folder "$PLOT_FOLDER" \\
                        --graphs_info_txt "$GRAPHS_INFO_TXT" \\
                        --pKas_for_frame_csv "$PKAS_FOR_FRAME_CSV" \\
                        --pair_distances_csv "$PAIR_DISTANCES_CSV" \\
                        --water_within_csv "$WATER_WITHIN_CSV" \\
                        --total_water_within_csv "$TOTAL_WATER_WITHIN_CSV" \\
                        --frame_to_time {plot["plot_frame_to_time"]} \\
                        --pmf_last_nth_frames {plot["pmf_last_nth_frames"]}"""
        if glob["res_id_label_shift"]:
            bash_script += f""" --res_id_label_shift '{glob["res_id_label_shift"]}'"""
        bash_script += ' >> "$LOGFILE" 2>&1\n'

    if len(graphs):
        bash_script += "\n# --- Module: DNet-Graphs ---"
        bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Starting DNet-Graph calculation..." >> "$LOGFILE"\n"""

        for graph in graphs:
            graph_type = "--residuewise"
            bash_script += (
                f"\n# --- Running graph calculation {graph['calcualtion_name']} ---"
            )
            bash_script += (
                "\nRES_FOLDER=${OUTPUT_FOLDER}/" + f"{graph['calcualtion_name']}"
            )
            bash_script += "\nmkdir -p $RES_FOLDER"

            if graph["nodes_colored_by"] == "Most frequently sampled pKa value":
                bash_script += (
                    "\ncp $PKA_FOLDER/${psf_name_no_ext}_pKa_data.txt $RES_FOLDER"
                )
                graph_type = "--residuewise"

            elif (
                graph["nodes_colored_by"]
                == "Avg. number of water molecules around the amino acid side chain"
            ):
                bash_script += "\ncp $DISTANCE_FOLDER/${psf_name_no_ext}_avg_waters_data.txt $RES_FOLDER"
                graph_type = "--residuewise"
                pass  # implement avg water coloring; copy file to the folder

            if (
                graph["edges_colored_by"]
                == "PN: population number (estimated number of conformations sampled by the connection)"
            ):
                graph_type = "--atomwise"
                bash_script += "\nEDGE_COLOR_FILE=${PLOT_FOLDER}/${psf_name_no_ext}_PN_per_edges.txt"

            if graph["component_search"] == "Path search between start and goal nodes":
                bash_script += (
                    f"""\nPATH_NAME=('{graph["start_node"]}' '{graph["goal_node"]}')"""
                )

            bash_script += f'''\npython3 -m dnet_graphs "$PSF_FILE" $DCD_FILES --output_folder "$RES_FOLDER" \\
                            {graph_type} \\
                            --step {glob["step"]} \\
                            --cut_angle {graph["angle_cut_off"]} \\
                            --distance {graph["distance_cut_off"]} \\
                            --max_water {graph["max_water"]} \\
                            --occupancy {graph["min_occupancy"]} \\
                            --selection "{graph['selection']}"\\
                            --plot_parameters "$PLOT_PARAMETERS" \\
                            --additional_donors "$DONORS" \\
                            --additional_acceptors "$ACCEPTORS"'''
            if glob["start"]:
                bash_script += f' --start {glob["start"]}'
            if glob["stop"]:
                bash_script += f' --stop {glob["stop"]}'
            if glob["no_label_plots"]:
                bash_script += f" --no_label_plots"
            if glob["dont_save_graph_objects"]:
                bash_script += f" --dont_save_graph_objects"
            if glob["collect_angles"]:
                bash_script += f" --collect_angles"
            if glob["res_id_label_shift"]:
                bash_script += (
                    f""" --res_id_label_shift '{glob["res_id_label_shift"]}'"""
                )
            if graph["nodes_colored_by"]:
                bash_script += f" --color_data"
            if (
                graph["nodes_colored_by"]
                == "Avg. number of water molecules around the amino acid side chain"
            ):
                bash_script += f' --node_color_map "Magma"'
            if graph["edges_colored_by"] == "Occupancy":
                bash_script += f" --color_edge_by_occupnacy"
            if (
                graph["edges_colored_by"]
                == "PN: population number (estimated number of conformations sampled by the connection)"
            ):
                bash_script += f' --color_edges_by "$EDGE_COLOR_FILE"'
            if graph["component_search"] == "Connected component of a root node":
                bash_script += f''' --root "{graph['root_node']}"'''
            if graph["component_search"] == "Path search between start and goal nodes":
                bash_script += ' --path "${PATH_NAME[@]}"'

            bash_script += ' >> "$LOGFILE" 2>&1 &\n'

    bash_script += "\nwait\n"
    bash_script += """\necho "$(date '+%Y-%m-%d %H:%M:%S'): Finished DNet calculation" >> "$LOGFILE"\n"""
    bash_script += "\ndeactivate"

    return bash_script


st.set_page_config(page_title="DNet Run File Generator", layout="wide")
st.title("DNet Run File Generator")


# col_env, col_params = st.columns([1, 2])

# with col_env:
st.header("Environment parameters")
st.info(
    "Only change the paths below if the run file is not located in the folder where DNet is installed."
)
c1, c2 = st.columns([1, 1])
with c1:
    venv = st.text_input(
        "Venv Path",
        "./dnet_env/bin/activate",
        help="If the run file is in the folder where DNet is installed, don't change this path",
    )
with c2:
    dnet_dir = st.text_input(
        "Path to the folder where DNET is installed, relative to this run file",
        "./DNet",
        help="If the run file is next to the DNet folder, don't change this path",
    )

st.divider()

st.header("Global parameters")

output_folder = st.text_input(
    "Output folder - required",
    help="Path to the folder where the log file and all results will be saved.",
)
if output_folder:
    st.caption(f"**Location of the outputs:** `{output_folder}`")
else:
    st.caption("*No output folder is selected.*")

dnet_run_name = st.text_input(
    "Name of the calculation:",
    value="dnet_calculation",
    help="The following name will be the name of the log file and this run file.",
)

dnet_run_name = re.sub(r"[^a-zA-Z0-9\s\-_]", "", dnet_run_name).replace(" ", "_")

if not dnet_run_name:
    dnet_run_name = f"dnet_calculation"

st.caption(f"**The log file will be saved in:** `{output_folder}/{dnet_run_name}.log`")


psf = st.text_input(
    "PSF File - required",
    help="Path to the PSF file required to load the molecular system.",
)
if psf:
    st.caption(f"**Selected .psf file:** `{psf}`")
else:
    st.caption("*No PSF file selected.*")

dcd = st.text_input(
    "DCD File(s) - required",
    help="Path(s) to trajectory files (supports wildcards). E.g: `path_to_dcd/9cis_m103a*.dcd | grep -v 'PBC.dcd'` --> will select all dcd files starting with 9cis_m103a, that does NOT end with PBC.dcd",
)
if dcd:
    st.caption(f"**Selected .dcd files:** `{dcd}`")
else:
    st.caption("*No DCD files are selected.*")

wrap_dcd = st.checkbox(
    "PBC wrap the trajectories",
    value=True,
    help="Apply periodic boundary condition wrapping to keep molecules inside the simulation box. Default is true.",
)


selection = st.text_input(
    "Selection to perform the analysis on, in a form of an MDAnalysis selection sting. ",
    "protein",
    help="Default is protein. Documentation of the selection language: `https://userguide.mdanalysis.org/stable/selections.html`",
)


st.subheader("Trajectory reading steps")

st.info(
    "The following variables set the start, stop and step size globally, but it can later be changed individually for each module."
)
c1, c2, c3 = st.columns([1, 1, 1])

with c1:
    start = c1.number_input(
        "Start",
        value=None,
        step=1,
        help="Starting frame index for trajectory analysis. If not provided, starts from the first frame.",
        key="start",
    )

with c2:
    stop = c2.number_input(
        "Stop",
        value=None,
        step=1,
        help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -20000 reads the last 20000 frames of the trajectory.",
        key="stop",
    )

with c3:
    step = c3.number_input(
        "Step",
        1,
        step=1,
        help="Step size for reading the trajectory frames. For example if its set to 10, only  every 10th frame is read, which reduces the computation time.",
        key="step",
    )

if not output_folder or not psf or not dcd:
    st.warning(
        "Please first fill out the required fields. Input fields for the additional parameters will only be shown after the required parameters are set up."
    )
    st.stop()

st.divider()

st.header("Module specific parameters")
st.subheader("Global H-bond network criteria")
st.info(
    "These parameters set up the global H-bond criteria for the H-bond calculation in the `DNet-dist` and `DNet-graphs` modules. Additional graph calculations with different criteria can be added under the `DNet-graphs` tab."
)


def H_bond_critera_element(key_index):

    c1, c2, c3, c4 = st.columns([1, 1, 1, 1])

    with c1:
        distance_cut_off = c1.number_input(
            "Distance cut off",
            value=3.5,
            step=0.1,
            min_value=0.1,
            max_value=15.0,
            key=f"distance_{key_index}",
            help="The distance criterion for the  H-bond search, measured between the heavy atoms. The default value is 3.5Å.",
        )

    with c2:
        angle_cut_off = c2.number_input(
            "Angle cut off",
            value=60,
            step=1,
            min_value=0,
            max_value=180,
            key=f"angle_{key_index}",
            help="Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. The default value is 60°.",
        )

    with c3:
        occupancy = c3.number_input(
            "Min H-bond occupancy",
            value=0.1,
            step=0.01,
            min_value=0.01,
            max_value=1.0,
            key=f"occupancy_{key_index}",
            help="Minimum H-bond occupancy required to include an edge in the graph (default: 0.1, which means 10% occupancy).",
        )
    with c4:
        max_water = c4.number_input(
            "Max number of waters allowed in the bridge",
            value=3,
            step=1,
            min_value=0,
            max_value=5,
            key=f"max_water_{key_index}",
            help="Maximum number of water molecules allowed in the water wire connections (default: 3). When it is set to 0, only direct H-bonds are considered.",
        )

    return distance_cut_off, angle_cut_off, occupancy, max_water


distance_cut_off, angle_cut_off, occupancy, max_water = H_bond_critera_element(0)


c1, c2 = st.columns([1, 1])

donors = c1.text_input(
    "Additional Donors",
    "",
    help="List of additional hydrogen bond donor atoms separated by a `,` e.g: `N, S`.",
)
donors = str([atom.strip() for atom in donors.split(",") if atom.strip()])


acceptors = c2.text_input(
    "Additional Acceptors",
    "",
    help="List of additional hydrogen bond donor atoms separated by a `,` e.g: `O, F`.",
)
acceptors = str([atom.strip() for atom in acceptors.split(",") if atom.strip()])


# c1, c2, c3 = st.columns([1, 1, 1])
# with c1:
no_label_plots = st.checkbox(
    "Generate additional plots without the labels",
    value=True,
    help="Creates all the graph plots without the nodes and edges labels as well. Useful for preparation additional figures, where the graph nodes need to be re-labeled or only a few graph nodes need to be labeled.",
)

# with c2:
dont_save_graph_objects = st.checkbox(
    "Don't save large data objects from the calculation",
    value=True,
    help="Don't save the metadata and full graph objects of the calculations. Use this flag if there is not enough space for the calculation results or when the graph objects are not needed for further calculations or analysis.",
)

# with c3:
collect_angles = st.checkbox(
    "Collect H-bond angles for additional analysis later",
    value=False,
    help="Create a csv file with the angles of all donor-acceptor pairs that are within the set H-bond distance criterion in each frame.",
)

res_id_label_shift = {}
shift_reid_labels = st.checkbox(
    "Shift residue ID labels by a given offset",
    help="One value per protein segment can to be provided.",
    value=False,
)

if shift_reid_labels:
    st.warning(
        "Use it with caution! Residue ID labels will have an offset only on the plots. In the additional data files residues are numbered according to the structure file. For root and path search the start and nodes have to be given according to their numbering in the structure file."
    )
    res_id_label_shift_input = st.text_area(
        "Enter segment_name: offset_value pairs (one per line). Example:  `PROA: 10`",
        value="segment1: 10\nsegment2: 20",
    )

    lines = res_id_label_shift_input.split("\n")

    for line in lines:
        if ":" in line:
            try:
                key, val = line.split(":", 1)
                key = key.strip()
                val = int(val.strip())
                res_id_label_shift.update({key: val})

            except ValueError:
                st.write(
                    f"Skipping {line}. Incorrect format. Please provide in a format of `PROA: 10`"
                )
                continue
    res_id_label_shift = json.dumps(res_id_label_shift)

with st.expander("Adjust Graph Plot Visualization", expanded=False):
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### Sizes")
        edge_w = st.number_input("Edge Width", min_value=1, value=2, step=1)
        node_s = st.number_input("Node Size", min_value=1, value=15, step=1)
        node_lbl_s = st.number_input("Node Label Size", min_value=1, value=15, step=1)
        edge_lbl_s = st.number_input("Edge Label Size", min_value=1, value=15, step=1)

        st.markdown("#### Colors")
        g_color = st.color_picker("Graph Color", "#808080")
        # w_color = st.color_picker("Water Node Color", "#db5c5c")
        # d_color = st.color_picker("Difference Graph Color", "#129fe6")
        np_color = st.color_picker("Color of Non-Protein Nodes", "#0000ff")

    with col2:
        st.markdown("#### Layout & Format")
        title_fs = st.number_input("Title Font Size", min_value=1, value=30)
        label_fs = st.number_input("Label Font Size", min_value=1, value=36)
        tick_fs = st.number_input("Tick Font Size", min_value=1, value=33)
        res = st.number_input("Resolution (DPI)", min_value=72, value=400)

        f_width = st.number_input("Fig Width", min_value=1, value=15)
        f_height = st.number_input("Fig Height", min_value=1, value=16)

        fmts = st.multiselect(
            "Export Formats", ["png", "pdf", "svg"], default=["png", "svg"]
        )
        show_chain = st.checkbox("Show Chain Label", value=False)

plot_parameters = {
    "edge_width": edge_w,
    "node_label_size": node_lbl_s,
    "edge_label_size": edge_lbl_s,
    "node_size": node_s,
    "graph_color": g_color,
    # "water_node_color": w_color,
    # "difference_graph_color": d_color,
    "non_prot_color": np_color,
    "plot_title_fontsize": title_fs,
    "plot_label_fontsize": label_fs,
    "plot_tick_fontsize": tick_fs,
    "plot_resolution": res,
    "figsize": (f_width, f_height),
    "formats": fmts,
    "show_chain_label": show_chain,
}


st.subheader("Set up calculation parameters for the DNet modules")
t1, t2, t3, t4 = st.tabs(
    [
        "DNet-pKa",
        "DNet-Dist",
        "DNet-Plot",
        "DNet-Graphs",
    ]
)

with t1:
    run_pKa = st.checkbox("Run DNet-pKa module", value=True)
    pKa_start = start
    pKa_stop = stop
    pKa_step = step
    pKa_selection = selection

    if run_pKa:
        custom_pKa_steps = st.checkbox(
            "Set up different `start`, `stop` and `step size` parameters for the DNet-pKa module, than the global values",
            value=False,
        )
        if custom_pKa_steps:
            c1, c2, c3 = st.columns([1, 1, 1])

            with c1:
                pKa_start = c1.number_input(
                    "pKa Start",
                    value=None,
                    step=1,
                    help="Starting frame index for trajectory analysis. If not provided, starts from the first frame.",
                    key="pKa_start",
                )

            with c2:
                pKa_stop = c2.number_input(
                    "pKa Stop",
                    value=None,
                    step=1,
                    help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -2000 reads the last 2000 frames of the trajectory.",
                    key="pKa_stop",
                )

            with c3:
                pKa_step = c3.number_input(
                    "pKa Step",
                    value=None,
                    step=1,
                    help="Step size for reading the trajectory frames. For example if its set to 10, only  every 10th frame is read, which reduces the computation time.",
                    key="pKa_steo",
                )

        custom_pKa_selection = st.checkbox(
            f"Perform pKa calculation on a different selection, than the global selection of: `{selection}`",
            value=False,
        )
        if custom_pKa_selection:
            pKa_selection = st.text_input(
                "Selection to perform the pKa calculation in a form of an MDAnalysis selection sting.",
                "protein",
                help="Default is protein Which calculates pKa values for each titratable amino acid side chains. Documentation of the selection language: `https://userguide.mdanalysis.org/stable/selections.html`",
            )


with t2:

    run_dist = st.checkbox("Run DNet-dist calculation", value=True)
    dist_start = start
    dist_stop = stop
    dist_step = step
    dist_max_water_distance = distance_cut_off

    custom_dist_steps = st.checkbox(
        "Set up different `start`, `stop` and `step size` parameters for the DNet-Dist module, than the global values",
        value=False,
    )

    if custom_dist_steps:
        c1, c2, c3 = st.columns([1, 1, 1])

        with c1:
            dist_start = c1.number_input(
                "Dist Start",
                value=None,
                step=1,
                help="Starting frame index for trajectory analysis. If not provided, starts from the first frame.",
                key="pKa_start",
            )

        with c2:
            dist_stop = c2.number_input(
                "Dist Stop",
                value=None,
                step=1,
                help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -2000 reads the last 2000 frames of the trajectory.",
                key="pKa_stop",
            )

        with c3:
            dist_step = c3.number_input(
                "Dist Step",
                value=None,
                step=1,
                help="Step size for reading the trajectory frames. For example if its set to 10, only  every 10th frame is read, which reduces the computation time.",
                key="pKa_steo",
            )

    custom_max_water_distance = st.checkbox(
        f"Search water molecules withing a different distance cut off than the H-bond search of `{distance_cut_off} Å`",
        value=False,
    )
    if custom_max_water_distance:
        dist_max_water_distance = st.number_input(
            "Maximum distance (in Å) from the side chain atoms within water molecules are considered for analysis",
            value=distance_cut_off,
            step=0.1,
            min_value=0.1,
            max_value=15.0,
        )

with t3:
    plot_run = st.checkbox("Run DNet-plot module", value=True)
    if not run_dist:
        st.error(
            "DNet-Plot module relies on the DNet-Dist module. Please enable the DNet-Dist calculation."
        )
        st.stop()
    plot_frame_to_time = st.number_input(
        "Convert frame to time. X frames corresponds to 1ns. X=",
        value=100,
        step=1,
        min_value=1,
        max_value=10000,
        help="Frame to time conversion. How many frames corresponds to 1ns. 100 if coordinates are written every 10 ps",
    )
    pmf_last_nth_frames = st.number_input(
        "Calculate the PMF profile from the last X number of frames X=",
        value=20000,
        step=1,
        min_value=1,
        max_value=10000000000,
        help="By default the PMF is calculated for the last 200ns of the trajectory. If the coordinates are written every 10ps it is the last 20000 frames.",
    )


with t4:
    if "graph_sets" not in st.session_state:
        st.session_state.graph_sets = []

    col_add, col_reset = st.columns([1, 1])

    def add_graph():
        if len(st.session_state.graph_sets) < 6:
            unique_id = str(time.time())
            st.session_state.graph_sets.append(
                {
                    "uid": unique_id,
                    "calcualtion_name": None,
                    "nodes_colored_by": None,
                    "edges_colored_by": None,
                    "component_search": None,
                    "distance_cut_off": None,
                    "angle_cut_off": None,
                    "min_occupancy": None,
                    "max_water": None,
                    "selection": None,
                }
            )
        else:
            st.warning("Maximum 6 additional calculations allowed.")

    def remove_graph(uid):
        st.session_state.graph_ids.remove(uid)
        for key in list(st.session_state.keys()):
            if uid in key:
                del st.session_state[key]

    if col_add.button("➕ Additional H-bond Calculation"):
        add_graph()

    for i, g_set in enumerate(st.session_state.graph_sets):

        uid = g_set["uid"]
        exp_col, btn_col = st.columns([8, 1])

        with exp_col:
            with st.expander(f"Calculation Set #{i+1}", expanded=True):

                calcualtion_name = st.text_input(
                    "Give a unique name to the calculation set:",
                    value=f"H-bond_graph_calc_{i + 1}",
                    help="The tool will create a folder with the calculation name, where the results will be saved.",
                )
                calcualtion_name = re.sub(
                    r"[^a-zA-Z0-9\s\-_]", "", calcualtion_name
                ).replace(" ", "_")

                if not calcualtion_name:
                    calcualtion_name = f"H-bond_graph_calc_{i + 1}"

                st.write(
                    f"Results will be saved in the `{output_folder}/{calcualtion_name}` folder"
                )

                nodes_colored_by = st.selectbox(
                    f"Nodes colored by (Set {i+1}):",
                    options=[
                        None,
                        "Most frequently sampled pKa value",
                        "Avg. number of water molecules around the amino acid side chain",
                    ],
                    key=f"node_col_{uid}",
                )

                edges_colored_by = st.selectbox(
                    f"Edges colored by (Set {i+1}):",
                    options=[
                        None,
                        "Occupancy",
                        "PN: population number (estimated number of conformations sampled by the connection)",
                    ],
                    key=f"edge_col_{uid}",
                )

                is_atomwise = (
                    edges_colored_by
                    == "PN: population number (estimated number of conformations sampled by the connection)"
                )

                if (
                    nodes_colored_by == "Most frequently sampled pKa value"
                    and is_atomwise
                ):
                    st.error(
                        """Coloring nodes by **Most frequently sampled pKa value** and coloring edges by **PN: population number** can't be performed in the same calculation set.
                        \nPlease set up separate calculation sets for each coloring types."""
                    )
                    st.stop()

                if (
                    nodes_colored_by
                    == "Avg. number of water molecules around the amino acid side chain"
                    and is_atomwise
                ):
                    st.error(
                        """Coloring nodes by **Avg. number of water molecules around the amino acid side chain** and coloring edges by **PN: population number** can't be performed in the same calculation set.
                        \nPlease set up separate calculation sets for each coloring types."""
                    )
                    st.stop()

                component_search = st.selectbox(
                    f"Component Search Mode (Set {i+1}):",
                    options=[
                        None,
                        "Connected component of a root node",
                        "Path search between start and goal nodes",
                    ],
                    key=f"comp_search_{uid}",
                )

                root_node = None
                start_node = None
                goal_node = None

                if component_search == "Connected component of a root node":
                    root_node = st.text_input(
                        f"Root Node (Set {i+1})",
                        placeholder=(
                            "e.g. PROA-ASP-32-OD1" if is_atomwise else "PROA-ASP-32"
                        ),
                        key=f"root_{uid}",
                        help="Format: segname-resname-resid.",
                    )
                    if root_node and len(root_node.split("-")) == 3 and is_atomwise:
                        st.error(
                            "Coloring edges by **PN: population number** requires and atomwise graph calculation. Please give the root node name in atomwise format, e.g: PROA-ASP-32-OD1 or change the edge coloring method."
                        )
                        st.stop()

                elif component_search == "Path search between start and goal nodes":
                    start_node = st.text_input(
                        f"Start Node (Set {i+1})",
                        placeholder=(
                            "e.g. PROA-ASP-32-OD1" if is_atomwise else "PROA-ASP-32"
                        ),
                        key=f"start_{uid}",
                    )
                    if start_node and len(start_node.split("-")) == 3 and is_atomwise:
                        st.error(
                            "Coloring edges by **PN: population number** requires and atomwise graph calculation. Please give the start node name in atomwise format, e.g: PROA-ASP-32-OD1 or change the edge coloring method."
                        )
                        st.stop()
                    goal_node = st.text_input(
                        f"Goal Node (Set {i+1})",
                        placeholder=(
                            "e.g. PROA-SER-22-OG" if is_atomwise else "PROA-SER-22"
                        ),
                        key=f"goal_{uid}",
                    )
                    if goal_node and len(goal_node.split("-")) == 3 and is_atomwise:
                        st.error(
                            "Coloring edges by **PN: population number** requires and atomwise graph calculation. Please give the goal node name in atomwise format, e.g: PROA-SER-22-OG or change the edge coloring method."
                        )
                        st.stop()

                    if (
                        start_node
                        and goal_node
                        and start_node.strip() == goal_node.strip()
                    ):
                        st.error(f"Start and Goal nodes are the same in Set {i+1}!")
                        st.stop()

                st.markdown(
                    f"""
                **The global H-bond criteria is set as:**
                * Distance cut off: `{distance_cut_off}`
                * Angle cut off: `{angle_cut_off}`
                * Min H-bond occupancy: `{occupancy}`
                * Max water: `{max_water}`
                """
                )
                c_distance_cut_off = distance_cut_off
                c_angle_cut_off = angle_cut_off
                c_occupancy = occupancy
                c_max_water = max_water
                c_selection = selection

                custom_criteria = st.checkbox(
                    "Use different H-bond criteria for this graph",
                    value=False,
                    key=f"crit_bool_{uid}",
                )
                if custom_criteria:
                    c_distance_cut_off, c_angle_cut_off, c_occupancy, c_max_water = (
                        H_bond_critera_element(i + 1)
                    )
                    c_selection = st.text_input(
                        "Selection to perform the analysis on, in a form of an MDAnalysis selection sting. ",
                        "protein",
                        help="Default is protein. Documentation of the selection language: `https://userguide.mdanalysis.org/stable/selections.html`",
                        key=f"selection_{uid}",
                    )

                st.session_state.graph_sets[i] = {
                    "uid": uid,
                    "calcualtion_name": calcualtion_name,
                    "nodes_colored_by": nodes_colored_by,
                    "edges_colored_by": edges_colored_by,
                    "component_search": component_search,
                    "root_node": root_node,
                    "start_node": start_node,
                    "goal_node": goal_node,
                    "distance_cut_off": c_distance_cut_off,
                    "angle_cut_off": c_angle_cut_off,
                    "min_occupancy": c_occupancy,
                    "max_water": c_max_water,
                    "selection": c_selection,
                }

        if i == len(st.session_state.graph_sets) - 1:
            with btn_col:
                if st.button(
                    "❌", key=f"remove_{uid}", help="Remove the last calculation set"
                ):
                    st.session_state.graph_sets.pop(i)
                    st.rerun()


# --- 3. THE PREVIEW (Updates instantly) ---
st.divider()
st.subheader("📄 Live Bash Script Preview")

env_params = {
    "venv_path": venv,
    "dnet_dir": dnet_dir,
    "output_folder": output_folder,
    "dnet_run_name": dnet_run_name,
    "psf": psf,
    "dcd": dcd,
    "wrap_dcd": wrap_dcd,
}

global_params = {
    "selection": selection,
    "start": start,
    "stop": stop,
    "step": step,
    "distance_cut_off": distance_cut_off,
    "angle_cut_off": angle_cut_off,
    "occupancy": occupancy,
    "max_water": max_water,
    "donors": donors,
    "acceptors": acceptors,
    "no_label_plots": no_label_plots,
    "dont_save_graph_objects": dont_save_graph_objects,
    "collect_angles": collect_angles,
    "res_id_label_shift": res_id_label_shift,
    "plot_parameters": plot_parameters,
}

DNet_pKa_params = {
    "run_pKa": run_pKa,
    "pKa_start": pKa_start,
    "pKa_stop": pKa_stop,
    "pKa_step": pKa_step,
    "pKa_selection": pKa_selection,
}

DNet_Dist_params = {
    "run_dist": run_dist,
    "dist_start": dist_start,
    "dist_stop": dist_stop,
    "dist_step": dist_step,
    "dist_max_water_distance": dist_max_water_distance,
}

DNet_Plot_params = {
    "plot_run": plot_run,
    "plot_frame_to_time": plot_frame_to_time,
    "pmf_last_nth_frames": pmf_last_nth_frames,
}

dnet_run_script = generate_bash(
    env_params,
    global_params,
    DNet_pKa_params,
    DNet_Dist_params,
    DNet_Plot_params,
    st.session_state.graph_sets,
)

st.code(dnet_run_script, language="bash")

st.warning(
    f"After downloading the run file, make the file executable with: `chmod +x run_{dnet_run_name}` "
)
st.download_button(
    label=f"Download DNet Run File",
    data=dnet_run_script,
    file_name=f"run_{dnet_run_name}",
    mime="text/x-shellscript",
)
