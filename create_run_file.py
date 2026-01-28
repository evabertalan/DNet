import streamlit as st


def generate_bash(env, pka, graph, dist, plot):
    # Header & Environment
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
LOGFILE="$OUTPUT_FOLDER/dnet_output.log"

echo "Starting DNet Workflow: $(date)" > "$LOGFILE"

# Global Files
PSF="{env['psf']}"


$(ls /Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/9cis_m103a/9cis_m103a*.dcd | grep -v 'PBC.dcd')
DCD=$(ls {env['dcd']})

psf_name_no_ext=$(basename "$PSF_FILE" | cut -d. -f1)

cd {env['dnet_dir']} || exit

"""

    # dnet_pka
    if pka["run_pKa"]:
        bash_script += "\n# --- Module: dnet_pka ---\n"
        out = pka["out"] if pka["out"] else "$OUTPUT_FOLDER/pKa"
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_pka "$PSF" $DCD --output_folder "{out}" --step {pka["step"]} --selection "{pka["sel"]}"'
        cmd += "cp $PKA_FOLDER/${psf_name_no_ext}_direct_data.txt $RESIDUEWISE_FOLDER"

        bash_script += f'{cmd} >> "$LOGFILE" 2>&1\n'

    # dnet_graphs
    if graph["run"]:
        bash_script += "\n# --- Module: dnet_graphs ---\n"
        out = graph["out"] if graph["out"] else "$OUTPUT_FOLDER/graphs"
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_graphs "$PSF" $DCD --output_folder "{out}" --max_water {graph["max_w"]} --occupancy {graph["occ"]} --step {graph["step"]}'
        if graph["atomwise"]:
            cmd += " --atomwise"
        if graph["coll_ang"]:
            cmd += " --collect_angles"
        bash_script += f'{cmd} >> "$LOGFILE" 2>&1\n'

    # dnet_dist
    if dist["run"]:
        bash_script += "\n# --- Module: dnet_dist ---\n"
        out = dist["out"] if dist["out"] else "$OUTPUT_FOLDER/distances"
        g_in = (
            dist["g_in"]
            if dist["g_in"]
            else "$OUTPUT_FOLDER/graphs/workfolder/*/*_info.txt"
        )
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_dist "$PSF" $DCD "{g_in}" --output_folder "{out}" --max_water_distance {dist["max_wd"]}'
        bash_script += f'{cmd} >> "$LOGFILE" 2>&1\n'

    # dnet_plot
    if plot["run"]:
        bash_script += "\n# --- Module: dnet_plot ---\n"
        out = plot["out"] if plot["out"] else "$OUTPUT_FOLDER/plots"
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_plot --plot_folder "{out}" --graphs_info_txt "{plot["g_in"]}" --pair_distances_csv "{plot["p_csv"]}" --water_within_csv "{plot["w_csv"]}" --total_water_within_csv "{plot["tw_csv"]}" --step {plot["step"]}'
        bash_script += f'\n{cmd} >> "$LOGFILE" 2>&1\n'

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
    st.caption("*No DCD files are selected.*")


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

wrap_dcc = st.checkbox(
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
        help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -2000 reads the last 2000 frames of the trajectory.",
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
    st.stop()

st.divider()

st.header("Module specific parameters")
st.subheader("Global H-bond network criteria")
st.info(
    "These parameters set up the global H-bond criteria for the H-bond calculation in the `DNet-dist` and `DNet-graphs` modules. Additional graph calculations with different criteria can be added under the `DNet-graphs` tab."
)

c1, c2, c3, c4 = st.columns([1, 1, 1, 1])

with c1:
    distance_off = c1.number_input(
        "Distance cut off",
        value=3.5,
        step=0.1,
        min_value=0.1,
        max_value=15.0,
        help="The distance criterion for the  H-bond search, measured between the heavy atoms. The default value is 3.5Å.",
    )

with c2:
    angle_cut_off = c2.number_input(
        "Angle cut off",
        value=60,
        step=1,
        min_value=0,
        max_value=180,
        help="Threshold value for the angle formed by the acceptor heavy atom, the H atom, and the donor heavy atom. The default value is 60°.",
    )

with c3:
    occupany = c3.number_input(
        "Min H-bond occupancy",
        value=0.1,
        step=0.01,
        min_value=0.01,
        max_value=1.0,
        help="Minimum H-bond occupancy required to include an edge in the graph (default: 0.1, which means 10% occupancy).",
    )
with c4:
    max_water = c4.number_input(
        "Max number of waters allowed in the bridge",
        value=3,
        step=1,
        min_value=0,
        max_value=5,
        help="Maximum number of water molecules allowed in the water wire connections (default: 3). When it is set to 0, only direct H-bonds are considered.",
    )

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
    value=True,
    help="Create a csv file with the angles of all donor-acceptor pairs that are within the set H-bond distance criterion in each frame.",
)

# res_id_label_shift
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

    res_id_label_shift = {}
    lines = res_id_label_shift_input.split("\n")

    for line in lines:
        if ":" in line:
            try:
                key, val = line.split(":", 1)
                key = key.strip()
                val = int(val.strip())

            except ValueError:
                st.write(
                    f"Skipping {line}. Incorrect format. Please provide in a format of `PROA: 10`"
                )
                continue


with st.expander("Adjust Graph & Plot Visuals", expanded=True):
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("#### Sizes")
        edge_w = st.number_input("Edge Width", min_value=1, value=2, step=1)
        node_s = st.number_input("Node Size", min_value=1, value=15, step=1)
        node_lbl_s = st.number_input("Node Label Size", min_value=1, value=15, step=1)
        edge_lbl_s = st.number_input("Edge Label Size", min_value=1, value=15, step=1)

    with col2:
        st.markdown("#### Colors")
        g_color = st.color_picker("Graph Color", "#808080")
        # w_color = st.color_picker("Water Node Color", "#db5c5c")
        d_color = st.color_picker("Difference Graph Color", "#129fe6")
        np_color = st.color_picker("Color of Non-Protein Nodes", "#0000ff")

    with col3:
        st.markdown("#### Layout & Format")
        title_fs = st.number_input("Title Font Size", min_value=1, value=30)
        label_fs = st.number_input("Label Font Size", min_value=1, value=36)
        tick_fs = st.number_input("Tick Font Size", min_value=1, value=33)
        res = st.number_input("Resolution (DPI)", min_value=72, value=400)

        f_width = st.number_input("Fig Width", min_value=1, value=15)
        f_height = st.number_input("Fig Height", min_value=1, value=16)

        fmts = st.multiselect(
            "Export Formats", ["png", "eps", "pdf", "svg"], default=["png", "eps"]
        )
        show_chain = st.checkbox("Show Chain Label", value=False)

plot_parameters = {
    "edge_width": edge_w,
    "node_label_size": node_lbl_s,
    "edge_label_size": edge_lbl_s,
    "node_size": node_s,
    "graph_color": g_color,
    # "water_node_color": w_color,
    "difference_graph_color": d_color,
    "non_prot_color": np_color,
    "plot_title_fontsize": title_fs,
    "plot_label_fontsize": label_fs,
    "plot_tick_fontsize": tick_fs,
    "plot_resolution": res,
    "figsize": (f_width, f_height),
    "formats": fmts,
    "show_chain_label": show_chain,
}


t1, t2, t3, t4 = st.tabs(
    [
        "DNet-pKa",
        "DNet-dist",
        "DNet-plot",
        "DNet-graphs",
    ]
)

with t1:
    run_pKa = st.checkbox("Run DNet-pKa module", value=True)
    pKa_start = start
    pKa_stop = stop
    pKa_selection = selection

    if run_pKa:
        custom_pKa_steps = st.checkbox(
            "Set up different start, stop and step size parameters for the DNet-pKa module, than the global values",
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
    dist_selection = selection
    # d_gin = st.text_input(
    #     "Graph Info Input", help="Blank = Auto-detect from Graph Step"
    # )
    # d_maxwd = st.number_input("Max Water Distance", 3.5)
    # d_out = st.text_input("Dist Output (leave blank for default)")


with t3:
    g_run = st.checkbox("Enable Graphs", value=True)
    g_max_w = st.slider("Max Water", 1, 10, 3)
    g_occ = st.slider("Min Occupancy", 0.0, 1.0, 0.1)
    g_step = st.number_input("Graph Step", 1, 1000, 100)
    g_atom = st.checkbox("Atomwise Mode", value=True)
    g_coll = st.checkbox("Collect Angles")
    g_out = st.text_input("Graph Output (leave blank for default)")


with t4:
    pl_run = st.checkbox("Enable Plotting")
    pl_in = st.text_input("Graphs Info TXT", key="pl_in")
    pl_pcsv = st.text_input("Pair Distances CSV")
    pl_wcsv = st.text_input("Water Within CSV")
    pl_twcsv = st.text_input("Total Water CSV")
    pl_out = st.text_input("Plot Output (leave blank for default)")
    pl_step = st.number_input("Plot Step", 1, 1000, 1)

# --- 3. THE PREVIEW (Updates instantly) ---
st.divider()
st.subheader("📄 Live Bash Script Preview")

env_params = {
    "venv_path": venv,
    "dnet_dir": dnet_dir,
    "output_folder": output_folder,
    "psf": psf,
    "dcd": dcd,
    "wrap_dcc": wrap_dcc,
}
pka_params = {
    "run_pKa": run_pKa,
    "pKa_start": pKa_start,
    "pKa_stop": pKa_stop,
    "pKa_selection": pKa_selection,
}

graph_d = {
    "run": g_run,
    "max_w": g_max_w,
    "occ": g_occ,
    "step": g_step,
    "atomwise": g_atom,
    "coll_ang": g_coll,
    "out": g_out,
}
dist_d = {"run": d_run, "g_in": d_gin, "max_wd": d_maxwd, "out": d_out}
plot_d = {
    "run": pl_run,
    "out": pl_out,
    "g_in": pl_in,
    "p_csv": pl_pcsv,
    "w_csv": pl_wcsv,
    "tw_csv": pl_twcsv,
    "step": pl_step,
}

# Generate the string
final_script = generate_bash(env_params, pka_params, graph_d, dist_d, plot_d)

# Display the code
st.code(final_script, language="bash")

# Download Button
st.download_button(
    label="Download run_dnet.sh",
    data=final_script,
    file_name="run_dnet.sh",
    mime="text/x-shellscript",
)
