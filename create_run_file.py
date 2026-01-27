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
        value=None,
        step=0.1,
        min_value=0.1,
        max_value=15,
        help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -2000 reads the last 2000 frames of the trajectory.",
    )

with c2:
    angle_cut_off = c2.number_input(
        "Angle cut off",
        value=60,
        step=1,
        min_value=0,
        max_value=180,
        help="Last frame index for trajectory analysis. If not provided, processes until the last frame. It can take a negative value, e.g: -2000 reads the last 2000 frames of the trajectory.",
    )

with c3:
    occupany = c3.number_input(
        "Min H-bond occupancy",
        0.1,
        step=0.01,
        min_value=0,
        max_value=1,
        help="Step size for reading the trajectory frames. For example if its set to 10, only  every 10th frame is read, which reduces the computation time.",
    )
with c4:
    max_water = c1.number_input(
        "Max number of waters allowed in the bridge",
        value=3,
        step=1,
        min_value=0,
        max_value=5,
        help="Maximum number of water molecules allowed in the water wire connections (default: 3). When it is set to 0, only direct H-bonds are considered.",
    )


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
