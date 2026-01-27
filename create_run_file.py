import streamlit as st
import json


def generate_bash(env, pka, graph, dist, plot):
    # Header & Environment
    bash_script = f"""#!/bin/bash
source {env['venv_path']}

# Setup Logging
BASE_DIR="{env['log_dir']}"
mkdir -p "$BASE_DIR"
LOGFILE="$BASE_DIR/dnet_output.log"

echo "Starting DNet Pipeline: $(date)" > "$LOGFILE"

# Global Files
PSF="{env['psf']}"
DCD="{env['dcd']}"
"""

    # dnet_pka
    if pka["run"]:
        bash_script += "\n# --- Module: dnet_pka ---\n"
        out = pka["out"] if pka["out"] else "$BASE_DIR/pKa"
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_pka "$PSF" $DCD --output_folder "{out}" --step {pka["step"]} --selection "{pka["sel"]}"'
        if pka["plot"]:
            cmd += " --plot"
        bash_script += f'{cmd} >> "$LOGFILE" 2>&1\n'

    # dnet_graphs
    if graph["run"]:
        bash_script += "\n# --- Module: dnet_graphs ---\n"
        out = graph["out"] if graph["out"] else "$BASE_DIR/graphs"
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
        out = dist["out"] if dist["out"] else "$BASE_DIR/distances"
        g_in = (
            dist["g_in"] if dist["g_in"] else "$BASE_DIR/graphs/workfolder/*/*_info.txt"
        )
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_dist "$PSF" $DCD "{g_in}" --output_folder "{out}" --max_water_distance {dist["max_wd"]}'
        bash_script += f'{cmd} >> "$LOGFILE" 2>&1\n'

    # dnet_plot
    if plot["run"]:
        bash_script += "\n# --- Module: dnet_plot ---\n"
        out = plot["out"] if plot["out"] else "$BASE_DIR/plots"
        bash_script += f'mkdir -p "{out}"\n'
        cmd = f'python3 -m dnet_plot --plot_folder "{out}" --graphs_info_txt "{plot["g_in"]}" --pair_distances_csv "{plot["p_csv"]}" --water_within_csv "{plot["w_csv"]}" --total_water_within_csv "{plot["tw_csv"]}" --step {plot["step"]}'
        bash_script += f'\n{cmd} >> "$LOGFILE" 2>&1\n'

    return bash_script


st.set_page_config(page_title="DNet Bash Script Generator", layout="wide")
st.title("DNet Bash Script Generator")


col_env, col_params = st.columns([1, 2])

with col_env:
    st.header("📂 Global Environment")
    venv = st.text_input("Venv Path", "~/.venvs/pKa_proj_env/bin/activate")
    log_dir = st.text_input("Project Root Directory", "/absolute/path/to/project")
    psf = st.text_input("PSF File", "protein.psf")
    dcd = st.text_input("DCD Files", "*.dcd")

with col_params:
    t1, t2, t3, t4 = st.tabs(["🧪 pKa", "📊 Graphs", "📏 Distances", "🎨 Plotting"])

    with t1:
        p_run = st.checkbox("Enable pKa", value=True)
        p_sel = st.text_input("pKa Selection", "protein")
        p_step = st.number_input("pKa Step", 1, 1000, 100)
        p_out = st.text_input("pKa Output (leave blank for default)")
        p_plot = st.checkbox("Generate pKa Plots", value=True)

    with t2:
        g_run = st.checkbox("Enable Graphs", value=True)
        g_max_w = st.slider("Max Water", 1, 10, 3)
        g_occ = st.slider("Min Occupancy", 0.0, 1.0, 0.1)
        g_step = st.number_input("Graph Step", 1, 1000, 100)
        g_atom = st.checkbox("Atomwise Mode", value=True)
        g_coll = st.checkbox("Collect Angles")
        g_out = st.text_input("Graph Output (leave blank for default)")

    with t3:
        d_run = st.checkbox("Enable Distances")
        d_gin = st.text_input(
            "Graph Info Input", help="Blank = Auto-detect from Graph Step"
        )
        d_maxwd = st.number_input("Max Water Distance", 3.5)
        d_out = st.text_input("Dist Output (leave blank for default)")

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

# Package inputs into dictionaries
env_d = {"venv_path": venv, "log_dir": log_dir, "psf": psf, "dcd": dcd}
pka_d = {"run": p_run, "sel": p_sel, "step": p_step, "out": p_out, "plot": p_plot}
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
final_script = generate_bash(env_d, pka_d, graph_d, dist_d, plot_d)

# Display the code
st.code(final_script, language="bash")

# Download Button
st.download_button(
    label="Download run_dnet.sh",
    data=final_script,
    file_name="run_dnet.sh",
    mime="text/x-shellscript",
)
