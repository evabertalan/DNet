# DNet (Dynamic Networks)

**DNet** is a tool for performing pKa calculations, generating atomwise and residuewise graphs, computing pairwise distances, calculating water interactions, and generating plots based on the analysis.

EXAMPLE final output from dummy data

## How to Set Up

1. **Download the Code**  
   Clone the repository from GitHub and copy it to the location where you wish to perform the analysis.

2. **Install Required Packages**  
   The repository includes a `requirements.txt` file and an `install_dnet` script to install all necessary dependencies with the required versions.  
   To install the dependencies, execute the following command in your terminal:
   ```bash
   ./install_dnet
   ```

## Basic Setup to Run Code from the Command Line
The easiest way to run the tool is to customize the `run_dnet` included in the repo.
Running this script performs the complete calculation. The parts marked with #ADD has to be customized in the file, such as the the paths to the PSF and DCD files.
For a more customized set up and to understand all the possible input parameters please look at the advanced usage.

   - **LOCATION_OF_LOGFILE**: Path to the folder where the logfile should be saved. Example: `LOCATION_OF_LOGFILE='/path/to/logfile'`
   - **PSF_FILE**: Path to the PSF (Protein Structure File). Example: `PSF_FILE='/path/to/psf/example.psf'`
   - **DCD_FILES**: Path(s) to the DCD (trajectory) file(s), supports wildcards. Example: `DCD_FILES=$(ls /path/to/dcdfiles/example_n*.dcd | grep -v 'PBC.dcd')`
   - **OUTPUT_FOLDER**: Directory for output files and plots. Defaults to the PSF file location if not provided. Example: `OUTPUT_FOLDER='/path/to/output'`
   - **SELECTION**: Atom selection string for graph calculation. Default: `'protein'`. Example: `SELECTION="protein or resname LYR or resname HSE"`
   - **ACCEPTORS**: Additional hydrogen bond acceptor atoms in Python list format. Example: `ACCEPTORS="['O', 'F']"`
   - **DONORS**: Additional hydrogen bond donor atoms in Python list format. Example: `DONORS="['N', 'S']"`
   - **MAX_WATER**: Maximum number of water molecules in water wire connections. Default: `3`. Example: `MAX_WATER=3`
   - **OCCUPANCY**: Minimum hydrogen bond occupancy to include an edge in the graph. Default: `0.1`. Example: `OCCUPANCY=0.1`
   - **RESID_LABEL_SHIFT**: Shift residue ID labels by a specified amount in plots. Default: `0`. Example: `RESID_LABEL_SHIFT=0`
   - **STEP**: Step size for iterating through the trajectory frames. In case of large trajectories it is recommended to only read every 10th frame for the pKa and distance calculation. Default: `1`. Example: `STEP=10`



## Advanced Usage

Detailed description for each of the modules.
Each module (dnet_pKa, dnet_graphs, dnet_dist, dnet_plot) supports the built in help function: `python3 -m dnet_dist --help` which will output the list of all supported arguments for the module.


### I. DNet-pKa
The dnet_pKa module is to analyze and compute pKa values from molecular dynamics (MD) trajectories. It utilizes [MDAnalysis](https://www.mdanalysis.org/) to process Protein Structure Files (PSF) and trajectory files (DCD), applying pKa calculations through the [PropkaTraj](https://github.com/Becksteinlab/propkatraj). 
* Processes MD trajectory files (PSF and DCD) to compute pKa values.
* Allows selection of specific residues using MDAnalysis selection syntax.
* Computes statistical summaries of pKa values.
* Generates plots for time-series and statistical distributions of pKa values.
* Exports computed pKa values and statistics to external files.
* Supports optional filtering of residues using a [C-Graphs](https://github.com/evabertalan/cgraphs)  `_info.txt ` file.


#### Command-line Arguments
##### Positional arguments:
| Argument | Description |
|----------|-------------|
| `psf` | Path to the Protein Structure File (PSF). This file defines the molecular topology of the system. |
| `dcd` | Path(s) to the trajectory DCD files. You can provide multiple DCD files separated by spaces or use a wildcard pattern (e.g., 'traj_*.dcd'). These files
                        contain the molecular dynamics trajectory data.|

##### Optional arguments:
| Argument           | Default Value | Description                                                                                                                                        |
|--------------------|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| `--selection`      | `protein`     | Atom selection for pKa calculation using MDAnalysis syntax. Default is 'protein'.                                                                  |
| `--start`          | -             | Starting frame index for trajectory analysis. If not provided, starts from the first frame.                                                        |
| `--stop`           | -             | Stopping frame index for trajectory analysis. If not provided, processes until the last frame.                                                     |
| `--step`           | -             | Step size for iterating through the trajectory frames. E.g., '--step 10' processes every 10th frame.                                               |
| `--output_folder`  | -             | Directory where output files (pKa data, statistics, and plots) will be saved. Defaults to the directory of the PSF file if not specified.          |
| `--plot`           | `False`       | If enabled, generates and saves pKa time series and statistical plots.                                                                             |
| `--cgraphs_input`  | -             | Path to a C-Graphs '_info.txt' file with precomputed residue connectivity information. Restricts pKa computation to residues found in the file.    |


#### Execution
```sh
python3 -m dnet_pka <psf_file> <dcd_files>
```

Example:
```sh
python3 -m dnet_pka '../protein1.psf' 'traj_*.dcd' --selection "protein or resname LYR" --output_folder './results'  --step 10
```


#### Output files:
The script generates the following files in the specified --output_folder:

* `pkas_for_frames_<base_name>.csv`: pKa values for each frame
* `pkas_stats_<base_name>.csv`: Statistical summary of pKa values
* `<base_name>_data.txt`: Most frequent pKa values for external visualization.
* `ts_<selection>_<base_name>.png`: Time-series plot (if --plot is enabled).
* `stats_<selection>_<base_name>.png:` Statistical distribution plot (if --plot is enabled).

Where `<base_name>` is the base file name of the PSF file.


### II. DNet-Graphs
The dnet-graphs module analyzes molecular dynamics (MD) trajectories to construct and visualize water wire networks based on hydrogen bond connectivity. It allows users to process MD simulation files, extract hydrogen bond networks, and generate graphs and plots.
* Processes MD trajectory files (PSF and DCD) to compute water meditated H-bond networks with [Bridge](https://github.com/maltesie/bridge) and [C-Graphs](https://github.com/evabertalan/cgraphs).
* Customize node and edge attributes for visualization
* Generate plots for water wire graphs

#### Command-line Arguments

##### Positional arguments:
| Argument | Description |
|----------|-------------|
| `psf` | Path to the **PSF (Protein Structure File)**. |
| `dcd` | Path(s) to the **DCD (trajectory) files**. Supports wildcards (`*.dcd`). |

##### Optional arguments:
| Argument                | Default Value  | Description                                                                                                                                                    |
|-------------------------|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--output_folder`        | -              | Directory where output files and plots will be saved. If not provided, defaults to the location of the PSF file.                                              |
| `--max_water`            | `3`            | Maximum number of water molecules allowed in water wire connections.                                                                                          |
| `--occupancy`            | `0.1`          | Minimum hydrogen bond occupancy required to include an edge in the graph.                                                                                     |
| `--selection`            | `protein`      | Atom selection string for defining the region of interest in the molecular system for graph calculation.                                                      |
| `--additional_donors`    | `[]`           | List of additional hydrogen bond donor atoms, formatted as a Python list (e.g., "['N', 'S']").                                                                 |
| `--additional_acceptors` | `[]`           | List of additional hydrogen bond acceptor atoms, formatted as a Python list (e.g., "['O', 'F']").                                                             |
| `--start`                | -              | Starting frame index for trajectory analysis. If not provided, starts from the first frame.                                                                   |
| `--stop`                 | -              | Stopping frame index for trajectory analysis. If not provided, processes until the last frame.                                                                |
| `--step`                 | `1`            | Step size for iterating through the trajectory frames. For example, '--step 10' processes every 10th frame to reduce computation time.                        |
| `--residuewise`          | `True`         | Calculate hydrogen bonds at the residue level instead of the atomic level.                                                                                   |
| `--atomewise`            | `False`        | Calculate hydrogen bonds at the atomic level instead of the residue level (overrides `--residuewise`).                                                       |
| `--wrap_dcd`             | `true`         | Apply periodic boundary condition wrapping to keep molecules inside the simulation box. Use 'true' or 'false'.                                                |
| `--res_id_label_shift`   | `0`            | Shift residue ID labels by a specified amount in plots.                                                                                                       |
| `--color_data`           | `False`        | Color nodes in the graph based on external data values.                                                                                                       |
| `--node_color_selection` | `protein`      | Selection criteria for which nodes to color in the graph.                                                                                                     |
| `--node_color_map`       | `coolwarm_r`   | Colormap used for node coloring.                                                                                                                              |
| `--plot_parameters`      | `{}`           | Dictionary of plot parameters formatted as a string (e.g., `{'graph_color': '#666666', 'formats': ['png', 'eps']}`).                                           |
| `--include_backbone`     | `False`        | Include interactions between backbone and sidechain atoms in the analysis.                                                                                   |


#### Execution
```sh
python3 -m dnet_graphs <psf_file> <dcd_files>
```

Example:
```sh
python3 -m dnet_graphs '../protein1.psf' 'traj_*.dcd' --selection "protein or resname LYR" --output_folder './results'  --plot_parameters '{"graph_color": "#0000FF", "formats": ["png", "svg"]}' --wrap_dcd false --atomwise
```


#### Output files:
The script creates a `workfolder` in the location of the `--output_folder` with the following subfolders and files:
* `workfolder/graph_objects/<max_water>_water_wires/<base_name>_water_wires_graph.pickle`:  Stores the computed water wire graph.
* `workfolder/<base_name>/<base_name>_max_<max_water>_water_bridges_min_occupancy_<occupancy>_graph_labeled.png`: Water wire graph plot with nodes and edges.
* `workfolder/<base_name>/<base_name>_max_<max_water>_water_bridges_min_occupancy_<occupancy>_water_wire_graph_info.txt`: Info file with nodes and edges of the graphs listed.


Where `<base_name>` is the base file name of the PSF file.


### III. DNet-Dist
The dnet-dist module analyzes the distances between hydrogen bonding atoms of the networks and water molecule distribution around the atoms in molecular dynamics (MD) simulations.
The script calculates:
* Distances between hydrogen-bonded residues or atoms over the simulation trajectory.
* Number of water molecules near these bonded groups.
* Saves results in CSV files for further analysis.

#### Command-line Arguments
##### Positional arguments:
| Argument | Description |
|----------|-------------|
| `psf` | Path to the **PSF (Protein Structure File)**. |
| `dcd` | Path(s) to the **DCD (trajectory) files**. Supports wildcards (`*.dcd`). |
| `graphs_input` | Path to the _info.txt DNet-Grpahs file containing graph edges and nodes for distance calculations. |


##### Optional arguments:
| Argument                | Default Value  | Description                                                                                                                                                    |
|-------------------------|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--output_folder`        | -              | Path to the folder where output files (CSV and plots) will be saved.                                                                                         |
| `--selection`            | -              | MDAnalysis selection string to restrict distance calculations to a subset of atoms from the cgraphs input.                                                     |
| `--start`                | -              | Starting frame index for trajectory analysis. If not provided, starts from the first frame.                                                                   |
| `--stop`                 | -              | Stopping frame index for trajectory analysis. If not provided, processes until the last frame.                                                                |
| `--step`                 | -              | Step size for iterating through the trajectory frames. For example, '--step 10' processes every 10th frame to reduce computation time.                        |
| `--max_water_distance`   | `3.5`          | Maximum distance (in Å) within which water molecules are considered for analysis. Default is 3.5 Å.                                                          |
| `--wrap_dcd`             | `True`         | Apply wrapping transformations to keep molecules within the simulation box.                                                                                   |


#### Execution
```sh
python3 -m dnet_graphs <psf_file> <dcd_files>  <graphs_input>
```

Example:
```sh
python3 -m dnet_graphs '../protein1.psf' 'traj_*.dcd' '../workfolder/protein1/protein1_max_3_water_bridges_min_occupancy_0.1_water_wire_graph_info.txt' --output_folder './results'  --step 10
```


#### Output files:
`<base_name>_pair_distances.csv`   Distances between H-bonded atoms over time. Each column of the table represents a hydrogen-bonded residue pair, and values are distances (in Angstroms) at different trajectory frames.
`<base_name>_waters_within_3_5_of_group.csv`   Number of water molecules near each hydrogen H-bonding group.
`<base_name>_total_waters_within_3_5_of_res.csv`  Total number of water molecules around each H-bonding residue

### IV. DNet-Plot

#### Command-line Arguments
Here is the updated markdown table with an additional column indicating whether each argument is required or optional:

| Argument                | Default Value  | Required/Optional | Description                                                                                                                                                    |
|-------------------------|----------------|--------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--plot_folder`          | -              | Required           | Folder to save the plots.                                                                                                                                     |
| `--graphs_info_txt`      | -              | Required           | Path to the graphs information text file.                                                                                                                     |
| `--pKas_for_frame_csv`   | -              | Required           | CSV file for pKa values for each frame.                                                                                                                      |
| `--pair_distances_csv`   | -              | Required           | CSV file for pair distances.                                                                                                                                  |
| `--water_within_csv`     | -              | Required           | CSV file for water molecules within a certain distance.                                                                                                      |
| `--total_water_within_csv`| -             | Required           | CSV file for total water molecules within a certain distance.                                                                                                |
| `--sim_name`             | -              | Optional           | Optional argument for the simulation name.                                                                                                                   |
| `--step`                 | `1`            | Optional           | Step size for iterating through frames (default is 1).                                                                                                        |
| `--frame_to_time`        | `100`          | Optional           | Frame to time conversion factor, used for converting frame indices to time. Default is 100.                                                                  |
| `--pmf_last_nth_frames`  | `20000`        | Optional           | The number of frames to consider for PMF calculation. Default is 20000.                                                                                      |
| `--plot_formats`         | `['png']`      | Optional           | List of plot formats to save (default is `['png']`).                                                                                                          |
| `--res_id_label_shift`   | `0`            | Optional           | Shift residue ID labels by a specified amount in plots (default is 0).                                                                                        |

#### Execution

#### Output files:


### Additional Notes:
* This script is designed for MD simulations, particularly those involving water networks
* For large trajectories, use a higher --step value to speed up processing
* Use selection strings only include necessary atoms
* Run the script on a high-performance machine for large trajectories

