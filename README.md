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


### 1. DNet-pKa
The dnet_pKa module is to analyze and compute pKa values from molecular dynamics (MD) trajectories. It utilizes [MDAnalysis](https://www.mdanalysis.org/) to process Protein Structure Files (PSF) and trajectory files (DCD), applying pKa calculations through the [PropkaTraj](https://github.com/Becksteinlab/propkatraj). 
* Processes MD trajectory files (PSF and DCD) to compute pKa values.
* Allows selection of specific residues using MDAnalysis selection syntax.
* Computes statistical summaries of pKa values.
* Generates plots for time-series and statistical distributions of pKa values.
* Exports computed pKa values and statistics to external files.
* Supports optional filtering of residues using a [C-Graphs](https://github.com/evabertalan/cgraphs)  `_info.txt ` file.


#### Command-line Arguments
positional arguments:
  psf                   Path to the Protein Structure File (PSF). This file defines the molecular topology of the system.
  dcd                   Path(s) to the trajectory DCD files. You can provide multiple DCD files separated by spaces or use a wildcard pattern (e.g., 'traj_*.dcd'). These files
                        contain the molecular dynamics trajectory data.

optional arguments:
  -h, --help            show this help message and exit
  --selection SELECTION
                        Atom selection for pKa calculation using MDAnalysis syntax. Default is 'protein'. You can use selections like 'resname ASP or resname GLU' to focus on
                        specific residues.
  --start START         Starting frame index for trajectory analysis. If not provided, starts from the first frame.
  --stop STOP           Stopping frame index for trajectory analysis. If not provided, processes until the last frame.
  --step STEP           Step size for iterating through the trajectory frames. For example, '--step 10' processes every 10th frame to reduce computation time.
  --output_folder OUTPUT_FOLDER
                        Directory where output files (pKa data, statistics, and plots) will be saved. If not specified, defaults to the directory of the PSF file.
  --plot                If enabled, generates and saves pKa time series and statistical plots.
  --cgraphs_input CGRAPHS_INPUT
                        Path to a C-Graphs `_info.txt` file containing precomputed residue connectivity information. If provided, pKa values will be computed only for residues
                        found in this file, which can be further restricted using the '--selection' argument.

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

* pkas_for_frames_<base_name>.csv: pKa values for each frame
* pkas_stats_<base_name>.csv: Statistical summary of pKa values
* <base_name>_data.txt: Most frequent pKa values for external visualization.
* ts_<selection>_<base_name>.png: Time-series plot (if --plot is enabled).
* stats_<selection>_<base_name>.png: Statistical distribution plot (if --plot is enabled).

Where <base_name> is the base file name of the PSF file.
