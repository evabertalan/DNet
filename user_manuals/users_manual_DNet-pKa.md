# Description
The `DNetPKa` script is designed to analyze and compute pKa values from molecular dynamics (MD) trajectories. It utilizes MDAnalysis to process Protein Structure Files (PSF) and trajectory files (DCD), applying pKa calculations through the `PropkaTraj` module. The script supports various functionalities including residue selection, statistical analysis, visualization of pKa values over time, and exporting data for further analysis.

# Features
- Processes MD trajectory files (PSF and DCD) to compute pKa values.
- Allows selection of specific residues using MDAnalysis selection syntax.
- Computes statistical summaries of pKa values.
- Generates plots for time-series and statistical distributions of pKa values.
- Exports computed pKa values and statistics to external files.
- Supports optional filtering of residues using a C-Graphs `_info.txt` file.

# User Manual

## Prerequisites
Ensure the following Python packages are installed:
- `MDAnalysis`
- `seaborn`
- `matplotlib`
- `pandas`
- `propkatraj`

## Usage
Run the script with the required input files and optional arguments:

### Basic Execution
```sh
python script.py <psf_file> <dcd_files>
```
Example:
```sh
python script.py protein.psf traj_1.dcd traj_2.dcd
```

### Command-line Arguments
- `psf` (required): Path to the Protein Structure File (PSF).
- `dcd` (required): Path(s) to the trajectory DCD files. Can be multiple files separated by spaces or a wildcard (e.g., `traj_*.dcd`).
- `--selection` (optional): Atom selection string (default: "protein").
- `--start` (optional): Starting frame index for trajectory analysis.
- `--stop` (optional): Stopping frame index.
- `--step` (optional): Step size for iterating frames.
- `--output_folder` (optional): Directory for saving output files (default: same directory as PSF file).
- `--plot` (optional): If set, generates and saves pKa plots.
- `--cgraphs_input` (optional): Path to a C-Graphs `_info.txt` file to filter residues.

### Output Files
Upon execution, the script generates the following files in the specified `--output_folder`:
- `pkas_for_frames_<base_name>.csv`: pKa values for each frame.
- `pkas_stats_<base_name>.csv`: Statistical summary of pKa values.
- `<base_name>_data.txt`: Most frequent pKa values for external visualization.
- `ts_<selection>_<base_name>.png`: Time-series plot (if `--plot` is enabled).
- `stats_<selection>_<base_name>.png`: Statistical distribution plot (if `--plot` is enabled).

### Example Usage with Selection and Output Folder
```sh
python script.py protein.psf traj_*.dcd --selection "resname ASP or resname GLU" --output_folder results --plot
```

This command computes pKa values for ASP and GLU residues in all `traj_*.dcd` files, saving results in the `results/` directory with plots enabled.

## Error Handling
- If no valid DCD files are found, the script raises a `FileNotFoundError`.
- If an invalid PSF file is provided, it raises a `ValueError`.
- If no pKa data is computed before analysis, the script warns the user to run `compute_pka_for_traj` first.

## Additional Functionalities
- `print_selection_info()`: Prints number of atoms, residues, and frames in the selection.
- `get_pka_for_frame(write_to_file)`: Saves pKa values per frame to a CSV file.
- `get_pka_statistic(selection, write_to_file)`: Computes and saves statistical summaries.
- `plot_pka_time_series_for_selection(selection, write_to_file, figsize)`: Generates a time-series plot.
- `plot_pka_statistic_for_selection(selection, write_to_file, figsize)`: Creates a statistical box plot.
- `write_pka_to_external_data_file(write_to_file)`: Saves most frequent pKa values for external visualization.

# Conclusion
This script provides a comprehensive tool for analyzing pKa changes in MD simulations, offering both computational and visualization capabilities. By leveraging its command-line arguments and output options, users can efficiently process trajectory data and gain insights into residue-specific pKa fluctuations.

