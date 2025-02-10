Descritopn: DNet (Dynamic Networks) is a tool to performs pKa calculations, generate atomwise and residuewise graphs, computes pair distances, calculates water interactions, and generates plots based on the analysis.

## How to set up:
1. download the code from GitHub and copy to the location where you would like to perform the analysis
2. install the required packages. The repo contains a requirements.txt and an install_dnet script to install all dependencies with the requred versions
   execute: `./install_dnet` in the terminal

## Basic set up to run code from command line
 Customize the run_dnet file. The repo contaisn a template run file, but the paths to the PSF and DCD files are missing. Before runninig the tool, all parts marked with `# ADD` has to be fileed out
* LOCATION_OF_LOGFILE= #Add relative path to the folder, where the logfile should be saved
* PSF_FILE='/p/project/test/dcdfiles/example.psf'
  Path to the PSF (Protein Structure File) used for molecular dynamics simulations.
* DCD_FILES=$(ls /p/project/test/dcdfiles/example_n*.dcd | grep -v 'PBC.dcd')
  Path(s) to the DCD (trajectory) file(s). Supports wildcard patterns (e.g., '*.dcd'). The expression above selects all dcd files which match the pattern, but exculed the oves ending for PBC.dcd.
* OUTPUT_FOLDER=/p/home/users/me/dnet_outputs
  Directory where output files and plots will be saved. If not provided, defaults to the location of the PSF file.

* SELECTION= # e.g: "protein or resname LYR or resname HSE"
  Atom selection string for defining the region of interest in the molecular system for graph calculation (default: 'protein').
* ACCEPTORS="[]" # e.g: "['NH']"
  List of additional hydrogen bond acceptor atoms, formatted as a Python list (e.g., "['O', 'F']").
* DONORS="[]" # e.g: "['NH']" 
  List of additional hydrogen bond donor atoms, formatted as a Python list (e.g., "['N', 'S']").
* MAX_WATER=3
   Maximum number of water molecules allowed in water wire connections (default: 3).
* OCCUPANCY=0.1
  Minimum hydrogen bond occupancy required to include an edge in the graph (default: 0.1).
* RESID_LABEL_SHIFT=0
  Shift residue ID labels by a specified amount in plots (default: 0)

## Advanced usage
Create your own run file:
(To be added, with detailed users manual for each module)
