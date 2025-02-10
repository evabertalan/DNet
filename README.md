Here’s a revised version of your README that improves clarity, corrects typos, and formats the information for better readability:

---

# DNet (Dynamic Networks)

**DNet** is a tool for performing pKa calculations, generating atomwise and residuewise graphs, computing pairwise distances, calculating water interactions, and generating plots based on the analysis.

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

1. **Customize the `run_dnet` File**  
   The repository includes a template `run_dnet` file, but the paths to the PSF and DCD files are missing. Before running the tool, fill in the placeholders marked with `# ADD`

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


## Advanced Usage

To be added:  
Detailed user manual for each module and instructions on creating custom `run_dnet` files.

---

This format should help users follow along with setting up and running the tool. It includes clean formatting, consistent terminology, and clear instructions for customization.
