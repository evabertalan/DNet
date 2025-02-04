# **User Manual: Water Wire Graph Analysis Script**  

## **1. Introduction**  
This script analyzes molecular dynamics (MD) trajectories to construct and visualize **water wire networks** based on hydrogen bond connectivity. It allows users to process MD simulation files, extract hydrogen bond networks, and generate informative graphs and plots.  

### **Features:**  
✅ Parse PSF and DCD files to extract molecular trajectory data  
✅ Identify hydrogen bond networks between residues and water molecules  
✅ Customize node and edge attributes for visualization  
✅ Generate plots for water wire graphs and linear lengths of components  

---

## **2. Installation and Dependencies**  
To use this script, you need Python 3 and the following dependencies:  

### **Required Python Libraries**  
Install the necessary packages using:  

```bash
pip install numpy MDAnalysis mdhbond matplotlib pathlib argparse
```

Ensure you have **MDAnalysis** installed for handling MD trajectory files.  

---

## **3. Usage Instructions**  

### **Basic Command Structure:**  
```bash
python script.py <psf_file> <dcd_files> [options]
```
Where:  
- `<psf_file>` → Path to the **PSF file** defining the molecular system  
- `<dcd_files>` → One or more **DCD trajectory files** (supports wildcards `*.dcd`)  

---

## **4. Command-Line Arguments**  

### **Required Arguments**
| Argument | Description |
|----------|-------------|
| `psf` | Path to the **PSF (Protein Structure File)**. |
| `dcd` | Path(s) to the **DCD (trajectory) files**. Supports wildcards (`*.dcd`). |

### **Optional Arguments**
| Argument | Default | Description |
|----------|---------|-------------|
| `--output_folder` | Same as PSF file's location | Directory where output files and plots are saved. |
| `--max_water` | `3` | Maximum number of water molecules allowed in a hydrogen bond bridge. |
| `--occupancy` | `0.1` | Minimum hydrogen bond occupancy required for an edge in the graph. |
| `--selection` | `"protein"` | Selection string for defining the region of interest. |
| `--additional_donors` | `"[]"` | List of extra hydrogen bond donors (e.g., `'["N", "S"]'`). |
| `--additional_acceptors` | `"[]"` | List of extra hydrogen bond acceptors (e.g., `'["O", "F"]'`). |
| `--start` | None | Starting frame index for trajectory analysis. |
| `--stop` | None | Ending frame index for trajectory analysis. |
| `--step` | `1` | Step size for iterating through trajectory frames. |
| `--residuewise` | `True` | Compute H-bonds at the **residue level**. |
| `--atomwise` | Overrides `--residuewise` | Compute H-bonds at the **atomic level** instead. |
| `--wrap_dcd` | `"true"` | Apply periodic boundary condition wrapping (`"true"` or `"false"`). |
| `--res_id_label_shift` | `0` | Adjust residue numbering in plots. |
| `--color_data` | `False` | Enable node coloring based on external data. |
| `--node_color_selection` | `"protein"` | Selection criteria for node coloring. |
| `--node_color_map` | `"coolwarm_r"` | Colormap for node visualization. |
| `--plot_parameters` | `"{}"` | Dictionary of plot parameters (e.g., `'{"graph_color": "#666666"}'`). |
| `--include_backbone` | `False` | Include backbone-sidechain interactions in calculations. |

---

## **5. Examples**  

### **Example 1: Basic Run**
```bash
python script.py system.psf trajectory.dcd
```
- Uses **default settings** to process `system.psf` and `trajectory.dcd`.  
- Outputs graphs and results in the same directory as `system.psf`.  

### **Example 2: Specifying an Output Directory**
```bash
python script.py system.psf trajectory.dcd --output_folder results/
```
- Saves results to the `results/` directory.  

### **Example 3: Analyzing Multiple DCD Files**
```bash
python script.py system.psf traj1.dcd traj2.dcd traj3.dcd
```
- Processes multiple **DCD files**.  

### **Example 4: Filtering by Occupancy**
```bash
python script.py system.psf traj.dcd --occupancy 0.3
```
- Only includes hydrogen bonds **with at least 30% occupancy**.  

### **Example 5: Including Additional Donors and Acceptors**
```bash
python script.py system.psf traj.dcd --additional_donors '["N", "S"]' --additional_acceptors '["O", "F"]'
```
- Expands hydrogen bond detection with extra donors (`N, S`) and acceptors (`O, F`).  

### **Example 6: Disabling Periodic Wrapping**
```bash
python script.py system.psf traj.dcd --wrap_dcd false
```
- Turns off **periodic boundary condition wrapping**.  

### **Example 7: Custom Plotting Parameters**
```bash
python script.py system.psf traj.dcd --plot_parameters '{"graph_color": "#0000FF", "formats": ["png", "svg"]}'
```
- Changes graph color to **blue (`#0000FF`)** and outputs **PNG & SVG formats**.  

---

## **6. Output Files**  

### **1. Processed Data Files**  
- **Graph object** (pickle file): Stores the computed water wire graph.  
- **JSON file**: Contains graph edge information (bond occupancy, distances, etc.).  

### **2. Visualization Outputs**  
- **Graph plots** (`.png`, `.svg`, etc.): Water wire graphs with nodes and edges.  
- **Linear length plots**: Shows continuous network components.  

---

## **7. Troubleshooting & FAQs**  

### **Q1: No DCD files found!**  
**Solution:**  
- Ensure correct file path and check for wildcards (`*.dcd`).  
- Run `ls` to verify file existence:  
  ```bash
  ls *.dcd
  ```

### **Q2: Unexpected `SyntaxError` when using `--additional_donors` or `--plot_parameters`?**  
**Solution:**  
- Ensure JSON-like arguments are enclosed in **single quotes `' '`**:  
  ```bash
  python script.py system.psf traj.dcd --additional_donors '["N", "S"]'
  ```

### **Q3: Plots are missing labels or colors look incorrect?**  
**Solution:**  
- Check if `--color_data` is correctly set.  
- Ensure external data files exist in the **same directory** as the PSF file.  

---

## **8. Additional Notes**  
- This script is designed for **MD simulations**, particularly those involving **water networks**.  
- For **large trajectories**, use a **higher `--step` value** to speed up processing.  

---

## **9. Credits & Acknowledgments**  
Developed using **MDAnalysis, Matplotlib, and mdhbond**.  
If you use this script for research, consider citing the relevant **MDAnalysis** and **mdhbond** papers.  

---

This manual ensures that users can **easily install, use, and troubleshoot** the script. Let me know if you need further refinements! 🚀