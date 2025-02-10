# **User's Manual: Hydrogen Bond Distance Analysis Script**

## **1. Overview**
This script analyzes **hydrogen bond networks** and **water molecule distribution** in **molecular dynamics (MD) simulations** using data from:
- **PSF files** (Protein Structure File)
- **DCD trajectories** (CHARMM/NAMD simulation format)
- **C-Graphs output files** containing hydrogen bond network information

The script calculates:
1. **Distances between hydrogen-bonded residues or atoms** over the simulation trajectory.
2. **Number of water molecules** near these bonded groups.
3. **Saves results in CSV files** for further analysis.

---

## **2. Prerequisites**
### **2.1. Required Software & Libraries**
Ensure you have the following installed:
- **Python 3.7+**
- **MDAnalysis** (`pip install MDAnalysis`)
- **NumPy** (`pip install numpy`)
- **Pandas** (`pip install pandas`)
- **Matplotlib** (`pip install matplotlib`)

### **2.2. Input Files**
- **PSF file (`.psf`)**: Defines molecular topology.
- **DCD file (`.dcd`)**: Molecular dynamics trajectory.
- **C-Graphs output (`_info.txt`)**: Lists hydrogen bond connections.

---

## **3. Usage Instructions**
### **3.1. Command-Line Syntax**
Run the script from the terminal:
```bash
python script.py <PSF_FILE> <DCD_FILE(S)> <CGRAPHS_FILE> --output_folder <OUTPUT_DIR>
```

### **3.2. Command-Line Arguments**
| Argument | Required | Description |
|----------|----------|-------------|
| `<PSF_FILE>` | Yes | Path to the **PSF** topology file. |
| `<DCD_FILE(S)>` | Yes | One or more **DCD** trajectory files. Wildcards (`*.dcd`) are supported. |
| `<CGRAPHS_FILE>` | Yes | Path to the **_info.txt** file from **C-Graphs** containing bond network information. |
| `--output_folder <OUTPUT_DIR>` | No | Directory for saving output files. Default: Same as PSF file location. |
| `--selection <STRING>` | No | MDAnalysis selection string to filter atoms (e.g., `"protein and name CA"`). |
| `--start <INT>` | No | First frame to analyze (default: first frame). |
| `--stop <INT>` | No | Last frame to analyze (default: last frame). |
| `--step <INT>` | No | Frame step size (e.g., `--step 10` analyzes every 10th frame). |
| `--max_water_distance <FLOAT>` | No | Distance cutoff for water molecules (default: 3.5 Å). |
| `--wrap_dcd <True/False>` | No | Wrap molecules in the periodic box (default: `True`). |

### **3.3. Example Usage**
#### **Basic Example**
```bash
python script.py protein.psf trajectory.dcd hbond_info.txt --output_folder results/
```

#### **Analyzing Every 5th Frame from Frame 100 to 1000**
```bash
python script.py protein.psf trajectory.dcd hbond_info.txt --start 100 --stop 1000 --step 5
```

#### **Analyzing Water Distribution Around Protein Backbone Only**
```bash
python script.py protein.psf trajectory.dcd hbond_info.txt --selection "protein and backbone"
```

---

## **4. Output Files & Interpretation**
The script generates the following CSV files in the output folder:

| File | Description |
|------|-------------|
| `<base_name>_pair_distances.csv` | Distances between bonded atoms over time. |
| `<base_name>_waters_within_3_5_of_group.csv` | Number of water molecules near each **hydrogen bonding group**. |
| `<base_name>_total_waters_within_3_5_of_res.csv` | Total number of water molecules around each **residue**. |

#### **Example: Distance CSV File (`protein_pair_distances.csv`)**
| Frame | ALA-12 - GLY-15 | LYS-30 - ASP-45 | ... |
|-------|---------------|---------------|-----|
| 0 | 3.2 | 4.1 | ... |
| 1 | 3.3 | 4.2 | ... |
| 2 | 3.1 | 4.0 | ... |

Each column represents a hydrogen-bonded **residue pair**, and values are distances (in Angstroms) at different trajectory frames.

---

## **5. Advanced Features & Notes**
### **5.1. Handling Multiple DCD Files**
If your simulation is split into multiple DCD files, use:
```bash
python script.py protein.psf *.dcd hbond_info.txt
```
The script will automatically merge and analyze them in order.

### **5.2. Plotting Distance Trends (Optional)**
Uncomment the `plot_results()` function in `main()` to generate a **distance vs. time plot**:
```python
# dist_traj.plot_results()
```
This will save **`dist_time_series_per_res.png`** showing distance fluctuations over time.

---

## **6. Troubleshooting**
### **6.1. Common Errors**
| Issue | Possible Cause | Solution |
|-------|---------------|-----------|
| `FileNotFoundError: No valid DCD files found.` | Incorrect file path or missing files. | Check file names and paths. Use absolute paths if needed. |
| `ValueError: The first argument has to be a PSF file.` | Incorrect first argument. | Ensure the **PSF file** is the first input. |
| `SyntaxError: unexpected EOF while parsing` | Malformed **_info.txt** file. | Verify that the **C-Graphs** output is correctly formatted. |

### **6.2. Performance Optimization**
- Use **`--step`** to reduce the number of analyzed frames.
- Ensure **selection strings** only include necessary atoms.
- Run the script on a **high-performance machine** for large trajectories.

---

## **7. Conclusion**
This script provides a **powerful tool for analyzing hydrogen bond networks** and **water distribution** in molecular dynamics simulations. By leveraging **MDAnalysis**, users can efficiently extract and visualize distance trends over time.

For further assistance, refer to the **MDAnalysis documentation**: https://userguide.mdanalysis.org/

**Happy Analyzing! 🚀**

