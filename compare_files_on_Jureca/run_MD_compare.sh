#!/bin/bash
# ========================================

source  ~/.venvs/pKa_proj_env/bin/activate

LOGFILE=cgraphs_output.log
echo "$(date '+%Y-%m-%d %H:%M:%S'): Starting C-Graphs calculation script..." >> "$LOGFILE"

# ========================================
# Under this line come the parameters to run the script with

#path to the psf file:
PSF_FILE1='/Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/9cis_m103a/read_protein_membrane_7_9cis_m103a_3_2x.psf'
PSF_FILE2='/Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/9cis_optimized/read_protein_membrane_7_opt_3_2x.psf'

#path to the dcd files. Files can be defined individually listed after each other. Or multiple files selected using the glob wildcards.
DCD_FILES1=$(ls /Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/9cis_m103a/9cis_m103a*.dcd | grep -v 'PBC.dcd')
DCD_FILES2=$(ls /Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/9cis_optimized/9cis*.dcd | grep -v 'PBC.dcd')

#path to the output folder, in this folder will be the calculation results saved
WORKFOLDER=/Users/evabertalan/Documents/projects/cgraphs/to_organize/test_cgraphs/cgraphs_test_for_script/jsr1_tests/compare_test_m103_vs_opt/

mkdir -p $WORKFOLDER

python3 -m compare_two_MD  "$PSF_FILE1"  "$PSF_FILE2" "$WORKFOLDER"  --dcd1 $DCD_FILES1 --dcd2 $DCD_FILES2 --color1 "yellow" --color2 "purple"  >> "$LOGFILE" 2>&1


# =====================================
# Don't change this part
echo "$(date '+%Y-%m-%d %H:%M:%S'): Finished C-Graphs calculation" >> "$LOGFILE"
deactivate