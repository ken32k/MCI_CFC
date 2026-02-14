#!/bin/bash
# This script runs the micapipe_cleanup command for a given subject in the MCI-CFC dataset. It takes one argument, which is the subject ID. 
# The script uses environment variables to specify the input BIDS directory and the output directory for the cleaned data. The micapipe_cleanup command is run with options to process DWI data and to perform skull stripping and bias field correction.

sub=$1
bids=${MCI_CFC_DIR}/mica_bids
out=${MCI_CFC_DIR}/mica_out
micapipe_cleanup -bids $bids -out $out -sub $sub -proc_dwi -SC