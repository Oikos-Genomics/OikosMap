#!/bin/bash
# Keiler Collier
# Script designed to test logic of pipeline without bioinformatic processing
# Uses nextflow's stub functions

####################################################################################
### SOURCING CONDA
eval "$(conda shell.bash hook)"


####################################################################################
### MAIN
# Run program
conda activate nextflow

nextflow run OikosMap.nf --indir ./testdata/readfiles --refseq ./testdata/refseq/Chlamydotis_macqueenii_Ctg1_15kbp.fa --prefix test -with-trace -with-report
