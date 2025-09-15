!/bin/bash
# Keiler Collier
# Script designed to test logic of pipeline without bioinformatic processing
# Uses nextflow's stub functions

####################################################################################
### SOURCING CONDA
eval "$(conda shell.bash hook)"


####################################################################################
### MAIN
# Run program
nextflow run OikosMap.nf -stub-run --indlist ./testdata/indlist.txt --indir ./testdata/readfiles --refseq /Users/keilercollier/Documents/Repos/OikosMap/testdata/refseq/Chlamydotis_macqueenii_Ctg1_15kbp.fa --prefix test #-with-trace -with-report
