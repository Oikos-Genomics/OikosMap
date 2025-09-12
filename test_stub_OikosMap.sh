!/bin/bash
# Keiler Collier
# Script designed to test logic of pipeline without bioinformatic processing
# Uses nextflow's stub functions

####################################################################################
### SOURCING CONDA
eval "$(conda shell.bash hook)"


####################################################################################
### MAIN

mkdir -p testdata/dummy_readfiles
touch testdata/dummy_readfiles/ind1_R{1,2}.fq.gz
touch testdata/dummy_readfiles/ind2_R{1,2}.fq.gz
touch testdata/dummy_readfiles/ind3_R{1,2}.fq.gz

# create dummy refseq
mkdir -p testdata/dummy_refseq
touch testdata/dummy_refseq/refseq.fa

# Generate new --indlist
echo -e "ind1\nind2\nind3" > testdata/dummy_indlist.txt

# Run program
nextflow run OikosMap.nf -stub-run --indlist testdata/dummy_indlist.txt --indir testdata/dummy_readfiles --refseq testdata/dummy_refseq/refseq.fa --prefix dummy_out -with-trace -with-report
if [ $? -eq 0 ]; then
    echo "OikosMap ran successfully for stubs. Cleaning up."
    rm -r testdata/dummy*
    rm -r dummy*
else
    echo "OikosMap failed for subs. Check logfiles."
fi