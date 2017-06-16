######################
# Set up directories #
######################

# create a dir in home dir as follows:
mkdir ~/H12_visualization

# create a dir for scripts
mkdir ~/H12_visualization/scripts/
# git clone scripts into this directory
git clone 

# create a dir for data
mkdir ~/H12_visualization/data
# put the genotype data in this dir:
~/H12_visualization/data/genotypeFiles

# create a dir for H12 output (i.e. computations)
mkdir ~/H12_visualization/H12_output

# create a dir for figures
mkdir ~/H12_visualization/figures

##############################################
# process data into comma-delimited format   #
##############################################

# now process the data by removing any fractions and invariants:
# remove any 'fractions' from the file
for x in {1..20}; do
    dir=~/H12_visualization/data/
    mkdir ${dir}/processedData/
    mkdir ${dir}/processedData/no_fractions
    zcat ${dir}/genotypeFiles/chr${x}_genotypes.txt.gz | grep -v '/' > ${dir}/processedData/no_fractions/chr${x}_genotypes.txt 
done

# remove any invariants
for x in {1..20}; do
    dir=~/H12_visualization/data
    mkdir ${dir}/processedData/no_invariants/
    python ~/H12_visualization/scripts/removeInvariants.py ${dir}/processedData/chr${x}_genotypes.txt ${dir}/processedData/no_invariants/chr${x}_no_invariants.txt 
done

############################
# visualize the data       #
############################

#example:
chr=1
coord=2901144
win_H12=201
win_viz=300
num_samples=29

bash ~/H12_visualization/scripts/visualize_genomic_data_wrapper.sh \
    $chr \
    $coord \
    $win_H12 \
    $win_viz \
    $num_samples
