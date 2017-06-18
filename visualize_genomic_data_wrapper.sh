# Nandita Garud and Pleuni Pennings
# ngarud@stanford.edu and pennings@sfsu.edu
# Stanford University and San Francisco State University
# November 2014
chr=$1
coord=$2
win_H12=$3
win_viz=$4
num_samples=$5

dir=~/H12_visualization

# run H12
echo 'Running H12 for' chr${chr} ${coord} '...'

python ${dir}/scripts/H12_H2H1.py \
    ${dir}/data/processedData/no_invariants/chr${chr}_no_invariants.txt \
    $num_samples \
    -o ~/H12_visualization/H12_output/chr${chr}_${coord}_H12_w${win_H12}.txt \
    -w ${win_H12} \
    -j 1 \
    -d 0 \
    -s $coord 

echo 'Done computing H12!'

# run visualization script
echo 'Visualizing...'

bash ~/westway-home/rats/documentation_rats/SelectionHapStats/scripts/visualizeGenomicData.sh \
    ${dir}/data/processedData/no_invariants/chr${chr}_no_invariants.txt \
    ${dir}/H12_output/chr${chr}_${coord}_H12_w${win_H12}.txt \
    $coord \
    $win_viz \
    $num_samples \
    ${dir}/figures/genomic_visualization_chr${chr}_coord${coord}_winH12_${win_H12}_winViz_${win_viz}.pdf

echo 'Done visualizing!'
