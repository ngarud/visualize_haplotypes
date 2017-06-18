# Nandita Garud and Pleuni Pennings
# ngarud@stanford.edu and pennings@sfsu.edu
# Stanford University and San Francisco State University
# November 2014

dataFile=$1
analysisWindow=$2
coord=$3
window=$4
sampleSize=$5
outFile=$6

lineNo=`cat $dataFile | cut -f1 -d',' | grep -w $coord -n | cut -f1 -d':'`
halfWindow=`echo $(( window/2 ))`
upperLineNo=`echo $((lineNo+ halfWindow))`
doubleHalf=`echo $(( halfWindow*2 +1))`

# make a dir for tmp intermediate files
mkdir ~/H12_visualization/tmp_intermediate

cat $dataFile | head -$upperLineNo | tail -$doubleHalf > ~/H12_visualization/tmp_intermediate/${coord}_data_tmp.txt

lineNo2=`cat $analysisWindow | cut -f1 | grep -w $coord -n | cut -f1 -d':'`
cat $analysisWindow | head -$lineNo2 | tail -1 > ~/H12_visualization/tmp_intermediate/${coord}_H12_H2H1_tmp.txt

script_dir=~/H12_visualization/scripts/
Rscript ${script_dir}/hapData_viz_annotate_strain.R ~/H12_visualization/tmp_intermediate/${coord}_H12_H2H1_tmp.txt ~/H12_visualization/tmp_intermediate/${coord}_data_tmp.txt $outFile $window $sampleSize

rm ~/H12_visualization/tmp_intermediate/${coord}_data_tmp.txt 
rm ~/H12_visualization/tmp_intermediate/${coord}_H12_H2H1_tmp.txt
