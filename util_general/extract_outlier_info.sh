#!/bin/bash

logpath="/Volumes/jcsf_data/data/downscaled_stations/ncamerica";

fname_part=$1
outlier_type=$2
outname=$3


fnames=$(ls $logpath/*${fname_part}*.log)

#echo "fnames:"
#printf "%s\n" $fnames
#echo

nfiles=$(echo $fnames | wc -w)
echo $nfiles files

grep -h "alt_mapval" $fnames | head -1 > $outname

#echo "grep $outlier_type $fnames >> $outname"

grep -h $outlier_type $fnames >> $outname

echo "done. outliers written to $outname"
echo


