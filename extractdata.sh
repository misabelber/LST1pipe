#! /bin/bash

HERE=`pwd`
DATA_PATH='/scratch/bernardos/LST1/Gamma/'

cd $DATA_PATH
files=(`ls`)

cd $HERE

for i in "${files[@]}"
do
    python LST1_Hillas.py Gamma $i
done
