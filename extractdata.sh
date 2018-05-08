#! /bin/bash

HERE=`pwd`

PARTICLE='Gamma'

DATA_PATH='/scratch/bernardos/LST1/'$PARTICLE

cd $DATA_PATH
files=(`ls`)

cd $HERE

for i in "${files[@]}"
do
    python LST1_Hillas.py $PARTICLE $i
done
