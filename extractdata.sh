#! /bin/bash

HERE=`pwd`

PARTICLE='Gamma'

DATA_PATH='/scratch/bernardos/LST1/'$PARTICLE

cd $DATA_PATH
files=(`ls *.gz`)

cd $HERE

for i in "${files[@]}"
do
    python  /afs/ciemat.es/user/b/bernardos/GitHub/cta-lstchain/reco/LST1_Hillas.py --filename=$DATA_PATH'/'$i --outdir=/scratch/bernardos/LST1/Events --filetype=hdf5
    rm $DATA_PATH'/'$i
done


