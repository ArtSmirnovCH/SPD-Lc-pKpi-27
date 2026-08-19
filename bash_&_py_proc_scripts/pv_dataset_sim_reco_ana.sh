#!/bin/sh

#SBATCH --tmp=10G # $TMPDIR space
#SBATCH --output=none


CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1

N_event=15000

SEED=${SLURM_ARRAY_TASK_ID}
Signal=1

cd ${TMPDIR}

cp /afs/jinr.ru/user/a/asmirnov/sim_27.cpp ./
cp /afs/jinr.ru/user/a/asmirnov/reco_27.cpp ./
cp /afs/jinr.ru/user/a/asmirnov/pv_dataset_27.cpp ./

singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "sim_27.cpp(${N_event}, ${SEED}, ${Signal})" > /dev/null
rm core.*
singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "reco_27.cpp(${SEED})" > /dev/null
rm core.*
singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "pv_dataset_27.cpp(${SEED})" > /dev/null
rm core.*

rm run_${SEED}.root
rm params_${SEED}.root
rm reco_full_${SEED}.root

mv ana_PV_27.root /eos/user/a/asmirnov/reco_data_lambda_c/sig_data/
mv info_analysis_PV.txt /afs/jinr.ru/user/a/asmirnov/

exit 0

