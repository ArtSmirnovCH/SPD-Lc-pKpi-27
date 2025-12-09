#!/bin/sh

#SBATCH --tmp=10G # $TMPDIR space
#SBATCH --output=none


CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1

N_event=3000

SEED=${SLURM_ARRAY_TASK_ID}
Signal=1

cd ${TMPDIR}

cp /afs/jinr.ru/user/a/asmirnov/sim_10.cpp ./
cp /afs/jinr.ru/user/a/asmirnov/reco_10.cpp ./
cp /afs/jinr.ru/user/a/asmirnov/PID_ana.cpp ./

singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.sif spdroot.py "sim_10.cpp(${N_event}, ${SEED}, ${Signal})" > /dev/null
rm core.*
singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.sif spdroot.py "reco_10.cpp(${SEED})" > /dev/null
rm core.*
singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.sif spdroot.py "PID_ana.cpp(${SEED})" > /dev/null
rm core.*

rm run_${SEED}.root
rm params_${SEED}.root
rm reco_full_${SEED}.root

mv PID_ana.root /eos/user/a/asmirnov/reco_data_lambda_c/sig_data/

exit 0

