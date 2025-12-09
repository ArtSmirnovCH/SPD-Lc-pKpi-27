#!/bin/sh

#SBATCH --tmp=10G # $TMPDIR space
#SBATCH --output=none


CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1

N_event=4500
Explanation="N_event in each .root file in directory"
Signal=1

SEED=${SLURM_ARRAY_TASK_ID}


cd ${TMPDIR}

cp /afs/jinr.ru/user/a/asmirnov/sim_lambda_new.cpp ./
cp /afs/jinr.ru/user/a/asmirnov/reco_event_lambda_new.C ./

singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "sim_lambda_new.cpp(${N_event}, ${SEED}, ${Signal})" > /dev/null
rm core.*
singularity run -H ./:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "reco_event_lambda_new.C(${SEED})" > /dev/null
rm core.*
rm run_${SEED}.root
rm params_${SEED}.root

if [ $Signal -eq 1 ]; then
    mv reco_full_${SEED}.root /eos/user/a/asmirnov/reco_data_lambda_c/sig_data/
    echo "$N_event" > /eos/user/a/asmirnov/reco_data_lambda_c/sig_data/N_event.txt
    echo "$Explanation" > /eos/user/a/asmirnov/reco_data_lambda_c/sig_data/N_event_exmplane.txt
else
    mv reco_full_${SEED}.root /eos/user/a/asmirnov/reco_data_lambda_c/bg_data/
    echo "$N_event" > /eos/user/a/asmirnov/reco_data_lambda_c/bg_data/N_event.txt
    echo "$Explanation" > /eos/user/a/asmirnov/reco_data_lambda_c/bg_data/N_event_exmplane.txt
fi

exit 0
