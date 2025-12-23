#!/bin/sh

#SBATCH --array=4
#SBATCH --job-name=test 		# Job name
#SBATCH --mail-type=NONE 		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=a_smirnov@jinr.ru 	# Where to send mail
#SBATCH --no-requeue 			# do not auto requeue on errors
#SBATCH --ntasks=1 			# Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb 			# Job memory request
#SBATCH --time=0-03:00:00 		# Time limit days-hrs:min:sec
#SBATCH --tmp=10G 			# $TMPDIR space
#SBATCH --output=none

CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1


cd ${TMPDIR}

cp /afs/jinr.ru/user/a/asmirnov/MC_script_pK.C ./

singularity run -H ./:/WORKDIR --bind /eos/user/a/asmirnov/reco_data_lambda_c:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "MC_script_pK.C(${SLURM_ARRAY_TASK_ID})" > /dev/null

rm core.*
mv check.txt MC_check_pK_${SLURM_ARRAY_TASK_ID}.txt
echo "job_id=${SLURM_ARRAY_TASK_ID}" >> MC_check_pK_${SLURM_ARRAY_TASK_ID}.txt

mv *.txt /afs/jinr.ru/user/a/asmirnov

exit 0
