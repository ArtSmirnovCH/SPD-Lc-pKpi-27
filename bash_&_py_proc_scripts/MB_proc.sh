#!/bin/sh

#SBATCH --array=9
#SBATCH --job-name=MB-proc 		# Job name
#SBATCH --mail-type=END,FAIL 		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=a_smirnov@jinr.ru 	# Where to send mail
#SBATCH --no-requeue 			# do not auto requeue on errors
#SBATCH --ntasks=1 			# Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb 			# Job memory request
#SBATCH --time=2-00:00:00 		# Time limit days-hrs:min:sec
#SBATCH --tmp=10G 			# $TMPDIR space
#SBATCH --output=none

# 4000 events per job (check)

out_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/background"

if [ ! -d "$out_dir/root" ]; then
    echo "Creating /root directory!"
    mkdir "$out_dir/root"
fi

if [ ! -d "$out_dir/txt" ]; then
    echo "Creating /txt directory!"
    mkdir "$out_dir/txt"
fi

# spdroot-dev-4.1.7.4.sif

# /afs/jinr.ru/user/a/asmirnov/MB_proc.py

CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1

cd ${TMPDIR}
python3 /afs/jinr.ru/user/a/asmirnov/MB_proc_2.py \
    --container-path=/cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif \
    --ana-script=/afs/jinr.ru/user/a/asmirnov/analyze_Lc_pKpi_27_MB.cpp \
    --reco-file-list=/afs/jinr.ru/user/a/asmirnov/txt_data_files/bg_data/reco_list.txt \
    --par-file-list=/afs/jinr.ru/user/a/asmirnov/txt_data_files/bg_data/par_list.txt \
    --files-in-job=40 \
    --output-file-dir=${out_dir}/root \
    --job-id=${SLURM_ARRAY_TASK_ID}

mv CutFlow.txt CutFlow_${SLURM_ARRAY_TASK_ID}.txt
mv info_analysis.txt info_analysis_${SLURM_ARRAY_TASK_ID}.txt
rm -rf *MC2025*.txt
mv info_analysis_*.txt ${out_dir}/txt/
mv CutFlow_*.txt ${out_dir}/txt/
# mv *.txt /afs/jinr.ru/user/a/asmirnov/

exit 0

