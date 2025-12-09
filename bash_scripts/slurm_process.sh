#!/bin/sh

#SBATCH --job-name=main              # Job name
#SBATCH --mail-user=a_smirnov@jinr.ru   # Where to send mail
#SBATCH --no-requeue                    # do not auto requeue on errors
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb                       # Job memory request
#SBATCH --time=2-00:00:00               # Time limit days-hrs:min:sec
#SBATCH --tmp=10G                       # $TMPDIR space
#SBATCH --output=none


CURLOG=$SLURM_ARRAY_TASK_ID.log
# Redirec stdout & strerr
exec > $CURLOG 2>&1

Signal=1 # signal --- 1 (0-470); background --- 0 ();

if [ $Signal -eq 1 ]; then
    eos_dir="/eos/user/a/asmirnov/reco_data_lambda_c/sig_data/"
    out_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/signal/"
    
else
    eos_dir="/eos/user/a/asmirnov/reco_data_lambda_c/bg_data/"
    out_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/background/"
fi

count=0
for file in ${eos_dir}/reco_full*.root; do
    if [ -f "$file" ]; then
        ((count=count+1))
    fi
done

echo "Number of files with extension .root: $count"

SEED=${SLURM_ARRAY_TASK_ID}

if [ ! -d "$out_dir/root" ]; then
    echo "Creating /root directory!"
    mkdir "$out_dir/root"
fi

if [ ! -d "$out_dir/txt" ]; then
    echo "Creating /txt directory!"
    mkdir "$out_dir/txt"
fi

FILE="${eos_dir}reco_full_${SEED}.root"
if [ ! -f "$FILE" ]; then
    echo "File $FILE doesn't exist, stop process!"
    exit
fi

echo "Analyse data from ${FILE}"

cd ${TMPDIR}

N_max=10000 # Max number of events in each file to process;

cp /afs/jinr.ru/user/a/asmirnov/analyze_Lc_pKpi_27_cluster.cpp ./

singularity run -H ./:/WORKDIR --bind /eos/user/a/asmirnov/reco_data_lambda_c:/WORKDIR /cvmfs/spd.jinr.ru/images/spdroot-dev-4.1.7.4.sif spdroot.py "analyze_Lc_pKpi_27_cluster.cpp(${Signal}, ${SEED}, ${N_max})" > /dev/null

rm core.*      

if [ $Signal -eq 1 ]; then
	mv analysed_signal.root analysed_signal_${SEED}.root
else
	mv analysed_background.root analysed_background_${SEED}.root
fi


mv analysed*.root ${out_dir}/root/

mv info_analysis.txt info_analysis_${SEED}.txt
mv CutFlow.txt CutFlow_${SEED}.txt
mv *.txt ${out_dir}/txt/

exit 0

