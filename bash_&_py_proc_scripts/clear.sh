#!/bin/sh

if ls slurm*.out 1> /dev/null 2>&1; then
    rm slurm*.out
    echo "Directory is clear!"
else
    echo "No slurm_info files in directory"
fi

if ls *.log 1> /dev/null 2>&1; then
    rm *.log
    echo "Directory is clear!"
else
    echo "No slurm_info files in directory"
fi

exit 0
