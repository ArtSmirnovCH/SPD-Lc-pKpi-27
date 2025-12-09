#! /bin/bash

data_type=1 # signal --- 1; background --- 0;

if [ $data_type -eq 1 ]; then
    data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/sig_data/"
else
    data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/bg_data/"
fi

cd $data_dir

j=0
for file in reco_full*.root; do
    if [ $j -gt $(echo "$file" | grep -o '[0-9]' | tr -d '\n') ]; then
	continue
    fi
    if echo "$file" | grep -o '[0-9]\+' | awk '{if ($1 < 500) exit 1}'; then       
	continue
    fi
    i=$((j))
    while [ $i -lt 500 ]; do
        new_name="reco_full_${i}.root"
        if [ "$file" = "$new_name" ]; then
	    j=$((i + 1))
            break
        fi
        if [ -e "$new_name" ]; then
            ((i++))
            new_name="reco_full_${i}.root"
        else
	    j=$((i + 1))
            mv "$file" "$new_name"
            break
        fi
    done
done 

j=500
for file in reco_full*.root; do
    if [ $j -gt $(echo "$file" | grep -o '[0-9]' | tr -d '\n') ]; then
        continue
    fi
    if echo "$file" | grep -o '[0-9]\+' | awk '{if ($1 > 499) exit 1}'; then
        continue
    fi
    i=$((j))
    while [ $i -lt 1000 ]; do
        new_name="reco_full_${i}.root"
        if [ "$file" = "$new_name" ]; then
            j=$((i + 1))
            break
        fi
        if [ -e "$new_name" ]; then
            ((i++))
            new_name="reco_full_${i}.root"
        else
            j=$((i + 1))
            mv "$file" "$new_name"
            break
        fi
    done
done

echo "Sorting DONE!"
echo "Max file number is ${i}."

exit
