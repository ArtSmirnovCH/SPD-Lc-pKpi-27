#!/bin/sh

# Calculate signal raw events
data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/sig_data/"

single_event_size=$(cat ${data_dir}N_event.txt)
files_amount=$(ls -1 ${data_dir}reco_full*.root | grep -v "^$" | wc -l)
total_sum=$(($single_event_size * $files_amount))
echo "Total number of raw signal events: $total_sum"

# Calculate background raw events
data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/bg_data/"

single_event_size=$(cat ${data_dir}N_event.txt)
files_amount=$(ls -1 ${data_dir}reco_full*.root | grep -v "^$" | wc -l)
total_sum=$(($single_event_size * $files_amount))
echo "Total number of raw background events: $total_sum"

# Calculate signal processed events
data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/signal/txt/"
#data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data_storage/OnCut_Lc_pkpi_27/signal//txt/"
total_sum=0
for file in $data_dir/info_analysis_*.txt; do
	if [ -f "$file" ]; then
	    while read number; do
		total_sum=$((total_sum + number))
	    done < "$file"
	fi
done

echo "Total number of processed(analysed) signal events: $total_sum"
touch "sig_events.txt"
rm sig_events.txt
echo "$total_sum" >> "sig_events.txt"

# Calculate background processed events
data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/background/txt/"
#data_dir="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data_storage/OnCut_Lc_pkpi_27/background/txt/"
total_sum=0
for file in $data_dir/info_analysis_*.txt; do
	if [ -f "$file" ]; then
	    while read number; do
		total_sum=$((total_sum + number))
	    done < "$file"
	fi
done

echo "Total number of processed(analysed) background events: $total_sum"
#touch "bg_events.txt"
rm bg_events.txt
echo "$total_sum" >> "bg_events.txt"

exit

