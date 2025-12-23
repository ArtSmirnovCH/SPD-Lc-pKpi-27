import os
import subprocess
import argparse
import hashlib


def get_job_files(file_list_path, files_in_job, job_n):
    file_list = []

    with open(file_list_path, "r") as f:
        file_list = f.read().splitlines()

    nfiles = len(file_list)

    return file_list[job_n * files_in_job : min((job_n + 1)*files_in_job, nfiles)]

# Get arguments
parser = argparse.ArgumentParser("MB_processor")
parser.add_argument("--container-path", type=str, help="path to the .sif container")
parser.add_argument("--ana-script", type=str, help="script in the .C format")
parser.add_argument("--reco-file-list", type=str, help="path to the text files with reco-filenames")
parser.add_argument("--par-file-list", type=str, help="path to the text files with par-filenames")
parser.add_argument("--files-in-job", type=int, help="maximum reco-files in job")
parser.add_argument("--output-file-dir", type=str, help="output file directory")
parser.add_argument("--job-id", type=int, default=0, help="job id (overrinden by SLURM_ARRAY_TASK_ID if defined")
args = parser.parse_args()
print(args)

job_id = args.job_id

# Get list of files for the job
reco_list = get_job_files(args.reco_file_list, args.files_in_job, job_id)
if not reco_list:
    print("no reco-files for these job, terminating the execution")
    quit()

par_file = get_job_files(args.par_file_list, 1, 0)[0]
if not par_file:
    print("param-file not found, terminating the execution")
    quit()

# Copy the par file
subprocess.call(["xrdcp", par_file, "params.root"])

# Copy script
subprocess.call(["cp", args.ana_script, "./"])
script_name = args.ana_script.split("/")[-1]

# Create temporary output directory
subprocess.call(["mkdir", "-p", "out"])

# Process reco files in loop
ana_files = []
txt_info_files = []
txt_cutFlow_files = []
for rfile in reco_list:
    filename = rfile.split("/")[-1]
    
    # Remove production version due to SpdRoot bug
    if filename[-2] == ".":
        filename = filename[:-2]

    subprocess.call(["xrdcp", rfile, filename])

    ofilename = "out/ana_" + filename

    # Generate integer ID from MD5 hash (first 8 bytes as integer)
    hash_bytes = hashlib.md5(rfile.encode()).digest()[:8]  # Get first 8 bytes
    file_id = int.from_bytes(hash_bytes, byteorder='big')  # Convert to integer
    
    # Optional: limit to reasonable size
    file_id = file_id % (10**5)  # 5-digit integer

    script_string = script_name + '({fid}, \\"{ifile}\\", \\"{ofile}\\")'.format(
        fid=file_id,
        ifile=filename, 
        ofile=ofilename
    )
    
    FNULL = open(os.devnull, 'w')
    result = subprocess.call(["singularity", "run", "-H", "./:/WORKDIR", args.container_path,
                     "spdroot.py", "-b", "-q", script_string], stdout=FNULL)

    # Check if info_analysis.txt and CutFlow.txt were created
    info_file = "info_analysis.txt"
    cutflow_file = "CutFlow.txt"
    
    if os.path.exists(info_file):
        new_info_file = f"info_analysis_{filename}.txt"
        subprocess.call(["mv", info_file, new_info_file])
        txt_info_files.append(new_info_file)
    else:
        print(f"Warning: {info_file} not created for {filename}")
        
    if os.path.exists(cutflow_file):
        new_cutflow_file = f"CutFlow_{filename}.txt"
        subprocess.call(["mv", cutflow_file, new_cutflow_file])
        txt_cutFlow_files.append(new_cutflow_file)
    else:
        print(f"Warning: {cutflow_file} not created for {filename}")
    
    ana_files.append(ofilename)
    subprocess.call(["rm", "-f", filename])

jobout_filename = "analysed_background_" + str(job_id) + ".root"

# Merge root files
if ana_files:
    subprocess.call(["hadd", "-f",  jobout_filename] + ana_files)
    # Remove left root files
    subprocess.call(["rm", "-f", "out/ana_*.root"])
else:
    print("Warning: No root files were created")

# Merging info_analysis.txt
if txt_info_files:
    sums = []
    for filename in txt_info_files:
        try:
            with open(filename, 'r') as f:
                for i, line in enumerate(f):
                    number = int(line.strip())
                    if i >= len(sums):
                        sums.append(number)
                    else:
                        sums[i] += number
        except FileNotFoundError:
            print(f"Warning: {filename} not found, skipping")
            continue

    with open('info_analysis.txt', 'w') as outfile:
        for s in sums:
            outfile.write("{}\n".format(s))
    
    subprocess.call(["rm", "-f", "info_analysis_*.txt"])
else:
    print("Warning: No info_analysis files were created")

# Merging CutFlow.txt
if txt_cutFlow_files:
    sums = []
    for filename in txt_cutFlow_files:
        try:
            with open(filename, 'r') as f:
                for i, line in enumerate(f):
                    number = int(line.strip())
                    if i >= len(sums):
                        sums.append(number)
                    else:
                        sums[i] += number
        except FileNotFoundError:
            print(f"Warning: {filename} not found, skipping")
            continue

    with open('CutFlow.txt', 'w') as outfile:
        for s in sums:
            outfile.write("{}\n".format(s))
    
    subprocess.call(["rm", "-f", "CutFlow_*.txt"])
else:
    print("Warning: No CutFlow files were created")

# Copy root result to eos
if os.path.exists(jobout_filename):
    subprocess.call(["cp", jobout_filename, args.output_file_dir])
else:
    print(f"Warning: {jobout_filename} not created, nothing to copy")
