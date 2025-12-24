# Feasibility Study of $\Lambda_c^{+}$ Production in pp Collisions with the SPD Stage II Detector

## Table of Contents  
[1. Project Description](.README.md#project-description)  
[2. Objectives](.README.md#objectives)  
[3. Data Description](.README.md#data-description)  
[4. Stages of Work on the Project](.README.md#stages-of-work-on-the-project) 
[5. Project Structure](.README.md#project-structure)  
[6. Conclusions](.README.md#conclusions) 

### Project Description    
The project involves the complete chain of generation, detector simulation, and event analysis for the full-scale simulation of the production of the charmed $\Lambda_c^+$ baryon under the conditions of the second stage of the SPD experiment at the NICA collider. The decay channel under study is $\Lambda_c^{+} \to p^{+}K^{-}\pi^{+}$. 

### Objectives:
* Conduct reliable simulation of the production of the charmed $\Lambda_c^+$ baryon.

* Perform accurate simulation of the SPD detector in the configuration of the second stage of the experiment.

* Create a dataset of signal and background events.

* Develop a signal selection algorithm.

* Evaluate the feasibility of observing and studying the properties of the $\Lambda_c^+$ baryon with the SPD detector.

* To assess the feasibility of measuring Transverse Single Spin Asymmetry of $\Lambda_c^+$ production by evaluating statistical uncertainties.


### Data Description:
Data was created at cluster JINR (CIVC).
Detailed description can be found in Data directory [here](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/Data/DataDescription.md).

### Stages of Work on the Project:
1. Set up a correct simulation of Lc production and the detector model.
[Simulation (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/tree/master/Simulation_and_Modeling)

2. Create an optimal algorithm for particle identification.
[PID analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/PID_analysis.ipynb)

3. Create an optimal algorithm for primary vertex reconstruction.
[PV analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/pv_analysis.ipynb) & 
[Lc decay length resolution analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/Lc_decay_length_analysis.ipynb)

4. Create a preselection function to remove data artifacts.
[EDA (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/eda.ipynb)

5. Explore the simulated and preselected data.
[EDA (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/eda.ipynb)

6. Create a cut-based selection algorithm.
[Cut-Based Analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/cut_based_analysis.ipynb)

7. Evaluate the performance of the cut-based approach.
[Cut-Based Analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/analysis_notebooks/cut_based_analysis.ipynb)

8. MVA... (Investigate/Implement a Multivariate Analysis approach)

### Project Structure:
```console
├── ./
│   ├── .gitignore
│   ├── README.md
│   └── requirements.txt
│   ├── bash_&_py_proc_scripts/
│   │   ├── CutFlow.sh
│   │   ├── MB_proc.sh
│   │   ├── MB_proc_2.py
│   │   ├── check_empty_analysis_files.sh
│   │   ├── check_missed_files.sh
│   │   ├── clear.sh
│   │   ├── count_events.sh
│   │   ├── merge_processed_root_files_tfile.py
│   │   ├── pid_dataset_sim_reco_ana.sh
│   │   ├── pv_dataset_sim_reco_ana.sh
│   │   ├── slurm_MC_check.sh
│   │   ├── slurm_process.sh
│   │   ├── slurm_simu_reco.sh
│   │   └── sorting_filenames.sh
│   ├── plots/
│   ├── analysis_notebooks/
│   │   ├── Lc_decay_length_analysis.ipynb
│   │   ├── PID_analysis.ipynb
│   │   ├── cut_based_analysis.ipynb
│   │   ├── eda.ipynb
│   │   ├── fits.ipynb
│   │   ├── pv_analysis.ipynb
│   │   └── temp.ipynb
│   ├── Data/
│   │   ├── DataDescription.md
│   ├── analysis_scripts/
│   │   ├── cfg.py
│   │   ├── draw_scripts.py
│   │   ├── estimate_scripts.py
│   │   ├── fit_scripts.py
│   │   └── selection_scripts.py
│   ├── Simulation_and_Modeling/
│   │   ├── analyze_Lc_pKpi_27.cpp
│   │   ├── analyze_Lc_pKpi_27_MB.cpp
│   │   ├── pid_dataset_27.cpp
│   │   ├── pv_dataset_27.cpp
│   │   ├── reco_27.cpp
│   │   └── sim_27.cpp

```

### Conclusions:

...