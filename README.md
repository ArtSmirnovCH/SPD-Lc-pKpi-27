# Feasibility Study of $\Lambda_c^{+}$ Production in pp Collisions with the SPD Stage II Detector

## Table of Contents  
[1. Project Description](.README.md#Project-Description)  
[2. Objectives](.README.md#Objectives)  
[3. Data Description](.README.md#Data-Description)  
[4. Stages of Work on the Project](.README.md#Stages-of-Work-on-the-Project)   
[5. Conclusions](.README.md#Conclusions) 

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
Detailed Description can be found here[LINK].


### Project Structure:
```console
├── ./
│   ├── .gitignore
│   ├── README.md
│   └── requirements.txt
│   ├── Analysis/
│   │   ├── EDA.ipynb
│   │   ├── PID_analysis.ipynb
│   │   ├── fits.ipynb
│   │   ├── main_analysis.ipynb
│   │   ├── scripts.py
│   │   └── temp.ipynb
│   ├── plots/
│   ├── Data/
│   │   ├── DataDescription.md
│   ├── Simulation_and_Modeling/
│   │   ├── PID_ana.cpp
│   │   ├── analyze_Lc_pKpi_27_cluster.cpp
│   │   ├── analyze_Lc_pKpi_27_cluster_tree_MB.cpp
│   │   ├── reco_event_lambda_new.cpp
│   │   └── sim_lambda_new.cpp
```

### Conclusions:

[Main Analysis (ipynb)](https://github.com/ArtSmirnovCH/SPD-Lc-pKpi-27/blob/master/Analysis/main_analysis.ipynb)