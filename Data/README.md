# Data Description

[1. analysed_signal and analysed_background](#analysed_signal-and-analysed_background)  
[2. ana_PV_27_var_1 and ana_PV_27_var_2](#ana_pv_27_var_1-and-ana_pv_27_var_2)  
[3. pid_ana](#pid_ana)


## analysed_signal and analysed_background:

### Main info:
* Number of Features: 71
* Number of Entries: 
    * Signal: 415460
    * Background: 1326210

### Features Description:

#### Event Information:
* n_event — Event identifier number from modeling (don't use. Contains inconsistent duplicates)

* id — Unique event identifier (contains duplicates for the same PVs)

* multiplicity — Number of reconstructed tracks in the event

#### Primary Particle Kinematics:
* mass_Lc — Reconstructed invariant mass of $\Lambda_c^+$ candidate

* P_Lc — Total momentum magnitude of reconstructed $\Lambda_c^+$

* Pt_Lc — Transverse momentum of reconstructed $\Lambda_c^+$

* eta_Lc — Pseudorapidity of reconstructed $\Lambda_c^+$

* xF — Feynman-x variable of $\Lambda_c^+$ (longitudinal momentum fraction)

* phi_angle — Azimuthal production angle ($\phi$) of $\Lambda_c^+$ in the laboratory frame

#### Daughter Particle Kinematics:

Proton ($p^+$) Candidate:

* P_p — Total momentum magnitude

* Pt_p — Transverse momentum

* eta_p — Pseudorapidity

* OA_p — Opening angle between proton and $\Lambda_c^+$ momentum direction at $\Lambda_c^+$ frame.

Pion ($\pi^+$) Candidate:

* P_pip — Total momentum magnitude

* Pt_pip — Transverse momentum

* eta_pip — Pseudorapidity

* OA_pip — Opening angle between pion and $\Lambda_c^+$ momentum direction at $\Lambda_c^+$ frame.

Kaon ($K^⁻$) Candidate:

* P_K — Total momentum magnitude

* Pt_K — Transverse momentum

* eta_K — Pseudorapidity

* OA_K — Opening angle between kaon and $\Lambda_c^+$ momentum direction at $\Lambda_c^+$ frame.

#### Decay Vertex Properties (Transverse Plane XY):
* lengthXY_Lc — $\Lambda_c^+$ decay length in the transverse plane

* dlengthXY_Lc — Uncertainty of $\Lambda_c^+$ decay length in transverse plane

* ctau_Lc — ($c\tau$) of $\Lambda_c^+$

* chi2_Lc — $\chi^2$ of $\Lambda_c^+$ decay vertex fit

#### Geometric Consistency Variables (XY-plane):
* cosAngle_r_Lc_momentum_Lc_xy — ~~Cosine~~Just angle between $\Lambda_c^+$ flight vector and reconstructed momentum in transverse plane

* cosAngle_r_Lc_sum_momentum_xy — ~~Cosine~~Just angle between $\Lambda_c^+$flight vector and sum of daughter momenta in transverse plane

* cosAngle_momentum_Lc_sum_momentum_xy — ~~Cosine~~Just angle between $\Lambda_c^+$ momentum and sum of daughter momenta in transverse plane

* ptOverE — $\Lambda_c^+$ transverse momentum over sum energy of tracks with $\theta$ > 0.174533 rad

#### Vertex Association Parameters (Transverse Plane XY):

$\Lambda_c^+$ - Primary Vertex Relations:

* chi2_Lc_PV_xy — Impact parameter $\chi^2$ of $\Lambda_c^+$ trajectory to PV in transverse plane

* dist_Lc_PV_xy — Distance of $\Lambda_c^+$ trajectory to PV in transverse plane

* dist_Lc_PV_xy_custom — Alternative distance calculation method for $\Lambda_c^+$-PV separation

Daughter - Primary Vertex Impact Parameters:

* chi2_p_PV_xy — Proton impact parameter $\chi^2$ to PV in transverse plane

* dist_p_PV_xy — Proton distance to PV in transverse plane

* dist_p_PV_xy_custom — Alternative proton-PV distance calculation

* chi2_K_PV_xy — Kaon impact parameter $\chi^2$ to PV in transverse plane

* dist_K_PV_xy — Kaon distance to PV in transverse plane

* dist_K_PV_xy_custom — Alternative kaon-PV distance calculation

* chi2_pip_PV_xy — Pion impact parameter $\chi^2$ to PV in transverse plane

* dist_pip_PV_xy — Pion distance to PV in transverse plane

* dist_pip_PV_xy_custom — Alternative pion-PV distance calculation

Daughter - $\Lambda_c^+$ Vertex Relations:

* chi2_p_Lc_xy — Proton impact parameter $\chi^2$ to $\Lambda_c^+$ decay vertex in transverse plane

* dist_p_Lc_xy — Proton distance to $\Lambda_c^+$ decay vertex in transverse plane

* dist_p_Lc_xy_custom — Alternative proton-$\Lambda_c^+$ vertex distance calculation

* chi2_K_Lc_xy — Kaon impact parameter $\chi^2$ to $\Lambda_c^+$ decay vertex in transverse plane

* dist_K_Lc_xy — Kaon distance to $\Lambda_c^+$ decay vertex in transverse plane

* dist_K_Lc_xy_custom — Alternative kaon-$\Lambda_c^+$ vertex distance calculation

* chi2_pip_Lc_xy — Pion impact parameter $\chi^2$ to $\Lambda_c^+$ decay vertex in transverse plane

* dist_pip_Lc_xy — Pion distance to $\Lambda_c^+$ decay vertex in transverse plane

* dist_pip_Lc_xy_custom — Alternative pion-$\Lambda_c^+$ vertex distance calculation

Track-Track Correlations (XY-plane):

* chi2_K_pip_xy — Impact parameter $\chi^2$ between kaon and pion in transverse plane

* dist_K_pip_xy — Distance between kaon and pion in transverse plane

* dist_K_pip_xy_custom — Alternative kaon-pion distance calculation

* chi2_p_K_xy — Impact parameter $\chi^2$ between proton and kaon in transverse plane

* dist_p_K_xy — Distance between proton and kaon in transverse plane

* dist_p_K_xy_custom — Alternative proton-kaon distance calculation

* chi2_p_pip_xy — Impact parameter $\chi^2$ between proton and pion in transverse plane

* dist_p_pip_xy — Distance between proton and pion in transverse plane

* dist_p_pip_xy_custom — Alternative proton-pion distance calculation

#### Monte Carlo Truth Information:

* true_decay — Binary flag indicating true $\Lambda_c^+$ decay reconstruction (1 = true signal, 0 = background)

#### Vertex Position Differences (MC Truth):

$\Lambda_c^+$ Vertex Resolution (only valid for true_decay == 1):

* Lc_diff_x — $\Delta x$ between true and reconstructed $\Lambda_c^+$ vertex positions

* Lc_diff_y — $\Delta y$ between true and reconstructed $\Lambda_c^+$ vertex positions

* Lc_diff_z — $\Delta z$ between true and reconstructed $\Lambda_c^+$ vertex positions

Primary Vertex Resolution:
Before PV Re-reconstruction:

* PV_diff_x — $\Delta x$ between true and reconstructed primary vertex

* PV_diff_y — $\Delta y$ between true and reconstructed primary vertex

* PV_diff_z — $\Delta z$ between true and reconstructed primary vertex

After Iterative PV Re-reconstruction:

* PV_diff_ES_x — $\Delta x$ between true and re-reconstructed primary vertex

* PV_diff_ES_y — $\Delta y$ between true and re-reconstructed primary vertex

* PV_diff_ES_z — $\Delta$ between true and re-reconstructed primary vertex

PV Reconstruction Quality Indicators:

* kf_pv_size_before_re_fit — Number of tracks used in primary vertex fit before re-reconstruction

* kf_pv_size_after_re_fit — Number of tracks used in primary vertex fit after re-reconstruction

## ana_PV_27_var_1 and ana_PV_27_var_2: 

### Main info:
* Number of Features: 14
* Number of Entries: ~255600

### Features Description:

#### Event Information:

* n_event — Event identifier number (don't use. contains inconsistent duplicates)

#### Primary Vertex Position Resolution (MC Truth Comparisons):

Original PV Reconstruction:

* PV_diff_x — Difference in x-coordinate between Monte Carlo truth primary vertex and originally reconstructed PV (true - reconstructed)

* PV_diff_y — Difference in y-coordinate between Monte Carlo truth primary vertex and originally reconstructed PV (true - reconstructed)

* PV_diff_z — Difference in z-coordinate between Monte Carlo truth primary vertex and originally reconstructed PV (true - reconstructed)

Enhanced PV Reconstruction (After Re-fitting):

* PV_diff_ES_x — Difference in x-coordinate between Monte Carlo truth primary vertex and enhanced/re-reconstructed PV (true - enhanced)

* PV_diff_ES_y — Difference in y-coordinate between Monte Carlo truth primary vertex and enhanced/re-reconstructed PV (true - enhanced)

* PV_diff_ES_z — Difference in z-coordinate between Monte Carlo truth primary vertex and enhanced/re-reconstructed PV (true - enhanced)

#### Vertex Fit Quality Parameters:

* chi2_PV — Total $\chi^2$ value of the primary vertex fit

* ndf_PV — Number of degrees of freedom in the primary vertex fit

* chi2_per_ndf_PV — Normalized $\chi^2$ ($\chi^2$/ndf) indicating goodness-of-fit for the primary vertex

#### Vertex Reconstruction Control Variables:

* max_chi2_cut — Maximum allowed $\chi^2$ contribution per track used as a cut in the iterative vertex fitting procedure

* max_chi2_observed — Maximum $\chi^2$ contribution from any single track observed during the vertex fitting

* num_removed_tracks — Number of tracks removed from the vertex fit during the iterative cleaning procedure

* tracks_per_vertex — Average number of tracks associated with the primary vertex per event


## pid_ana: 

### Main info:
* Number of Features: 22
* Number of Entries: 37418

### Features Description:

#### Monte Carlo Truth Information:

* mc_pid — Monte Carlo truth particle ID code

#### Likelihood-Based PID Variables:

Time-of-Flight (TOF) Likelihoods:

* ll_tof_pi — Log-likelihood for pion hypothesis from Time-of-Flight detector

* ll_tof_k — Log-likelihood for kaon hypothesis from Time-of-Flight detector

* ll_tof_p — Log-likelihood for proton hypothesis from Time-of-Flight detector

dE/dx Likelihoods:

* ll_ts_pi — Log-likelihood for pion hypothesis from dE/dx measurements

* ll_ts_k — Log-likelihood for kaon hypothesis from dE/dx measurements

* ll_ts_p — Log-likelihood for proton hypothesis from dE/dx measurements

#### Track Kinematics:

* px — Momentum component in x-direction

* py — Momentum component in y-direction

* pz — Momentum component in z-direction

* p — Total momentum magnitude

* pt — Transverse momentum

* eta — Pseudorapidity

* phi — Azimuthal angle in transverse plane

* theta — Polar angle relative to beam axis

#### Track Quality Parameters:

* chi2overndf — $\chi^2$ per degree of freedom of the track fit (normalized track fit quality)

* nhitsits — Number of hits in the Inner Tracking System (ITS)

* nhitstsb — Number of hits in Time-of-Flight barrel

* nhitstsec — Number of hits in Time-of-Flight end-cups

* isfitted — Boolean flag indicating whether the track was successfully fitted

* isgood — Boolean flag indicating overall track quality

* convergency — Flag or value indicating convergence status of the track fitting algorithm