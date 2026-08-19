# Scripts Key Points Description

## Table of Contents  
[1. SoftWare](#software)  
[2. sim_27](#sim_27)    
[3. reco_27](#reco_27)  
[4. analyze_Lc_pKpi_27 and analyze_Lc_pKpi_27_MB](#analyze_lc_pkpi_27-and-analyze_lc_pkpi_27_mb)    
[5. pid_dataset_27](#pid_dataset_27)    
[6. pv_dataset_27](#pv_dataset_27)


## SoftWare
* SpdRoot 4.1.7.4
* KFParticle
* Pythia8

## sim_27

### Generator Set Up:
* Beams parameters:
    ```cpp
    P8gen -> SetBeam(2212, 2212, 27.);
    primGen->SetBeam(0., 0., 0.1, 0.1);
    primGen->SmearGausVertexXY(kTRUE);
    primGen->SetTarget(0., 30.);
    primGen->SmearGausVertexZ(kTRUE);
    ```
* Subpprocesses:
    * Signal:
        ```cpp
        P8gen -> SetParameters("HardQCD:gg2ccbar = on");
        P8gen -> SetParameters("HardQCD:qqbar2ccbar = on");
        P8gen->SetParameters("PhaseSpace:pTHatMin = 1.");
        ```
    * Background: 
        ```cpp
        P8gen -> SetParameters("SoftQCD:all = on");
        ```
* Forced decay: $\Lambda_c^{+} \to p^{+}K^{-}\pi^{+}$
* Additional requirements in case of signal:
    * Select events only with $\Lambda_c^+$ production

### Detector Set Up:
* Standard II Phase Configuration:
    * Vertex Detector:
        * DSSD
    * Particle Identification:
        * Time-of-Flight System

## reco_27

* Preselection for tracks parameters:
    ```cpp
    track_finder->CheckMinItsHits(true,0);
    track_finder->CheckMinTsHits(true,6);
    track_finder->CheckMinHits(true,7);
    track_finder->CheckMaxPartGeneration(true, 3);
    track_finder->CheckMinPartPt(true,0.1);
    track_finder->CheckMinPartMomentum(true,0.15);
    ```

## analyze_Lc_pKpi_27 and analyze_Lc_pKpi_27_MB

### PID:

* TOF PID with SoftMax transformation and minimum classification probability threshold.
* If no suitable TOF PID is found, skip the track.

### Tracks for PV selection:

```cpp
    if ( track -> GetNHitsIts() < 3) continue;
    if (trkPos.Perp() > 0.3 || (PV_pos_prefit-trkPos).Mag() > 0.4) continue;
```

### PV reconstruction:

* Improving PV reconstruction using a refitting algorithm

### Tracks for SV selection:

```cpp
    track -> GetNHitsIts() < 3 ... continue
```

## pid_dataset_27
...
## pv_dataset_27
...




