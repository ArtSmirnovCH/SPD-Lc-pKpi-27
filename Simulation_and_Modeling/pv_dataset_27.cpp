#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numeric>


#define HomogeneousField 


std::vector<double> softmax(const std::vector<double>& array) {
    if (array.empty()) {
        return {};
    }
    
    // Find maximum value
    double max_val = *std::max_element(array.begin(), array.end());
    
    // Calculate exponentials after subtracting max (for numerical stability)
    std::vector<double> exp_vals(array.size());
    std::transform(array.begin(), array.end(), exp_vals.begin(),
                   [max_val](double x) { return std::exp(x - max_val); });
    
    // Calculate sum of exponentials
    double sum_exp = std::accumulate(exp_vals.begin(), exp_vals.end(), 0.0);
    
    // Normalize by the sum
    std::vector<double> result(exp_vals.size());
    std::transform(exp_vals.begin(), exp_vals.end(), result.begin(),
                   [sum_exp](double exp_val) { return exp_val / sum_exp; });
    
    return result;
}


int choice_with_thresholds(const std::vector<double>& probabilities,
                          const std::vector<double>& thresholds = {}) {

    // Handle default thresholds if empty
    std::vector<double> actual_thresholds = thresholds;
    if (actual_thresholds.empty()) {
        actual_thresholds = std::vector<double>(probabilities.size(), 0.0);
    }

    // Check that vectors have the same size
    if (probabilities.size() != actual_thresholds.size()) {
        return 0;
    }

    // Check for empty input
    if (probabilities.empty()) {
        return 0;
    }

    // Check for invalid inputs
    for (size_t i = 0; i < probabilities.size(); ++i) {
        // Check for NaN or infinity
        if (std::isnan(probabilities[i]) || std::isinf(probabilities[i])) {
            return 0;
        }
        // Check for negative probabilities (if not allowed)
        if (probabilities[i] < 0.0) {
            return 0;
        }
    }

    // Apply threshold mask - check if any probability passes its threshold
    bool passes_threshold = false;
    for (size_t i = 0; i < probabilities.size(); ++i) {
        if (probabilities[i] > actual_thresholds[i]) {
            passes_threshold = true;
            break;
        }
    }

    // Return 0 if no threshold is passed
    if (!passes_threshold) {
        return 0;
    }

    // Create a working copy to avoid modifying the original
    std::vector<double> adjusted_probs = probabilities;

    // Subtract thresholds from probabilities that pass
    for (size_t i = 0; i < adjusted_probs.size(); ++i) {
        if (adjusted_probs[i] > actual_thresholds[i]) {
            adjusted_probs[i] -= actual_thresholds[i];
        } else {
            // Set probabilities below threshold to a very low value
            // to ensure they don't get selected as maximum
            adjusted_probs[i] = -1.0;
        }
    }

    // Find the maximum probability
    auto max_it = std::max_element(adjusted_probs.begin(), adjusted_probs.end());

    // Additional safety check - ensure the maximum value is valid
    if (max_it == adjusted_probs.end() || *max_it < 0.0) {
        return 0;
    }

    size_t max_index = std::distance(adjusted_probs.begin(), max_it);

    // Return the corresponding particle type - now dynamic based on index
    // You might want to adjust this mapping for your specific use case
    switch (max_index) {
        case 0: return 211;   // π⁺
        case 1: return 321;   // K⁺
        case 2: return 2212;  // p
        default: return 0;
    }
}


const TClonesArray *mcparticles{};
const TClonesArray *mctracks{};
const TClonesArray *rcvertices{};
const TClonesArray *mcvertices{};
const TClonesArray *particles_tof{}; 
const TClonesArray *particles_ts{};

//======================================================================================================================
// Selecting tracks that are used for primary vertex reconstruction 
//======================================================================================================================
TVector3 truePVPosition;                                        // True postion of generated vertex
Bool_t findPVCheck{0};                                          // If we found true generated PV        
Bool_t onlyPrim{0};                                       	// Use only particles from primary vertex

std::vector<KFParticle> KFParticles_PV;                         // Selected tracks for PV
TVector3 PV_pos_prefit{};                                       // PV RC position first approximation.

void SelectTracksForPV( TVector3 PV_pos_prefit )
{
    KFParticles_PV.clear();

    findPVCheck = 0;

    SpdTrackPropagatorGF fPropagator{};
    fPropagator.Init();

    for (Int_t i{}; i < mctracks -> GetEntriesFast(); ++i)
    {
        SpdTrackMC *track{ (SpdTrackMC*) mctracks -> At(i) };
        if ( !track ) continue;
        Int_t charge{ track -> GetParticlePdg() / abs(track -> GetParticlePdg()) };										
        SpdMCParticle *particle{ (SpdMCParticle*) mcparticles -> At( track -> GetParticleId() ) };
        if ( particle -> GetMCTrackId() < 0 ) continue;
        if ( particle -> GetGeneration() == 1 && !findPVCheck )   // GetGeneration() == 1 if track starts in PV
        {
                truePVPosition  = particle -> GetStartPos();
                findPVCheck = 1;
        }
        if ( onlyPrim && particle -> GetGeneration() != 1 ) continue;   // Use only tracks from PV
        if ( !track -> GetIsFitted() ) continue;
        SpdTrackFitPar *tpars{ track -> GetFitPars() };
        if ( !tpars ) continue;
        const SpdTrackState *state{ tpars -> GetFirstState() };
        if ( !state ) continue;
        // if ( tpars -> GetChi2overNDF() > 12.) continue;
        if ( track -> GetNHitsIts() < 3) continue;
        // if ( track -> GetNHitsTsB() + track -> GetNHitsTsEC() < 6) continue;
        // if ( track -> GetNHitsIts() + track -> GetNHitsTsB() + track -> GetNHitsTsEC() < 8) continue;
        //======================================================================================================================
        // PID
        //======================================================================================================================
        Int_t pdgTrk{211};
					
		// TOF
		std::vector<Double_t> Likelihoods_tof;              
		if ( particle -> GetTofParticleId() != -1 ){
			SpdTofParticle *tofparticle = dynamic_cast<SpdTofParticle*>( particles_tof -> At( particle -> GetTofParticleId() ) );
			if ( !tofparticle ) continue;

			Likelihoods_tof = tofparticle -> GetLikelihoods();  
		}
		else continue;

		Likelihoods_tof.pop_back();

		std::vector<Double_t> thresholds_tof = {0.4, 0.4, 0.4};

		std::vector<Double_t> softmax_output_tof = softmax(Likelihoods_tof);
		
		Int_t tof_pid{ static_cast<Int_t>(choice_with_thresholds(softmax_output_tof, thresholds_tof)) };

		if ( !tof_pid ) continue;

		pdgTrk = tof_pid;

		// vector <Double_t> Likelihoods_ts;       
		// if ( particle -> GetTsParticleId() != -1 ){   
		//     SpdTsParticle *ftsparticle = dynamic_cast<SpdTsParticle*>( particles_ts -> At( particle -> GetTsParticleId() ) );
		//     if ( !ftsparticle ) continue; 
		//     
		//     Likelihoods_ts = ftsparticle -> GetLikelihoods();
		// }
		// else continue;

		pdgTrk *= charge;
        //======================================================================================================================
        SpdTrackState  stateOut;
        fPropagator.InitTrack( pdgTrk, 0 );	// Propagation direction (-1, 0, 1) -> (backward, auto, forward) ???
        Double_t dist{ fPropagator.ExtrapolateToPoint(PV_pos_prefit, *state, stateOut) };// extrapolate track to 1-st estimation of PV
        if (dist == 0.0) continue;
        TVector3 trkPos = stateOut.GetPosition();
        TVector3 trkMom = stateOut.GetMomentum();
        const TMatrixDSym& trkCov = stateOut.GetCovariance();
        // Remove track which is faraway from vertex and beam line
        // trkPos.Perp() --- distanse to beam line(z axis) [cm]
        // (PV_pos_prefit-trkPos).Mag() --- distanse between PV prefit position and track ?position? [cm]
        if (trkPos.Perp() > 0.3 || (PV_pos_prefit-trkPos).Mag() > 0.4) continue;				// Sel
        KFPTrack kfTrack;
        kfTrack.SetParameters( trkPos.X(), trkPos.Y(), trkPos.Z(), trkMom.X(), trkMom.Y(), trkMom.Z() );

        Double_t C[21] = {
        trkCov(0,0),
        trkCov(1,0),trkCov(1,1),
        trkCov(2,0),trkCov(2,1),trkCov(2,2),
        trkCov(3,0),trkCov(3,1),trkCov(3,2),trkCov(3,3),
        trkCov(4,0),trkCov(4,1),trkCov(4,2),trkCov(4,3),trkCov(4,4),
        trkCov(5,0),trkCov(5,1),trkCov(5,2),trkCov(5,3),trkCov(5,4),trkCov(5,5) };
        kfTrack.SetCovarianceMatrix(C);
        kfTrack.SetNDF( tpars -> GetNDF() );
        kfTrack.SetChi2( tpars -> GetChi2() );
        kfTrack.SetCharge( charge );

        KFParticle p1(kfTrack, pdgTrk);
        KFParticles_PV.push_back(p1);
    }
}

//======================================================================================================================
// Main analysis
//======================================================================================================================
void analyse(Int_t N, Int_t seed) // N - max event number to analyse
{
    SpdMCDataIterator* IT = new SpdMCDataIterator();

    TString inFile = Form("reco_full_%d.root", seed);
    IT -> AddSourceFile(inFile);
    
    IT -> ActivateBranch("all");
    
    IT -> Init();
    
	mcparticles = IT -> GetParticles();    
	mctracks    = IT -> GetTracks();         
	rcvertices  = IT -> GetVerticesRC();
	mcvertices  = IT -> GetVertices();    
	particles_ts = IT -> GetTsParticles();
    particles_tof = IT -> GetTofParticles();
	
	//======================================================================================================================
	// Process variables
	const Int_t events_max{ N };		// Upper limit for events to process
	Int_t n_event{};			// Current event number
	
    Double_t PV_diff_x{};
    Double_t PV_diff_y{};
    Double_t PV_diff_z{};

    Double_t max_chi2_cut;            // Current chi2 cut value being tested
    Int_t num_removed_tracks;         // Number of tracks removed for this cut
    Int_t tracks_per_vertex;          // Number of tracks used in final vertex fit
    Double_t chi2_PV;                 // Chi2 of the primary vertex fit
    Double_t ndf_PV;                  // Number of degrees of freedom for PV
    Double_t chi2_per_ndf_PV;         // Chi2/NDF of PV fit
    Double_t PV_diff_ES_x;            // Difference in x: reco - true PV
    Double_t PV_diff_ES_y;            // Difference in y: reco - true PV
    Double_t PV_diff_ES_z;            // Difference in z: reco - true PV
    Double_t max_chi2_observed;       // Maximum chi2 value observed in final iteration

    Int_t NHits_pi{};
    Int_t NHitsIts_pi{};
    Int_t NHitsTsB_pi{};
    Int_t NHitsTsEC_pi{};
    Int_t Chi2OverNDF_pi{};

    Int_t NHits_pi_true{};
    Int_t NHitsIts_pi_true{};
    Int_t NHitsTsB_pi_true{};
    Int_t NHitsTsEC_pi_true{};
    Int_t Chi2OverNDF_pi_true{};

    Int_t NHits_p{};
    Int_t NHitsIts_p{};
    Int_t NHitsTsB_p{};
    Int_t NHitsTsEC_p{};
    Int_t Chi2OverNDF_p{};

    Int_t NHits_p_true{};
    Int_t NHitsIts_p_true{};
    Int_t NHitsTsB_p_true{};
    Int_t NHitsTsEC_p_true{};
    Int_t Chi2OverNDF_p_true{};

    Int_t NHits_K{};
    Int_t NHitsIts_K{};
    Int_t NHitsTsB_K{};
    Int_t NHitsTsEC_K{};
    Int_t Chi2OverNDF_K{};

    Int_t NHits_K_true{};
    Int_t NHitsIts_K_true{};
    Int_t NHitsTsB_K_true{};
    Int_t NHitsTsEC_K_true{};
    Int_t Chi2OverNDF_K_true{};

    Double_t max_chi2{};                // Temporary variable for max chi2 in loop

    Int_t PV_count{};
    Int_t PV_extra_sel_count{};

	//======================================================================================================================
	// Output File
	TFile* file;
    file = new TFile("ana_PV_27.root", "RECREATE");
	//======================================================================================================================
	// TTree
	TTree *tree = new TTree("tree", "tree");

    tree -> Branch("n_event", &n_event, "n_event/I");

    tree -> Branch("PV_diff_x", &PV_diff_x, "PV_diff_x/D");
    tree -> Branch("PV_diff_y", &PV_diff_y, "PV_diff_y/D");
    tree -> Branch("PV_diff_z", &PV_diff_z, "PV_diff_z/D");

    tree->Branch("max_chi2_cut", &max_chi2_cut, "max_chi2_cut/D");
    tree->Branch("num_removed_tracks", &num_removed_tracks, "num_removed_tracks/I");
    tree->Branch("tracks_per_vertex", &tracks_per_vertex, "tracks_per_vertex/I");
    tree->Branch("chi2_PV", &chi2_PV, "chi2_PV/D");
    tree->Branch("ndf_PV", &ndf_PV, "ndf_PV/D");
    tree->Branch("chi2_per_ndf_PV", &chi2_per_ndf_PV, "chi2_per_ndf_PV/D");
    tree->Branch("PV_diff_ES_x", &PV_diff_ES_x, "PV_diff_ES_x/D");
    tree->Branch("PV_diff_ES_y", &PV_diff_ES_y, "PV_diff_ES_y/D");
    tree->Branch("PV_diff_ES_z", &PV_diff_ES_z, "PV_diff_ES_z/D");
    tree->Branch("max_chi2_observed", &max_chi2_observed, "max_chi2_observed/D");

    //======================================================================================================================
    while ( IT -> NextEvent() && n_event < events_max )
	{
		++n_event;
		Int_t multiplicity{ mctracks -> GetEntriesFast() };
		if ( mctracks -> GetEntriesFast() < 5 ) continue;					// Check if more then 4 tracks in the event.

        //=========================================================================================================================
		// Primary vertex reconstruction
		//=========================================================================================================================
		// PV fit using standard SPDroot approach
		// Then is used for KFP Vertex first approximation
		SpdVertexRC* prim_vtx = (SpdVertexRC*) rcvertices -> At(0);
		if ( !prim_vtx || !prim_vtx -> IsPrimary() ) continue;
		SpdPrimVertexFitPar* prim_vtx_fit = dynamic_cast <SpdPrimVertexFitPar*> (prim_vtx -> GetFitPars()); // PV RC first approximation.
		if (!prim_vtx_fit) continue;
		TVector3 PV_pos{};                                                                      // PV position after fit.
		PV_pos_prefit = prim_vtx_fit -> GetVertex();                                            // PV RC position first approximation.  
        
        ++PV_count;

        SelectTracksForPV( PV_pos_prefit );                                                     // Tracks selection for PV reconstruction via KFP

		if ( !findPVCheck )     continue;                                                       // Skip if no PV.
		if ( KFParticles_PV.size() < 4 ) continue;						// Check if have more than 4 tracks fo PV reconstruction.
		//=========================================================================================================================
		// Set field in PV position 
		SpdTrackPropagatorGF fPropagator;
		fPropagator.Init();
		Float_t fieldBz{};
		SpdField* field = fPropagator.GetField();
		if (field) fieldBz = field->GetBz(PV_pos_prefit.X(),PV_pos_prefit.Y(),PV_pos_prefit.Z());
#ifdef HomogeneousField 
		KFParticle::SetField(fieldBz);
#endif
		//=========================================================================================================================
		// Reconstruct primary vertex with KFParticle
		
        KFParticle primVtx{};
		
        Int_t nDaughters_temp = KFParticles_PV.size();
		
        const KFParticle* vDaughters_temp[nDaughters_temp];
		
        for ( Int_t i{}; i < KFParticles_PV.size(); ++i ) vDaughters_temp[i] = &KFParticles_PV[i];

		primVtx.Construct(vDaughters_temp, nDaughters_temp, 0, -1);     // 0 --- no parant particle, 1 --- no mass hypothesis

		PV_diff_x = primVtx.GetX() - truePVPosition.X();
		PV_diff_y = primVtx.GetY() - truePVPosition.Y();
		PV_diff_z = primVtx.GetZ() - truePVPosition.Z();

        //=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (MAIN ANALYSIS) VAR-1
        //=========================================================================================================================

        // THE CODE FOR MAIN ANALYSIS

        // // Use only ONE optimal cut value (e.g., 5.0 based on your optimization)
        // Double_t optimal_cut = 5.0;  // Change this to your chosen value
        // max_chi2_cut = optimal_cut;
 
        // num_removed_tracks = 0;
        // Bool_t continue_removal = true;
 
        // while (continue_removal && KFParticles_PV.size() > 2) {
        //     Int_t worst_track_idx = -1;
        //     Double_t max_chi2 = 0.0;
        //     
        //     // Find track with maximum chi2 to current vertex
        //     for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
        //         const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
        //         if (chi2 > max_chi2) {
        //             worst_track_idx = i;
        //             max_chi2 = chi2;
        //         }
        //     }
        //     
        //     // Check if we should remove this track
        //     if (worst_track_idx >= 0 && max_chi2 > optimal_cut) {
        //         // Remove the worst track
        //         KFParticles_PV.erase(KFParticles_PV.begin() + worst_track_idx);
        //         
        //         // Recalculate vertex from remaining tracks
        //         primVtx = KFParticle(); // Reset vertex
        //         for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
        //             primVtx += KFParticles_PV[i]; // Add all remaining tracks to rebuild vertex
        //         }
        //         
        //         ++num_removed_tracks;
        //     } else {
        //         continue_removal = false;
        //     }
        // }

        // Now primVtx contains the updated vertex after removing bad tracks
        // KFParticles_PV contains only the tracks that passed the chi2 cut

        //=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (MAIN ANALYSIS) VAR-2
        //=========================================================================================================================

        // THE CODE FOR MAIN ANALYSIS

        // Double_t optimal_cut = 5.0;  // Change this to your chosen value
        // max_chi2_cut = optimal_cut;
 
        // num_removed_tracks = 0;
        // Bool_t continue_removal = true;
 
        // while (continue_removal && KFParticles_PV.size() > 2) {
        //     Int_t worst_track_idx = -1;
        //     Double_t max_chi2 = 0.0;
        //     
        //     // Find track with maximum chi2 to current vertex
        //     for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
        //         const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
        //         if (chi2 > max_chi2) {
        //             worst_track_idx = i;
        //             max_chi2 = chi2;
        //         }
        //     }
        //     
        //     // Check if we should remove this track
        //     if (worst_track_idx >= 0 && max_chi2 > optimal_cut) {
        //         // Remove track contribution from vertex using SubtractFromVertex
        //         KFParticles_PV[worst_track_idx].SubtractFromVertex(primVtx);
        //         
        //         // Remove the track from the collection
        //         KFParticles_PV.erase(KFParticles_PV.begin() + worst_track_idx);
        //         
        //         ++num_removed_tracks;
        //     } else {
        //         continue_removal = false;
        //     }
        // }

        // Now primVtx contains the updated vertex after removing bad tracks
        // KFParticles_PV contains only the tracks that passed the chi2 cut

		//=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (CUT ANALYSIS) VAR-1
        //=========================================================================================================================
        // UNCOMMENT THIS PART FOR CUT ANALYSIS

        // // Save the original state
        // std::vector<KFParticle> original_KFParticles_PV = KFParticles_PV;
        // KFParticle original_primVtx = primVtx;
        //  
        // // Create vector of indices for efficient tracking
        // std::vector<Int_t> track_indices(original_KFParticles_PV.size());
        // for (Int_t i = 0; i < original_KFParticles_PV.size(); ++i) {
        //     track_indices[i] = i;
        // }
 
        // std::vector<Double_t> max_chi2_cut_vector = {30., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2.};
 
        // for (Double_t current_cut : max_chi2_cut_vector) {
        //     // Start with all tracks
        //     std::vector<Int_t> current_indices = track_indices;
        //     primVtx = original_primVtx; // Start with original vertex
 
        //     max_chi2_cut = current_cut;
        //     num_removed_tracks = 0;
        //     Bool_t continue_removal = true;
 
        //     while (continue_removal && current_indices.size() > 2) {
        //         Int_t worst_idx_pos = -1; // Position in current_indices
        //         Int_t worst_original_idx = -1; // Original track index
        //         max_chi2 = 0.0;
 
        //         // Find worst track
        //         for (Int_t i = 0; i < current_indices.size(); ++i) {
        //             Int_t track_idx = current_indices[i];
        //             const Double_t chi2 = original_KFParticles_PV[track_idx].GetDeviationFromVertex(primVtx);
        //             if (chi2 > max_chi2) {
        //                 worst_idx_pos = i;
        //                 worst_original_idx = track_idx;
        //                 max_chi2 = chi2;
        //             }
        //         }
 
        //         // Check cut
        //         if (worst_idx_pos >= 0 && max_chi2 > current_cut) {
        //             // Remove track index
        //             current_indices.erase(current_indices.begin() + worst_idx_pos);
        //     
        //             // Rebuild vertex
        //             primVtx = KFParticle(); // Reset
        //             for (Int_t idx : current_indices) {
        //                 primVtx += original_KFParticles_PV[idx];
        //             }
 
        //             ++num_removed_tracks;
        //         } else {
        //             continue_removal = false;
        //         }
        //     }
 
        //     // Store results
        //     PV_diff_ES_x = primVtx.GetX() - truePVPosition.X();
        //     PV_diff_ES_y = primVtx.GetY() - truePVPosition.Y();
        //     PV_diff_ES_z = primVtx.GetZ() - truePVPosition.Z();
        //     
        //     chi2_PV = primVtx.GetChi2();
        //     ndf_PV = primVtx.GetNDF();
        //     chi2_per_ndf_PV = (ndf_PV > 0) ? chi2_PV / ndf_PV : -1.0;
        //     tracks_per_vertex = current_indices.size();
        //     max_chi2_observed = max_chi2;
        //     
        //     tree->Fill();
        // }
 
        // Restore for any code after the loop
        // KFParticles_PV = original_KFParticles_PV;
        // primVtx = original_primVtx;

        //=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (CUT ANALYSIS) VAR-2
        //=========================================================================================================================

        // // Save the original state
        std::vector<KFParticle> original_KFParticles_PV = KFParticles_PV;
        KFParticle original_primVtx = primVtx;

        std::vector<Double_t> max_chi2_cut_vector = {30., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2.};

        for (Double_t current_cut : max_chi2_cut_vector) {
            // Restore original state for this iteration
            KFParticles_PV = original_KFParticles_PV;
            primVtx = original_primVtx;

            max_chi2_cut = current_cut;
            num_removed_tracks = 0;
            Bool_t continue_removal = true;
            Int_t min_tracks_for_vertex = 2; // Need at least 2 tracks

            while (continue_removal && KFParticles_PV.size() > min_tracks_for_vertex) {
                Int_t worst_track_idx = -1;
                max_chi2 = 0.0;

                // Find track with maximum chi2 to current vertex
                for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
                    const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
                    if (chi2 > max_chi2) {
                        worst_track_idx = i;
                        max_chi2 = chi2;
                    }
                }

                // Check if we should remove this track
                if (worst_track_idx >= 0 && max_chi2 > current_cut) {
                    // Remove track contribution using SubtractFromVertex
                    KFParticles_PV[worst_track_idx].SubtractFromVertex(primVtx);
                    
                    // Remove track from collection
                    KFParticles_PV.erase(KFParticles_PV.begin() + worst_track_idx);
                    
                    ++num_removed_tracks;
                } else {
                    continue_removal = false;
                }
            }

            // Store results
            PV_diff_ES_x = primVtx.GetX() - truePVPosition.X();
            PV_diff_ES_y = primVtx.GetY() - truePVPosition.Y();
            PV_diff_ES_z = primVtx.GetZ() - truePVPosition.Z();
            
            chi2_PV = primVtx.GetChi2();
            ndf_PV = primVtx.GetNDF();
            chi2_per_ndf_PV = (ndf_PV > 0) ? chi2_PV / ndf_PV : -1.0;
            tracks_per_vertex = KFParticles_PV.size();
            max_chi2_observed = max_chi2;
 
            tree->Fill();
        }

        // Restore for any code after the loop
        // KFParticles_PV = original_KFParticles_PV;
        // primVtx = original_primVtx;

        //=========================================================================================================================

        const KFParticle pr_vtx = primVtx;
        ++PV_extra_sel_count;
    }
	//===============================================================================================

    std::ofstream file_info("info_analysis_PV.txt");
	if (file_info.is_open())
	{
		file_info << PV_count << std::endl;
        file_info << PV_extra_sel_count << std::endl;
		file_info.close();
	}

	file -> Write();
	file -> Close();
}				

//======================================================================================================================
// Launcher function
//======================================================================================================================
void pv_dataset_27(Int_t SEED = 0 ) 
{   
	// SpdMCDataIterator* IT = 0;
	analyse(20000, SEED); 
}
