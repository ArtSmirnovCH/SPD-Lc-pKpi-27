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

		if ( !tof_pid ) tof_pid = 211;

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
void analyse(SpdMCDataIterator* IT, Int_t N, Int_t seed) // N - max event number to analyse
{
    IT = new SpdMCDataIterator();

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

    Double_t PV_diff_ES_x{};
    Double_t PV_diff_ES_y{};
    Double_t PV_diff_ES_z{};

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

    Double_t chi2_PV{};
    Int_t num_removed_tracks{};
    Int_t tracks_per_vertex{};
    Int_t PV_accept{};
    Double_t max_chi2{};
    Double_t max_chi2_cut{};

    Int_t PV_count{};
    Int_t PV_extra_sel_count{};

	//======================================================================================================================
	// Output File
	TFile* file;
    file = new TFile("ana_PV_10.root", "RECREATE");
	//======================================================================================================================
	// TTree
	TTree *tree = new TTree("tree", "tree");

    tree -> Branch("n_event", &n_event, "n_event/I");

    tree -> Branch("PV_diff_x", &PV_diff_x, "PV_diff_x/D");
    tree -> Branch("PV_diff_y", &PV_diff_y, "PV_diff_y/D");
    tree -> Branch("PV_diff_z", &PV_diff_z, "PV_diff_z/D");

    tree -> Branch("PV_diff_ES_x", &PV_diff_ES_x, "PV_diff_ES_x/D");
    tree -> Branch("PV_diff_ES_y", &PV_diff_ES_y, "PV_diff_ES_y/D");
    tree -> Branch("PV_diff_ES_z", &PV_diff_ES_z, "PV_diff_ES_z/D");

    tree -> Branch("chi2_PV", &chi2_PV, "chi2_PV/D");
    tree -> Branch("num_removed_tracks", &num_removed_tracks, "num_removed_tracks/I");
    tree -> Branch("tracks_per_vertex", &tracks_per_vertex, "tracks_per_vertex/I");
    tree -> Branch("max_chi2", &max_chi2, "max_chi2/D");
    tree -> Branch("max_chi2_cut", &max_chi2_cut, "max_chi2_cut/D");

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
		if ( KFParticles_PV.size() < 5 ) continue;						// Check if have more than 4 tracks fo PV reconstruction.
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
		// Re-reconstructing primary vertex with KFParticle

        max_chi2 = 0.0;
        for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
            const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
            if (chi2 > max_chi2) {
                max_chi2 = chi2;
            }
        }

        // Save the original state before the loop
        std::vector<KFParticle> original_KFParticles_PV = KFParticles_PV;
        KFParticle original_primVtx = primVtx;

        std::vector<Double_t> max_chi2_cut_vector = {9., 10.};
        for (Double_t current_cut : max_chi2_cut_vector) {
            // Restore original state for each iteration
            KFParticles_PV = original_KFParticles_PV;
            primVtx = original_primVtx;
            max_chi2_cut = current_cut;  // Store which cut we're using
            
            num_removed_tracks = 0;
            Bool_t flag = true;
            Int_t limitRemovedTracks = static_cast<Int_t>(original_KFParticles_PV.size()) - 1;
            
            while (flag) {
                Int_t badTrk = -1;
                Double_t maxChi2 = 0.0;
                
                // Find track with maximum chi2
                for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
                    const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
                    if (chi2 > maxChi2) {
                        badTrk = i;
                        maxChi2 = chi2;
                    }
                }
                
                // Check if we should remove this track
                if (maxChi2 > current_cut && num_removed_tracks < limitRemovedTracks && badTrk >= 0) {
                    KFParticles_PV[badTrk].SubtractFromVertex(primVtx);
                    KFParticles_PV.erase(KFParticles_PV.begin() + badTrk);
                    ++num_removed_tracks;
                } else {
                    flag = false;
                }
            }
            
            // Calculate differences with final vertex
            PV_diff_ES_x = primVtx.GetX() - truePVPosition.X();
            PV_diff_ES_y = primVtx.GetY() - truePVPosition.Y();
            PV_diff_ES_z = primVtx.GetZ() - truePVPosition.Z();
            
            chi2_PV = primVtx.GetChi2();
            tracks_per_vertex = KFParticles_PV.size();
            
            tree->Fill();
        }

        // Restore for any code after the loop
        KFParticles_PV = original_KFParticles_PV;
        primVtx = original_primVtx;

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
	SpdMCDataIterator* IT = 0;
	analyse(IT, 20000, SEED); 
}
