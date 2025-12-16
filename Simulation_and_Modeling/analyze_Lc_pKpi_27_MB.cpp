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
        // Check for negative probabilities
        if (probabilities[i] < 0.0) {
            return 0;
        }
    }

    // Apply threshold mask
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
// Selecting tracks that are used for secondary vertex reconstruction 
//======================================================================================================================
std::vector<KFParticle> KFparticles; 					// vector with selected particles 
std::vector<Int_t> pdgs;
std::vector <SpdMCParticle*> mc_particles; 				// vector with mc particles
std::vector<Int_t> Daughters = { 2212, -321, 211 }; 			// Daughters of Lc+(4122) decay 

void SelectTracks(KFParticle PV)
{
	KFparticles.clear();    
	pdgs.clear();
	mc_particles.clear();

	for (Int_t i{}; i < mctracks -> GetEntriesFast(); ++i) 
   	{
		SpdTrackMC *track{ (SpdTrackMC*) mctracks -> At(i) };
		if ( !track ) continue;											// Sel
		SpdMCParticle* particle{ (SpdMCParticle*) mcparticles -> At( track -> GetParticleId() ) };	
		Int_t charge{ track -> GetParticlePdg() / abs(track -> GetParticlePdg()) };	
		SpdTrackFitPar *track_pars{ track -> GetFitPars() };
		if ( !track_pars ) continue;										// Sel
		const SpdTrackState *state{ track_pars -> GetFirstState() };
		if ( !state ) continue;											// Sel
		TVector3 trkPos = state -> GetPosition();
		TVector3 trkMom = state -> GetMomentum();
		const TMatrixDSym& trkCov = state -> GetCovariance();
					       
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
		kfTrack.SetNDF(track_pars -> GetNDF());
		kfTrack.SetChi2(track_pars -> GetChi2());
		kfTrack.SetCharge(charge);

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
		if ( std::find( Daughters.begin(), Daughters.end(), pdgTrk ) == Daughters.end() ) continue;
		KFParticle p1(kfTrack, pdgTrk);        
		KFparticles.push_back(p1);
		pdgs.push_back(pdgTrk);
		mc_particles.push_back(particle);					
	}
}
//======================================================================================================================
// Checking if track fits to corresponding cuts
//======================================================================================================================
constexpr Int_t cutsNum{ 24 };
std::bitset<cutsNum> bitsetCutFlow{};			// Describes every event. Each bit corresponds to selection criteria. If the n-th bit on any event process step 
							// become bit == 1, then the current event don't fit to n-th cut ( numeration starts from 0).
void CheckTracks( SpdMCParticle* (&particles)[3], Int_t (&pdgs)[3] )
{
	std::bitset<3> bitFlagsCriteria_I{}; 				// Check first criteria(0 -- first track, 1 -- second track...) if bit == 1 then track was rejected
	std::bitset<3> bitFlagsCriteria_II{};				// Check second criteria(0 -- first track, 1 -- second track...) if bit == 1 then track was rejected
	std::bitset<3> bitFlagsCriteria_III{};				// Check third criteria(0 -- first track, 1 -- second track...) if bit == 1 then track was rejected
    std::bitset<3> bitFlagsCriteria_IV{};				// Check fourth criteria(0 -- first track, 1 -- second track...) if bit == 1 then track was rejected
	for ( Int_t j{}; j < 3; ++j ){
		SpdTrackMC* track = dynamic_cast<SpdTrackMC*>( mctracks->At( particles[j] -> GetTrackId() ) );
		if (track -> GetNHitsIts() < 3) bitFlagsCriteria_I.set(j);
		// if ( !track -> GetIsFitted() ) bitFlagsCriteria_II.set(j);
		// SpdTrackFitPar *track_pars{ track -> GetFitPars() };
		// if ( !track_pars ) continue;
		// if ( track_pars -> GetConvergency() == 0 ) bitFlagsCriteria_III.set(j);	// 0 (not converged), -1 (partially converged), 1(fully converged).
		// if ( track_pars -> GetChi2overNDF() > 12 ) bitFlagsCriteria_IV.set(j);
	}
	//If all the elements of bitFlagsCriteria_i == 0 ( bitFlagsCriteria_i.none() == 1 ) then an event goes through the selection criteria_i.
	if ( !bitFlagsCriteria_I.none() ) bitsetCutFlow.set(0); 					// An event stuck at the 1st cut
	// if ( !bitFlagsCriteria_II.none() ) bitsetCutFlow.set(1); 					// An event stuck at the 2nd cut
	// if ( !bitFlagsCriteria_III.none() ) bitsetCutFlow.set(2);					// An event stuck at the 3rd cut 
	// if ( !bitFlagsCriteria_IV.none() ) bitsetCutFlow.set(3);					// An event stuck at the 4th cut
}

//======================================================================================================================
// Main analysis
//======================================================================================================================
void sv_reconstruction(Int_t N, std::string_view inputFile, std::string outputFile) // N - max event number to analyse
{
	SpdMCDataIterator* IT = new SpdMCDataIterator();

	IT -> AddSourceFile( inputFile );

	IT -> ActivateBranch("all");
	IT -> Init();

    mcparticles = IT -> GetParticles();    
    mctracks    = IT -> GetTracks();         
    rcvertices  = IT -> GetVerticesRC();
    mcvertices  = IT -> GetVertices();    
    particles_tof = IT -> GetTofParticles();
    particles_ts = IT -> GetTsParticles();
	
	//======================================================================================================================
	// Process variables
	const Int_t events_max{ N };        // Upper limit for events to process
	Int_t n_event{};                    // Current event number
	Int_t Total_count{};                // Number of events under the Lc+(4212) peak (2.24763; 2.32497)GeV for Bg samples and only MC Truth real Lc+ for Sig samples.
	Int_t CutFlow_total[cutsNum] = {};  // Number of events under the Lc+(4212) peak (2.24763; 2.32497)GeV for Bg samples and only MC Truth real Lc+ for Sig samples.
                                        // Increasing i-th position if an event reaches i + 1 cuts step(passes through i + 1 cuts).
	// Other Variables
	Double_t mass_Lc{};
	
	Double_t P_p{};
	Double_t P_pip{};
	Double_t P_K{};
	Double_t P_Lc{};

	Double_t eta_p{};
	Double_t eta_pip{};
	Double_t eta_K{};
	Double_t eta_Lc{};

	Double_t Pt_Lc{};
	Double_t Pt_p{};
	Double_t Pt_K{};
	Double_t Pt_pip{};

	Double_t lengthXY_Lc{};
	Double_t dlengthXY_Lc{};
	Double_t ctau_Lc{};

	Double_t OA_p{};
	Double_t OA_K{};
	Double_t OA_pip{};
	Double_t ptOverE{};

	Double_t chi2_Lc_PV_xy{};
	Double_t dist_Lc_PV_xy{};
	Double_t dist_Lc_PV_xy_custom{};

	Double_t chi2_p_PV_xy{};
	Double_t dist_p_PV_xy{};
	Double_t dist_p_PV_xy_custom{};

	Double_t chi2_K_PV_xy{};
	Double_t dist_K_PV_xy{};
	Double_t dist_K_PV_xy_custom{};

	Double_t chi2_pip_PV_xy{};
	Double_t dist_pip_PV_xy{};
	Double_t dist_pip_PV_xy_custom{};
	
	Double_t chi2_p_Lc_xy{};
	Double_t dist_p_Lc_xy{};
	Double_t dist_p_Lc_xy_custom{};

	Double_t chi2_K_Lc_xy{};
	Double_t dist_K_Lc_xy{};
	Double_t dist_K_Lc_xy_custom{};

	Double_t chi2_pip_Lc_xy{};
	Double_t dist_pip_Lc_xy{};
	Double_t dist_pip_Lc_xy_custom{};
	
	Double_t chi2_Lc{};

	Double_t chi2_K_pip_xy{};
	Double_t dist_K_pip_xy{};
	Double_t dist_K_pip_xy_custom{};

	Double_t chi2_p_K_xy{};
	Double_t dist_p_K_xy{};
	Double_t dist_p_K_xy_custom{};

	Double_t chi2_p_pip_xy{};
	Double_t dist_p_pip_xy{};
	Double_t dist_p_pip_xy_custom{};

	Double_t cosAngle_r_Lc_momentum_Lc_xy{};
	Double_t cosAngle_r_Lc_sum_momentum_xy{};
	Double_t cosAngle_momentum_Lc_sum_momentum_xy{};

	Double_t xF{};
	Double_t phi_angle{};
	
	Int_t multiplicity{};

	// GetIsFittedOk?
	// Int_t IsFitted_p{};
	// Int_t IsGood_p{};
	// Int_t Convergency_p{};
	// Double_t Chi2overNDF_p{};
	// Int_t NHitsIts_p{};
	// Int_t NHitsTsB_p{};
	// Int_t NHitsTsEC_p{};

	// Int_t IsFitted_K{};
	// Int_t IsGood_K{};
	// Int_t Convergency_K{};
	// Double_t Chi2overNDF_K{};
	// Int_t NHitsIts_K{};
	// Int_t NHitsTsB_K{};
	// Int_t NHitsTsEC_K{};

	// Int_t IsFitted_pip{};
	// Int_t IsGood_pip{};
	// Int_t Convergency_pip{};
	// Double_t Chi2overNDF_pip{};
	// Int_t NHitsIts_pip{};
	// Int_t NHitsTsB_pip{};
	// Int_t NHitsTsEC_pip{};

	Double_t Lc_diff_x{};
	Double_t Lc_diff_y{};
	Double_t Lc_diff_z{};

	Int_t kf_pv_size_before_re_fit{};	// Amount of tracks for PV reco before re-fit loop
	Int_t kf_pv_size_after_re_fit{};	// Amount of tracks for PV reco after re-fit loop

	Double_t PV_diff_x{};
	Double_t PV_diff_y{};
	Double_t PV_diff_z{};

	Double_t PV_diff_ES_x{};
	Double_t PV_diff_ES_y{};
	Double_t PV_diff_ES_z{};

	Int_t true_decay{};
	//======================================================================================================================
	// Output File
	TFile* file = new TFile(outputFile.c_str(), "RECREATE");
    //======================================================================================================================
	// TTree
	TTree *tree = new TTree("tree", "tree");
	
	tree -> Branch("n_event", &n_event, "n_event/I");

	tree -> Branch("mass_Lc", &mass_Lc, "#Lambda_c^{+} mass/D");
			
	tree -> Branch("P_p", &P_p, "P momentum/D");
	tree -> Branch("P_pip", &P_pip, "#Pi+ momentum/D");
	tree -> Branch("P_K", &P_K, "K momentum/D");
	tree -> Branch("P_Lc", &P_Lc, "#Lambda_c momentum/D");

	tree -> Branch("eta_p", &eta_p, "P momentum/D");
	tree -> Branch("eta_pip", &eta_pip, "#Pi+ momentum/D");
	tree -> Branch("eta_K", &eta_p, "K momentum/D");
	tree -> Branch("eta_Lc", &eta_Lc, "#eta_{#Lambda_c}/D");

	tree -> Branch("Pt_Lc", &Pt_Lc, "P momentum/D");
	tree -> Branch("Pt_p", &Pt_p, "#Pi+ momentum/D");
	tree -> Branch("Pt_K", &Pt_K, "K momentum/D");
	tree -> Branch("Pt_pip", &Pt_pip, "#Lambda_c momentum/D");

	tree -> Branch("lengthXY_Lc", &lengthXY_Lc, "K momentum/D");
	tree -> Branch("dlengthXY_Lc", &dlengthXY_Lc, "#Lambda_c momentum/D");
	tree -> Branch("ctau_Lc", &ctau_Lc, "#Lambda_c momentum/D");

	tree -> Branch("OA_p", &OA_p, "P momentum/D");
	tree -> Branch("OA_K", &OA_K, "#Pi+ momentum/D");
	tree -> Branch("OA_pip", &OA_pip, "K momentum/D");
	tree -> Branch("ptOverE", &ptOverE, "#Lambda_c momentum/D");

	tree -> Branch("chi2_Lc_PV_xy", &chi2_Lc_PV_xy, "chi2_Lc_PV_xy/D");
    tree -> Branch("dist_Lc_PV_xy", &dist_Lc_PV_xy, "dist_Lc_PV_xy/D");
	tree -> Branch("dist_Lc_PV_xy_custom", &dist_Lc_PV_xy_custom, "dist_Lc_PV_xy_custom/D");

    tree -> Branch("chi2_p_PV_xy", &chi2_p_PV_xy, "chi2_p_PV_xy/D");
    tree -> Branch("dist_p_PV_xy", &dist_p_PV_xy, "dist_p_PV_xy/D");
	tree -> Branch("dist_p_PV_xy_custom", &dist_p_PV_xy_custom, "dist_p_PV_xy_custom/D");

    tree -> Branch("chi2_K_PV_xy", &chi2_K_PV_xy, "chi2_K_PV_xy/D");
    tree -> Branch("dist_K_PV_xy", &dist_K_PV_xy, "dist_K_PV_xy/D");
	tree -> Branch("dist_K_PV_xy_custom", &dist_K_PV_xy_custom, "dist_K_PV_xy_custom/D");

    tree -> Branch("chi2_pip_PV_xy", &chi2_pip_PV_xy, "chi2_pip_PV_xy/D");
    tree -> Branch("dist_pip_PV_xy", &dist_pip_PV_xy, "dist_pip_PV_xy/D");
	tree -> Branch("dist_pip_PV_xy_custom", &dist_pip_PV_xy_custom, "dist_pip_PV_xy_custom/D");

    tree -> Branch("chi2_p_Lc_xy", &chi2_p_Lc_xy, "chi2_p_Lc_xy/D");
    tree -> Branch("dist_p_Lc_xy", &dist_p_Lc_xy, "dist_p_Lc_xy/D");
	tree -> Branch("dist_p_Lc_xy_custom", &dist_p_Lc_xy_custom, "dist_p_Lc_xy_custom/D");

    tree -> Branch("chi2_K_Lc_xy", &chi2_K_Lc_xy, "chi2_K_Lc_xy/D");
    tree -> Branch("dist_K_Lc_xy", &dist_K_Lc_xy, "dist_K_Lc_xy/D");
	tree -> Branch("dist_K_Lc_xy_custom", &dist_K_Lc_xy_custom, "dist_K_Lc_xy_custom/D");

    tree -> Branch("chi2_pip_Lc_xy", &chi2_pip_Lc_xy, "chi2_pip_Lc_xy/D");
    tree -> Branch("dist_pip_Lc_xy", &dist_pip_Lc_xy, "dist_pip_Lc_xy/D");
	tree -> Branch("dist_pip_Lc_xy_custom", &dist_pip_Lc_xy_custom, "dist_pip_Lc_xy_custom/D");

    tree -> Branch("chi2_Lc", &chi2_Lc, "chi2_Lc/D");

    tree -> Branch("chi2_K_pip_xy", &chi2_K_pip_xy, "chi2_K_pip_xy/D");
    tree -> Branch("dist_K_pip_xy", &dist_K_pip_xy, "dist_K_pip_xy/D");
	tree -> Branch("dist_K_pip_xy_custom", &dist_K_pip_xy_custom, "dist_K_pip_xy_custom/D");

    tree -> Branch("chi2_p_K_xy", &chi2_p_K_xy, "chi2_p_K_xy/D");
    tree -> Branch("dist_p_K_xy", &dist_p_K_xy, "dist_p_K_xy/D");
	tree -> Branch("dist_p_K_xy_custom", &dist_p_K_xy_custom, "dist_p_K_xy_custom/D");

    tree -> Branch("chi2_p_pip_xy", &chi2_p_pip_xy, "chi2_p_pip_xy/D");
    tree -> Branch("dist_p_pip_xy", &dist_p_pip_xy, "dist_p_pip_xy/D");
	tree -> Branch("dist_p_pip_xy_custom", &dist_p_pip_xy_custom, "dist_p_pip_xy_custom/D");

	tree -> Branch("cosAngle_r_Lc_momentum_Lc_xy", &cosAngle_r_Lc_momentum_Lc_xy, "#Lambda_c momentum_xy/D");
	tree -> Branch("cosAngle_r_Lc_sum_momentum_xy", &cosAngle_r_Lc_sum_momentum_xy, "#Lambda_c momentum_xy/D");
	tree -> Branch("cosAngle_momentum_Lc_sum_momentum_xy", &cosAngle_momentum_Lc_sum_momentum_xy, "#Lambda_c momentum_xy/D");
	
	tree -> Branch("xF", &xF, "#Lambda_c^{+} xF/D");
	tree -> Branch("phi_angle", &phi_angle, "#Lambda_c^{+} #Phi angle/D");
	
	tree -> Branch("multiplicity", &multiplicity, "Multiplicity/I");

	// tree -> Branch("NHits_pi", &NHits_pi, "Multiplicity/I");
	// tree -> Branch("NHitsIts_pi", &NHitsIts_pi, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_pi", &NHitsTsB_pi, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_pi", &NHitsTsEC_pi, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_pi", &Chi2OverNDF_pi, "Multiplicity/I");

	// tree -> Branch("NHits_pi_true", &NHits_pi_true, "Multiplicity/I");
	// tree -> Branch("NHitsIts_pi_true", &NHitsIts_pi_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_pi_true", &NHitsTsB_pi_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_pi_true", &NHitsTsEC_pi_true, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_pi_true", &Chi2OverNDF_pi_true, "Multiplicity/I");

	// tree -> Branch("NHits_p", &NHits_p, "Multiplicity/I");
	// tree -> Branch("NHitsIts_p", &NHitsIts_p, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_p", &NHitsTsB_p, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_p", &NHitsTsEC_p, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_p", &Chi2OverNDF_p, "Multiplicity/I");

	// tree -> Branch("NHits_p_true", &NHits_p_true, "Multiplicity/I");
	// tree -> Branch("NHitsIts_p_true", &NHitsIts_p_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_p_true", &NHitsTsB_p_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_p_true", &NHitsTsEC_p_true, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_p_true", &Chi2OverNDF_p_true, "Multiplicity/I");

	// tree -> Branch("NHits_K", &NHits_K, "Multiplicity/I");
	// tree -> Branch("NHitsIts_K", &NHitsIts_K, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_K", &NHitsTsB_K, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_K", &NHitsTsEC_K, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_K", &Chi2OverNDF_K, "Multiplicity/I");

	// tree -> Branch("NHits_K_true", &NHits_K_true, "Multiplicity/I");
	// tree -> Branch("NHitsIts_K_true", &NHitsIts_K_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsB_K_true", &NHitsTsB_K_true, "Multiplicity/I");
	// tree -> Branch("NHitsTsEC_K_true", &NHitsTsEC_K_true, "Multiplicity/I");
	// tree -> Branch("Chi2OverNDF_K_true", &Chi2OverNDF_K_true, "Multiplicity/I");

	tree -> Branch("Lc_diff_x", &Lc_diff_x, "Lc_diff_x/D");
	tree -> Branch("Lc_diff_y", &Lc_diff_y, "Lc_diff_y/D");
	tree -> Branch("Lc_diff_z", &Lc_diff_z, "Lc_diff_z/D");

	tree->Branch("kf_pv_size_before_re_fit", &kf_pv_size_before_re_fit, "kf_pv_size_before_re_fit/I");
	tree->Branch("kf_pv_size_after_re_fit", &kf_pv_size_after_re_fit, "kf_pv_size_after_re_fit/I");

	tree -> Branch("PV_diff_x", &PV_diff_x, "PV_diff_x/D");
	tree -> Branch("PV_diff_y", &PV_diff_y, "PV_diff_y/D");
	tree -> Branch("PV_diff_z", &PV_diff_z, "PV_diff_z/D");

	tree -> Branch("PV_diff_ES_x", &PV_diff_ES_x, "PV_diff_ES_x/D");
	tree -> Branch("PV_diff_ES_y", &PV_diff_ES_y, "PV_diff_ES_y/D");
	tree -> Branch("PV_diff_ES_z", &PV_diff_ES_z, "PV_diff_ES_z/D");

	tree -> Branch("true_decay", &true_decay, "#Lambda_c^{+} true decay/I");

	//======================================================================================================================
    while ( IT -> NextEvent() && n_event < events_max ) 
	{
		++n_event;		
		multiplicity = mctracks -> GetEntriesFast();
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
		PV_pos_prefit = prim_vtx_fit -> GetVertex();        // PV RC position first approximation.  
		SelectTracksForPV( PV_pos_prefit );                 // Tracks selection for PV reconstruction
		if ( !findPVCheck )     continue;                   // Skip if no PV.
		
		// if ( KFParticles_PV.size() < 5 ) continue;		// Check if have more than 5 tracks for PV reconstruction.
		kf_pv_size_before_re_fit = KFParticles_PV.size();

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
		KFParticle primVtx;
		Int_t nDaughters_temp = KFParticles_PV.size();
		const KFParticle* vDaughters_temp[nDaughters_temp];
		for ( Int_t i{}; i < KFParticles_PV.size(); ++i ) vDaughters_temp[i] = &KFParticles_PV[i];

		primVtx.Construct(vDaughters_temp,nDaughters_temp, 0, -1); 				// 0 --- no parant particle, 1 --- no mass hypothesis

		PV_diff_x = primVtx.GetX() - truePVPosition.X();
		PV_diff_y = primVtx.GetY() - truePVPosition.Y();
		PV_diff_z = primVtx.GetZ() - truePVPosition.Z();
		//=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (MAIN ANALYSIS) VAR-1
        //=========================================================================================================================

        Double_t optimal_cut = 5.0;
 
        Bool_t continue_removal = true;
 
        while (continue_removal && KFParticles_PV.size() > 2) {
            Int_t worst_track_idx = -1;
            Double_t max_chi2 = 0.0;
            
            // Find track with maximum chi2 to current vertex
            for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
                const Double_t chi2 = KFParticles_PV[i].GetDeviationFromVertex(primVtx);
                if (chi2 > max_chi2) {
                    worst_track_idx = i;
                    max_chi2 = chi2;
                }
            }
            
            // Check if we should remove this track
            if (worst_track_idx >= 0 && max_chi2 > optimal_cut) {
                // Remove the worst track
                KFParticles_PV.erase(KFParticles_PV.begin() + worst_track_idx);
                
                // Recalculate vertex from remaining tracks
                primVtx = KFParticle(); // Reset vertex
                for (Int_t i = 0; i < KFParticles_PV.size(); ++i) {
					primVtx += KFParticles_PV[i]; // Add all remaining tracks to rebuild vertex
				}
                
            } else {
                continue_removal = false;
            }
        }

        // Now primVtx contains the updated vertex after removing bad tracks
        // KFParticles_PV contains only the tracks that passed the chi2 cut

		// if ( KFParticles_PV.size() < 5 ) continue;		 // Min number of tracks to PV reco after re-fit loop.
		kf_pv_size_after_re_fit = KFParticles_PV.size();

		PV_diff_ES_x = primVtx.GetX() - truePVPosition.X();
		PV_diff_ES_y = primVtx.GetY() - truePVPosition.Y();
		PV_diff_ES_z = primVtx.GetZ() - truePVPosition.Z();
		
		const KFParticle pr_vtx = primVtx;            					// init primary vertex
		//=========================================================================================================================
		// Futher analysis
		//=========================================================================================================================
		// Selecting tracks that are used for secondary vertex reconstruction
		SelectTracks(pr_vtx);
		//=========================================================================================================================
		// Loop over tracks and fill the histograms
		Double_t sumEnergy{};	// Sum of energy of all the selected tracks in event with theta > 0.174533
		for (Int_t i{}; i < mc_particles.size(); ++i ) {
            // SpdTrackMC* track = dynamic_cast<SpdTrackMC*>( mctracks -> At( mc_particles[i] -> GetTrackId() ) );
			// if ( !track ) continue; // Maybe skip?
			// SpdTrackFitPar *track_pars{ track -> GetFitPars() };
			// if ( !track_pars ) continue; // Maybe skip?

			KFParticle particle{KFparticles[i]};	
	
			const Double_t energy{particle.GetE()};	
			const Double_t theta{particle.GetTheta()};
			if ( theta > 0.174533 ) sumEnergy += energy;		        
        }
		//=======================================================================================================
		// Skip event if number of tracks selected < number of tracks needed for decay
		if (KFparticles.size() < Daughters.size()) continue;  						
		//=======================================================================================================
		SpdVertexCombiFinder* vfinder = new SpdVertexCombiFinder();
		vfinder -> InitParticles(pdgs);              			// Init list of particles
		vfinder -> InitVertex(Daughters);           	 		// Init decay particles in SV 

		SpdMCParticle* d_particles[3];		
		Int_t pdg[3];

		std::vector <Int_t> vc;

		while ( vfinder -> Next(vc) ) 					// Each iteration takes new particles(tracks) combination (of 3 particles for the current decay)
		{
			bitsetCutFlow.reset();
	
			Int_t up_iter{ static_cast<Int_t>( vc.size() ) };
			Int_t i_p{}, i_pip{}, i_K{};  				// Containing indices for proton and pion+- in the three tracks combination.
			for (Int_t i{}; i < up_iter; ++i) 
			{ 
				d_particles[i] = mc_particles[vc[i]];
				pdg[i] = pdgs[vc[i]];
				if ( pdgs[vc[i]] == 2212 ) i_p = i;
				if ( pdgs[vc[i]] == -321 ) i_K = i;
				if ( pdgs[vc[i]] == 211 ) i_pip = i;	
			} 
			//============================================================================================
			// Tracks selection
			CheckTracks( d_particles, pdg ); 

			// for ( Int_t j{}; j < 3; ++j ){
			// 	SpdTrackMC *track = dynamic_cast<SpdTrackMC*>( mctracks->At( d_particles[j] -> GetTrackId() ) );
            //     SpdTrackFitPar *track_pars{ track -> GetFitPars() };
            //     if ( !track_pars ) continue;

            //     if ( j == i_p ){
            //         IsFitted_p = track -> GetIsFitted();
            //         IsGood_p = track_pars -> GetIsGood();
            //         Convergency_p = track_pars -> GetConvergency();
            //         Chi2overNDF_p = track_pars -> GetChi2overNDF();
            //         NHitsIts_p = track -> GetNHitsIts();
            //         NHitsTsB_p = track -> GetNHitsTsB();
            //         NHitsTsEC_p = track -> GetNHitsTsEC();
            //     }

            //     if ( j == i_K ){
            //         IsFitted_K = track -> GetIsFitted();
            //         IsGood_K = track_pars -> GetIsGood();
            //         Convergency_K = track_pars -> GetConvergency();
            //         Chi2overNDF_K = track_pars -> GetChi2overNDF();
            //         NHitsIts_K = track -> GetNHitsIts();
            //         NHitsTsB_K = track -> GetNHitsTsB();
            //         NHitsTsEC_K = track -> GetNHitsTsEC();
            //     }

            //     if ( j == i_pip ){
            //         IsFitted_pip = track -> GetIsFitted();
            //         IsGood_pip = track_pars -> GetIsGood();
            //         Convergency_pip = track_pars -> GetConvergency();
            //         Chi2overNDF_pip = track_pars -> GetChi2overNDF();
            //         NHitsIts_pip = track -> GetNHitsIts();
            //         NHitsTsB_pip = track -> GetNHitsTsB();
            //         NHitsTsEC_pip = track -> GetNHitsTsEC();
            //     }
			// }
			//============================================================================================
			// KFParticles reconstruction
			KFParticle p = KFparticles[vc[i_p]];
			KFParticle pip = KFparticles[vc[i_pip]];
			KFParticle K = KFparticles[vc[i_K]];
			KFParticle Lc( p, K, pip );

			Lc.SetProductionVertex(pr_vtx);
			K.SetProductionVertex(Lc);
			p.SetProductionVertex(Lc);
			pip.SetProductionVertex(Lc);
			//==============================================================================================
			mass_Lc = Lc.GetMass();

			const std::vector<Double_t> Mass = {2.18, 2.4};
			// const std::vector<Double_t> Mass = {1.6, 3.};
            if (mass_Lc < Mass[0] || mass_Lc > Mass[1]) continue; 
			
			Double_t minMass_Lc{};
			for (Int_t i{}; i < up_iter; ++i)
			{	
				const Double_t m{ KFparticles[vc[i]].GetMass() };
				if ( m > minMass_Lc ) minMass_Lc = m;
			}

			if (mass_Lc <= minMass_Lc) continue;

			p.TransportToProductionVertex();
			pip.TransportToProductionVertex();
			K.TransportToProductionVertex();
			Lc.TransportToProductionVertex();
			
			P_p = p.GetMomentum();
			P_pip = pip.GetMomentum();
			P_K = K.GetMomentum();
			P_Lc = Lc.GetMomentum();

			eta_Lc = Lc.GetEta();
			eta_p = p.GetEta();
			eta_K = K.GetEta();
			eta_pip = pip.GetEta();

			Pt_Lc = Lc.GetPt();
			Pt_p = p.GetPt();
			Pt_K = K.GetPt();
			Pt_pip = pip.GetPt();

			lengthXY_Lc = Lc.GetDecayLengthXY();
			dlengthXY_Lc = Lc.GetErrDecayLengthXY();
			ctau_Lc = Lc.GetLifeTime();

			phi_angle = Lc.GetPhi();
		
			xF = 2 * Lc.GetPz() / 27.;   
			
			TLorentzVector P_kaon( K.GetPx(), K.GetPy(), K.GetPz(), K.GetE() );
			TLorentzVector P_proton( p.GetPx(), p.GetPy(), p.GetPz(), p.GetE() );
			TLorentzVector P_pion( pip.GetPx(), pip.GetPy(), pip.GetPz(), pip.GetE() );
			const TLorentzVector P4_Lc( Lc.GetPx(), Lc.GetPy(), Lc.GetPz(), Lc.GetE() );
			P_proton.Boost(-P4_Lc.BoostVector());
			P_kaon.Boost(-P4_Lc.BoostVector());
			P_pion.Boost(-P4_Lc.BoostVector());
			const TVector3 P_kaon_frame_Lc( P_kaon[0], P_kaon[1], P_kaon[2] );
			const TVector3 P_pion_frame_Lc( P_pion[0], P_pion[1], P_pion[2] );
			const TVector3 P_proton_frame_Lc( P_proton[0], P_proton[1], P_proton[2] );
			OA_p = P_proton_frame_Lc.Angle( -P4_Lc.BoostVector() );
			OA_K = P_kaon_frame_Lc.Angle( -P4_Lc.BoostVector() );
			OA_pip = P_pion_frame_Lc.Angle( -P4_Lc.BoostVector() );

			ptOverE = Lc.GetPt() / sumEnergy;
			
			//==================================================================================
			// True decay
			// Int_t true_decay{};	// Indicates if we reconstructed real Lc+(4212) decay
			true_decay = 0;
			SpdMCParticle* mc_m_proton = (SpdMCParticle*) mcparticles -> At( d_particles[i_p] -> GetMotherId() );
			SpdMCParticle* mc_m_pip = (SpdMCParticle*) mcparticles -> At( d_particles[i_pip] -> GetMotherId() );
			SpdMCParticle* mc_m_K = (SpdMCParticle*) mcparticles -> At( d_particles[i_K] -> GetMotherId() );
			if 
			(
				mc_m_proton -> GetPdgCode() == 4122 && mc_m_pip -> GetPdgCode() == 4122 && mc_m_K -> GetPdgCode() == 4122 &&
                d_particles[i_p]-> GetMotherId() == d_particles[i_pip]-> GetMotherId() &&
                d_particles[i_pip]-> GetMotherId() == d_particles[i_K]-> GetMotherId() &&
				pdg[i_p] == 2212 && pdg[i_K] == -321 && pdg[i_pip] == 211
			)
			{
				true_decay = 1;
			}
            if ( true_decay ) continue; // Only BG from MB
			//==================================================================================
			// Lc vertex fit check
			TVector3 Lc_production_position{Lc.GetX(), Lc.GetY(), Lc.GetZ()};
			Lc.TransportToDecayVertex();
            TVector3 PV_position{ pr_vtx.GetX(), pr_vtx.GetY(), pr_vtx.GetZ() };
			TVector3 Lc_position{Lc.GetX(), Lc.GetY(), Lc.GetZ()};
            TVector3 Lc_true_position{ d_particles[0] -> GetStartPos() }; // (MC_truth) Start point of Lc daughter track(NO guarantee that this is true Lc decay!!!).
			
			std::vector<TVector3> tracksPositions_Lc{}; // Tracks positions after extrapolation to Lc+(4122).
			std::vector<TVector3> tracksPositions_PV{}; // Tracks positions after extrapolation to PV.
			SpdTrackPropagatorGF Propagator_PV{};
			SpdTrackPropagatorGF Propagator_Lc{};
            Propagator_Lc.Init();
			Propagator_PV.Init();
                
			Bool_t check_vertexGood{ 1 };
		    for ( Int_t j{}; j < 3; ++j ){
                SpdTrackMC* track = dynamic_cast<SpdTrackMC*>( mctracks->At( d_particles[j] -> GetTrackId() ) );
       			Int_t charge{ track -> GetParticlePdg() / abs(track -> GetParticlePdg()) };
          		SpdTrackFitPar *tpars{ track -> GetFitPars() };
         		const SpdTrackState *state{ tpars -> GetFirstState() };
        		// Extrapolation to Lc+(4122)
				Propagator_Lc.InitTrack( pdg[j], 0 ); // propagation direction (-1, 0, 1) -> (backward, auto, forward) ???
				SpdTrackState  stateOut_Lc;
				const Double_t dist_Lc{ Propagator_Lc.ExtrapolateToPoint(Lc_position, *state, stateOut_Lc) };
                if (dist_Lc == 0.0) check_vertexGood = 0;
                TVector3 trkPos_Lc = stateOut_Lc.GetPosition();
                tracksPositions_Lc.push_back(trkPos_Lc);
				//==========================================================================================
                // Extrapolation to PV
				Propagator_PV.InitTrack( pdg[j], 0 ); // propagation direction (-1, 0, 1) -> (backward, auto, forward) ???
				SpdTrackState  stateOut_PV;
				const Double_t dist_PV{ Propagator_PV.ExtrapolateToPoint(PV_position, *state, stateOut_PV) };
                if (dist_PV == 0.0) check_vertexGood = 0;
                TVector3 trkPos_PV = stateOut_PV.GetPosition();
                tracksPositions_PV.push_back(trkPos_PV);	
			}
			if ( !check_vertexGood ) continue;

			// XY-plane
			Lc_production_position.SetZ(0.0);
			Lc_position.SetZ(0.0);
			PV_position.SetZ(0.0);

			for (auto& v : tracksPositions_PV) {
    			v.SetZ(0.0);
			}

			for (auto& v : tracksPositions_Lc) {
    			v.SetZ(0.0);
			}

			Lc.TransportToProductionVertex();
            chi2_Lc_PV_xy = Lc.GetDeviationFromVertexXY(pr_vtx);
			dist_Lc_PV_xy = Lc.GetDistanceFromVertexXY(pr_vtx);
			dist_Lc_PV_xy_custom = ( Lc_production_position - PV_position ).Mag();

            p.TransportToProductionVertex();
            chi2_p_PV_xy = p.GetDeviationFromVertexXY(pr_vtx);
			dist_p_PV_xy = p.GetDistanceFromVertexXY(pr_vtx);
			dist_p_PV_xy_custom = (PV_position - tracksPositions_PV[i_p]).Mag();

            K.TransportToProductionVertex();
            chi2_K_PV_xy = K.GetDeviationFromVertexXY(pr_vtx);
			dist_K_PV_xy = K.GetDistanceFromVertexXY(pr_vtx);
            dist_K_PV_xy_custom = (PV_position - tracksPositions_PV[i_K]).Mag();

            pip.TransportToProductionVertex();
            chi2_pip_PV_xy = pip.GetDeviationFromVertexXY(pr_vtx);
			dist_pip_PV_xy = pip.GetDistanceFromVertexXY(pr_vtx);
			dist_pip_PV_xy_custom = (PV_position - tracksPositions_PV[i_pip]).Mag();
			
            Lc.TransportToDecayVertex();
            chi2_p_Lc_xy = p.GetDeviationFromVertexXY(Lc);
			dist_p_Lc_xy = p.GetDistanceFromVertexXY(Lc);
			dist_p_Lc_xy_custom = (Lc_position - tracksPositions_Lc[i_p]).Mag();

            chi2_K_Lc_xy = K.GetDeviationFromVertexXY(Lc); 
			dist_K_Lc_xy = K.GetDistanceFromVertexXY(Lc);
			dist_K_Lc_xy_custom = (Lc_position - tracksPositions_Lc[i_K]).Mag();

            chi2_pip_Lc_xy = pip.GetDeviationFromVertexXY(Lc); 
			dist_pip_Lc_xy = pip.GetDistanceFromVertexXY(Lc);
			dist_pip_Lc_xy_custom = (Lc_position - tracksPositions_Lc[i_pip]).Mag();

            chi2_Lc = Lc.GetChi2();

            chi2_K_pip_xy = K.GetDeviationFromParticleXY(pip);
			dist_K_pip_xy = K.GetDistanceFromParticleXY(pip);
            dist_K_pip_xy_custom = ( tracksPositions_Lc[i_K] - tracksPositions_Lc[i_pip] ).Mag();

			chi2_p_K_xy = p.GetDeviationFromParticleXY(K);
			dist_p_K_xy = p.GetDistanceFromParticleXY(K);
            dist_p_K_xy_custom = ( tracksPositions_Lc[i_p] - tracksPositions_Lc[i_K] ).Mag();

			chi2_p_pip_xy = p.GetDeviationFromParticleXY(pip);
			dist_p_pip_xy = p.GetDistanceFromParticleXY(pip);
            dist_p_pip_xy_custom = ( tracksPositions_Lc[i_p] - tracksPositions_Lc[i_pip] ).Mag();
		
			//==================================================================================
			// Momentum direction check (XY - plate])
            Lc.TransportToProductionVertex();
			const TVector3 startPosLc{ Lc.GetX(), Lc.GetY(), 0 };
			Lc.TransportToDecayVertex();	
			const TVector3 lastPosLc{ Lc.GetX(), Lc.GetY(), 0 };
			Lc.TransportToProductionVertex();
			const TVector3 r_Lc{ lastPosLc - startPosLc };
			const TVector3 momentum_Lc{Lc.GetPx(), Lc.GetPy(), 0 };
			const TVector3 momentum_p{p.GetPx(), p.GetPy(), 0 };
			const TVector3 momentum_K{K.GetPx(), K.GetPy(), 0 };
			const TVector3 momentum_pip{pip.GetPx(), pip.GetPy(), 0 };
			const TVector3 sum_momentum{momentum_p + momentum_K + momentum_pip}; // Sum momentum of p, K, pi from Lc

			cosAngle_r_Lc_momentum_Lc = r_Lc.Angle(momentum_Lc);
			cosAngle_r_Lc_sum_momentum = r_Lc.Angle(sum_momentum);
			cosAngle_momentum_Lc_sum_momentum = momentum_Lc.Angle(sum_momentum);
			//==================================================================================
			// Selection
			//==================================================================================
			// If the n-th bit at any event process step
			// become a bit == 1, then the current event don't fit to n-th cut ( numeration starts from 0).
			// Example: bitsetCutFlow.set(7) ---> current event don't fit to 8th applied cut (starting from 1).

			// if (
			// maxChi2ToPV_d > 5. || minChi2ToPV_d > 2. || maxDistToPV_d > 0.4 || minDistToPV_d > 0.06 || maxChi2ToLc_d > 13. || minChi2ToLc_d > 3. ||
			// maxDistToLc_d > 0.3 || minDistToLc_d > 0.1 || maxChi2TrkToTrk > 2.5 || minChi2TrkToTrk > 0.15 || maxDistTrkToTrk > 0.5 || minDistTrkToTrk > 0.14 ||
			// Chi2ToPV_Lc > 2. || distLcToPV > 0.02 || Chi2LcVtx > 45.
			// )
			// { bitsetCutFlow.set(4); }
			//==================================================================================
			Lc_diff_x = Lc_position.X() - Lc_true_position.X();
			Lc_diff_y = Lc_position.Y() - Lc_true_position.Y();
			Lc_diff_z = Lc_position.Z() - Lc_true_position.Z();
			//==================================================================================
			// Check cuts + Cut Flow
            if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 )
                Total_count += 1;
            std::bitset<cutsNum> checkCuts{};   		// checkCuts.test(i) == 1 if an event fits well to the first i + 1 cuts.
            for ( Int_t i{}; i < cutsNum; ++i )
            {
            	Bool_t cutCheck{1};			// True if an event passes through all the cuts to the current step i. 
            						// If i == 1 then cutCheck is true only when the event passes through cut 1 & cut 2;
            	for ( Int_t j{}; j < i + 1; ++j )
            	{
            		if ( bitsetCutFlow.test(j) ) cutCheck = 0;
            	}
            	if ( cutCheck )				// If an event fits to all the cuts up to selection step i(i.o. i + 1 cuts) then increasing CutFlow_total.
            	{
            		checkCuts.set(i);
            		if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 )
            			CutFlow_total[i] += 1; 	// Increasing because an event reaches i + 1 cuts step(passes through i + 1 cuts)
            	}
            	else
            		break;
            }
			//==================================================================================
			// If you use selection(Check cuts + Cut Flow) in this script
			// DON'T FORGET TO USE HERE SOMETHING LIKE
			//if ( cutCheck ):
			//	tree -> fill...	

			if ( checkCuts.test(0) ){
				tree -> Fill();
			}
			//==================================================================================
			// Fill histograms
			// checkCuts.test(i) == 1 if an event fits well to the first i + 1 cuts.
		}
    } 
	//===============================================================================================
	std::ofstream file_info("info_analysis.txt");
	if (file_info.is_open())
	{
		file_info << n_event << std::endl;
		file_info.close();
	}
	std::ofstream file_CutFlow("CutFlow.txt");
	if ( file_CutFlow.is_open() )
	{
		file_CutFlow << Total_count << std::endl;
		for ( Int_t i{}; i < cutsNum; ++i ) file_CutFlow << CutFlow_total[i] << std::endl;
		file_CutFlow.close();
	}
	//===============================================================================================

	file -> Write();
	file -> Close();
}

//======================================================================================================================
// Launcher function
//======================================================================================================================
void analyze_Lc_pKpi_27_MB( std::string_view inputFile, std::string outputFile ) 
{   
	// SpdMCDataIterator* IT = 0;
	Int_t nMax{ 100000 };
	sv_reconstruction(nMax, inputFile, outputFile); 
}
