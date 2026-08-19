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
    
    double max_val = *std::max_element(array.begin(), array.end());
    
    std::vector<double> exp_vals(array.size());
    std::transform(array.begin(), array.end(), exp_vals.begin(),
                   [max_val](double x) { return std::exp(x - max_val); });
    

    double sum_exp = std::accumulate(exp_vals.begin(), exp_vals.end(), 0.0);
    
    // Normalize by the sum
    std::vector<double> result(exp_vals.size());
    std::transform(exp_vals.begin(), exp_vals.end(), result.begin(),
                   [sum_exp](double exp_val) { return exp_val / sum_exp; });
    
    return result;
}


int choice_with_thresholds(const std::vector<double>& probabilities,
                          const std::vector<double>& thresholds = {}) {

    std::vector<double> actual_thresholds = thresholds;
    if (actual_thresholds.empty()) {
        actual_thresholds = std::vector<double>(probabilities.size(), 0.0);
    }

    if (probabilities.size() != actual_thresholds.size()) {
        return 0;
    }

    if (probabilities.empty()) {
        return 0;
    }

    for (size_t i = 0; i < probabilities.size(); ++i) {

        if (std::isnan(probabilities[i]) || std::isinf(probabilities[i])) {
            return 0;
        }

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

    if (!passes_threshold) {
        return 0;
    }

    std::vector<double> adjusted_probs = probabilities;

    // Subtract thresholds from probabilities that pass
    for (size_t i = 0; i < adjusted_probs.size(); ++i) {
        if (adjusted_probs[i] > actual_thresholds[i]) {
            adjusted_probs[i] -= actual_thresholds[i];
        } else {
            adjusted_probs[i] = -1.0;
        }
    }

    auto max_it = std::max_element(adjusted_probs.begin(), adjusted_probs.end());

    if (max_it == adjusted_probs.end() || *max_it < 0.0) {
        return 0;
    }

    size_t max_index = std::distance(adjusted_probs.begin(), max_it);

    switch (max_index) {
        case 0: return 211;   // pip
        case 1: return 321;   // k
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
std::vector<std::vector<Double_t>> softmax_outputs_tof;

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
		softmax_outputs_tof.push_back(softmax_output_tof);			
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
void sv_reconstruction(Int_t N, std::string_view tag, Int_t SEED) // N - max event number to analyse
{
    SpdMCDataIterator* IT = new SpdMCDataIterator();

	std::string source{};
    if (tag == "signal") {
        source = "/WORKDIR/sig_data/reco_full_" + std::to_string(SEED) + ".root";
        IT->AddSourceFile(source.c_str());
    }
    else if (tag == "background") {
        source = "/WORKDIR/bg_data/reco_full_" + std::to_string(SEED) + ".root";
        IT->AddSourceFile(source.c_str());
    }

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
	
	Double_t Px_p{};
	Double_t Px_k{};
	Double_t Px_pip{};

	Double_t Py_p{};
	Double_t Py_k{};
	Double_t Py_pip{};

	Double_t Pz_p{};
	Double_t Pz_k{};
	Double_t Pz_pip{};

	Int_t true_decay{};

    Int_t id{};

	//======================================================================================================================
	// Output File
	TFile* file;
        if ( tag == "signal" )
                file = new TFile("analysed_signal.root", "RECREATE");
        else
                file = new TFile("analysed_background.root", "RECREATE");
	//======================================================================================================================
	// TTree
	TTree *tree = new TTree("tree", "tree");
	
	tree -> Branch("n_event", &n_event, "n_event/I");

    tree -> Branch("id", &id, "id/I");

	tree -> Branch("mass_Lc", &mass_Lc, "#Lambda_c^{+} mass/D");
			
	tree -> Branch("Px_p", &Px_p, "P momentum/D");
	tree -> Branch("Px_k", &Px_k, "P momentum/D");
	tree -> Branch("Px_pip", &Px_pip, "P momentum/D");

	tree -> Branch("Py_p", &Py_p, "P momentum/D");
	tree -> Branch("Py_k", &Py_k, "P momentum/D");
	tree -> Branch("Py_pip", &Py_pip, "P momentum/D");

	tree -> Branch("Pz_p", &Pz_p, "P momentum/D");
	tree -> Branch("Pz_k", &Pz_k, "P momentum/D");
	tree -> Branch("Pz_pip", &Pz_pip, "P momentum/D");

	tree -> Branch("true_decay", &true_decay, "#Lambda_c^{+} true decay/I");

	//======================================================================================================================
    while ( IT -> NextEvent() && n_event < events_max ) 
	{
		++n_event;		

        id = SEED * 10000 + n_event;

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
		PV_pos_prefit = prim_vtx_fit -> GetVertex();    // PV RC position first approximation.  
		SelectTracksForPV( PV_pos_prefit );             // Tracks selection for PV reconstruction
		if ( !findPVCheck )     continue;               // Skip if no PV.
		
		// if ( KFParticles_PV.size() < 5 ) continue;		// Check if have more than 5 tracks for PV reconstruction.

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

		primVtx.Construct(vDaughters_temp,nDaughters_temp, 0, -1);	// 0 --- no parant particle, 1 --- no mass hypothesis
		//=========================================================================================================================
		// Re-reconstructing primary vertex with KFParticle (MAIN ANALYSIS) VAR-1
        //=========================================================================================================================

        Double_t optimal_cut{5.0};
 
        Bool_t continue_removal{true};
 
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
		
		const KFParticle pr_vtx = primVtx;            					// init primary vertex
		//=========================================================================================================================
		// Futher analysis
		//=========================================================================================================================
		// Selecting tracks that are used for secondary vertex reconstruction
		SelectTracks(pr_vtx);
		//=========================================================================================================================

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
				if ( pdgs[vc[i]] == 2212 )
				{
					i_p = i;
				}
				if ( pdgs[vc[i]] == -321 )
				{
					i_K = i;
				}
				if ( pdgs[vc[i]] == 211 )
				{
					i_pip = i;
				}
			} 
			//============================================================================================
			// Tracks selection
			CheckTracks( d_particles, pdg ); 

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
			
			Px_p = p.GetPx();
			Px_k = K.GetPx();
			Px_pip = pip.GetPx();

			Py_p = p.GetPy();
			Py_k = K.GetPy();
			Py_pip = pip.GetPy();

			Pz_p = p.GetPz();
			Pz_k = K.GetPz();
			Pz_pip = pip.GetPz();

			const TVector3 momentum_Lc{Lc.GetPx(), Lc.GetPy(), 0 };
			
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

			//==================================================================================
			// Check cuts + Cut Flow
			if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 && tag == "background" && !true_decay )
                        	Total_count += 1;
			if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 && tag == "signal" && true_decay )
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
					if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 && tag == "background" && !true_decay)
						CutFlow_total[i] += 1; 	// Increasing because an event reaches i + 1 cuts step(passes through i + 1 cuts)
					if ( 2.24763 < mass_Lc && mass_Lc < 2.32497 && tag == "signal" && true_decay )
                        			CutFlow_total[i] += 1;
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
void analyze_Lc_pKpi_27_dalitz( Bool_t signal = 1, Int_t SEED = 0, Int_t N_max = 500 ) 
{   
	std::string_view tag{ signal ? "signal" : "background" };
	sv_reconstruction(N_max, tag, SEED);
}
