#include <SpdMCParticle.h>
#include <SpdTrackMC.h>
#include <SpdMCTrack.h>

#include <SpdTrackFitPar.h>
#include <SpdVertexRC.h>

#include <SpdTrackFitPar.h>

#include <SpdAegMCHit.h>

#include "TFile.h"
#include "TTree.h"

#include "TVector3.h"

void PID_ana() 
{    
    SpdMCDataIterator* IT = 0;
    const SpdSetParSet* Pars_ = 0;

    const TClonesArray*   mcparticles     = 0;     
    const TClonesArray*   mctracks        = 0;     
    const TClonesArray*   rcecalpartilces = 0;    
    const TClonesArray*   rcvertices = 0; 
    const TClonesArray*   particles_aeg = 0;
    const TClonesArray*   mcvertices = 0;   
    const TClonesArray *particles_tof = 0; 
    const TClonesArray *particles_ts = 0;

    IT = new SpdMCDataIterator();
 
    IT -> AddSourceFile("reco_full_34.root");    
    
    IT -> ActivateBranch("all");
         
    IT -> Init();
          
    Pars_ = IT -> GetParameters(); 
    
    mcparticles = IT -> GetParticles();    
    mctracks = IT -> GetTracks();         
    rcecalpartilces = IT -> GetEcalParticlesRC();
    rcvertices = IT -> GetVerticesRC();
    mcvertices = IT -> GetVertices();
    particles_tof = IT -> GetTofParticles();
    particles_ts = IT -> GetTsParticles();

    Int_t events_max{5000}; 
    Int_t n_event{};

    TFile* file = new TFile("PID_ana.root", "RECREATE");

    TTree *tree = new TTree("tree", "tree");

    // vars
    Int_t mc_pid{};
    
    // likelihoods (gaussian): 0 - pion,  1 - kaon, 2 - proton
    Double_t ll_tof_pi{};
    Double_t ll_tof_k{};
    Double_t ll_tof_p{};

    Double_t ll_ts_pi{};
    Double_t ll_ts_k{};
    Double_t ll_ts_p{};

    Double_t px{};
    Double_t py{};
    Double_t pz{};

    Double_t p{};
    Double_t pt{};
    Double_t eta{};
    Double_t phi{};
    Double_t theta{};

    Double_t chi2overndf{};
    Int_t nhitsits{};
    Int_t nhitstsb{};
    Int_t nhitstsec{};
    Int_t isfitted{};
    Int_t isgood{};
    Int_t convergency{};

    tree -> Branch("mc_pid", &mc_pid, "mc_pid/I");

    tree->Branch("ll_tof_pi", &ll_tof_pi, "ll_tof_pi/D");
    tree->Branch("ll_tof_k", &ll_tof_k, "ll_tof_k/D");
    tree->Branch("ll_tof_p", &ll_tof_p, "ll_tof_p/D");

    tree->Branch("ll_ts_pi", &ll_ts_pi, "ll_ts_pi/D");
    tree->Branch("ll_ts_k", &ll_ts_k, "ll_ts_k/D");
    tree->Branch("ll_ts_p", &ll_ts_p, "ll_ts_p/D");

    tree->Branch("px", &px, "px/D");
    tree->Branch("py", &py, "py/D");
    tree->Branch("pz", &pz, "pz/D");

    tree->Branch("p", &p, "p/D");
    tree->Branch("pt", &pt, "pt/D");
    tree->Branch("eta", &eta, "eta/D");
    tree->Branch("phi", &phi, "phi/D");
    tree->Branch("theta", &theta, "theta/D");

    tree->Branch("chi2overndf", &chi2overndf, "chi2overndf/D");
    tree->Branch("nhitsits", &nhitsits, "nhitsits/I");
    tree->Branch("nhitstsb", &nhitstsb, "nhitstsb/I");
    tree->Branch("nhitstsec", &nhitstsec, "nhitstsec/I");
    tree->Branch("isfitted", &isfitted, "isfitted/I");
    tree->Branch("isgood", &isgood, "isgood/I");
    tree->Branch("convergency", &convergency, "convergency/I");


    while ( (IT->NextEvent()) && (n_event < events_max) )
    {   
        
        std::cout << "=======================================================================" << '\n';
        std::cout << "NEW EVENT:  " << n_event << '\n';

        for (Int_t i{}; i < mctracks -> GetEntriesFast(); ++i){
        
            SpdTrackMC* track{ dynamic_cast<SpdTrackMC*>( mctracks->At( i ) ) };
            
            SpdMCParticle* particle{ dynamic_cast<SpdMCParticle*>( mcparticles -> At( track -> GetParticleId() ) ) };

            TVector3 momentum = particle -> GetStartMom();

            px = momentum.X();     
            py = momentum.Y();     
            pz = momentum.Z();     
            pt = momentum.Pt();
            p = momentum.Mag();   
            
            theta = particle -> GetStartPtheta();
            phi = particle -> GetStartPphi();
            eta = -log(tan(theta / 2.0));

            nhitsits = track -> GetNHitsIts();
            nhitstsb = track -> GetNHitsTsB();
            nhitstsec = track -> GetNHitsTsEC();
            isfitted = track -> GetIsFitted();

            SpdTrackFitPar *track_pars{ track -> GetFitPars() };
            if ( !track_pars ) continue;

            chi2overndf = track_pars -> GetChi2overNDF();
            isgood = track_pars -> GetIsGood();            
            convergency = track_pars -> GetConvergency();            

            // PID

            mc_pid = track -> GetParticlePdg();
            
            // TOF
            
            vector <Double_t> Likelihoods_tof;              
            if ( particle -> GetTofParticleId() != -1 ){
                SpdTofParticle *tofparticle = dynamic_cast<SpdTofParticle*>( particles_tof -> At( particle -> GetTofParticleId() ) );
                if ( !tofparticle ) continue;

                Likelihoods_tof = tofparticle -> GetLikelihoods();  
            }
            else continue;
            
            // dE / dx
            
            vector <Double_t> Likelihoods_ts;       
            if ( particle -> GetTsParticleId() != -1 ){   
                SpdTsParticle *ftsparticle = dynamic_cast<SpdTsParticle*>( particles_ts -> At( particle -> GetTsParticleId() ) );
                if ( !ftsparticle ) continue; 
                
                Likelihoods_ts = ftsparticle -> GetLikelihoods();
            }
            else continue;

            ll_tof_pi = Likelihoods_tof[0];
            ll_tof_k = Likelihoods_tof[1];
            ll_tof_p = Likelihoods_tof[2];

            ll_ts_pi = Likelihoods_ts[0];
            ll_ts_k = Likelihoods_ts[1];
            ll_ts_p = Likelihoods_ts[2];            
            
            // Tree Fill 
            tree -> Fill();
        }
        ++n_event;
    
    }

    file -> Write();
    file -> Close();

    std::cout << "NEvents: " << n_event << std::endl;

}
