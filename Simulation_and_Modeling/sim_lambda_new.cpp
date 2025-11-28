#if defined(_COMPILE_MACRO_) 

#include <TRint.h>
#include <TStopwatch.h>
#include "FairParRootFileIo.h"
#include "SpdRunSim.h"
#include "SpdRootFileSink.h"
#include "SpdMCEventHeader.h"
#include "SpdFields.hh"
#include "SpdGenerators.hh"
#include "SpdCommonGeoMapper.h"
#include "SpdCave.h"
#include "SpdPipe.h"
#include "SpdHybMagnet.h"
#include "SpdIts.h"
#include "SpdTsTB.h"
#include "SpdTsTEC.h"
#include "SpdEcalTB2.h"
#include "SpdEcalTEC2.h"
#include "SpdRsTB2.h"
#include "SpdRsTEC2.h"
#include "SpdTofB.h"
#include "SpdTofEC.h"
#include "SpdZdc.h"
#include "SpdItsGeoMapperX.h"
#include "SpdTsTBGeoMapper.h"
#include "SpdTsTBGeoBuilder.h"
#include "SpdTsTECGeoMapper.h"

#endif


class MyFilter : public SpdPythia8Generator::P8EventFilter 
{
	virtual bool CheckEvent(SpdPythia8* pythia)
	{			
		Int_t particle_pdg{};
		TClonesArray* particles = new TClonesArray("TParticle");
		pythia -> ImportParticles(particles);		
		if ( !particles )
            		return 0;  	
		for (Int_t i{}; i < particles -> GetEntriesFast(); ++i) 
		{
            		TParticle* particle = (TParticle*) particles -> At(i);
			particle_pdg = particle -> GetPdgCode();		
			if ( particle_pdg == 4122 ) 
                		return 1;
		}
        	return 0;
	}
};


void set_decayer()
{
    SpdDecayer* spd_decayer = SpdDecayer::Instance();
    if (!spd_decayer->GetDecayer()) 
    spd_decayer->SetDecayer("SpdPythia8Decayer",1);
    SpdPythia8Decayer* decayer_lambda_c = dynamic_cast<SpdPythia8Decayer*>( spd_decayer -> GetDecayer() );

    decayer_lambda_c -> Init();
        
    //cout << "================== Forced decay mode for 4122 ===================" << endl;
    // Var I (Lc -> pK0s -> p(pip+pim))
    // decayer_lambda_c -> SelectForcedDecay(4122, 22); // DON'T CHANGE FOR ANTIPARTICLE!!!
    // decayer_lambda_c -> SelectForcedDecay(311, 1);
    // decayer_lambda_c -> SelectForcedDecay(310, 0);
	// decayer_lambda_c -> PrintParticleDecayChannels(4122);

	//cout << "================== Forced decay mode for 4122 ===================" << endl;
	// Var II (Lc -> p+K-pi+)
	decayer_lambda_c -> SelectForcedDecay(4122, 35);
   	
	decayer_lambda_c -> ForceSelectedDecays();
}


void sim_lambda_new( Int_t nEvents = 50, Int_t seed = -1, Bool_t signal = 1 )
{
    // std::string seed_str = "Random:seed = " + std::to_string(seed);
    TString outFile = Form("run_%d.root", seed);
    TString parFile = Form("params_%d.root", seed);

    gRandom->SetSeed(seed);
    SpdRunSim* run = new SpdRunSim();

    run->SetName("TGeant4");
    run->SetMaterials("media.geo");

    SpdRootFileSink* OutputFile = new SpdRootFileSink(outFile);
    run->SetSink(OutputFile);

    TString decayer_config = "DecayConfigPythia8.C";
    run->SetPythiaDecayer(decayer_config);

    run->SetMCEventHeader(new SpdMCEventHeader);

    /* +++++++++ GEOMETRY (QSOLENOID) ++++++++ */
    SpdCommonGeoMapper::Instance()->DefineQslGeometrySet();

    /* ++++++++++++++++++ CAVE ++++++++++++++++++ */
    FairModule* cave = new SpdCave("CAVE");
    run->AddModule(cave);

    /* ++++++++++++++++++ PIPE ++++++++++++++++++ */
    SpdPipe* Pipe = new SpdPipe();
    Pipe->SetPipeMaterial("beryllium",0);  // Central beampipe segment
    run->AddModule(Pipe);

    /* +++++++++ MAGNET +++++++++ */
    SpdSolMagnet2 *Magnet = new SpdSolMagnet2();
    run->AddModule(Magnet);

    /* +++++++++++++++++ DETECTORS ++++++++++++++ */
    SpdTsTB*  ts_barrel   = new SpdTsTB();
    SpdTsTEC* ts_ecaps    = new SpdTsTEC();
    SpdTofB* tof_barrel   = new SpdTofB();
    SpdTofEC* tof_ecaps   = new SpdTofEC();
    SpdEcalB *ecal_barrel = new SpdEcalB();
    SpdEcalEC *ecal_ecaps = new SpdEcalEC();
    SpdRsTB2 *rs_barrel   = new SpdRsTB2();
    SpdRsTEC2 *rs_ecaps   = new SpdRsTEC2();
    SpdBbc *bbc = new SpdBbc();
    SpdZdc *zdc = new SpdZdc();
    //SpdAeg* aeg = new SpdAeg();
    SpdFarich *farich = new SpdFarich(); /* +++++++++ FARICH  ++++++++++ */
    farich->setopticalphysics(true); // false - Farich is as material , true - in Farich create Cherenkov photons

    run->AddModule(ts_barrel);
    run->AddModule(ts_ecaps);
    run->AddModule(tof_barrel);
    run->AddModule(tof_ecaps);
    run->AddModule(ecal_barrel);
    run->AddModule(ecal_ecaps);
    run->AddModule(rs_barrel);
    run->AddModule(rs_ecaps);
    run->AddModule(bbc);
    run->AddModule(zdc);
    //run->AddModule(aeg);
    run->AddModule(farich);

    /* ===== Vertex detector ===== */
    SpdDssd* vd = new SpdDssd();
    run->AddModule(vd);

    //--------OLD magnetic field in 4.1.3-------------------------
    //SpdConstField* MagField = new SpdConstField(); 
    //MagField->SetField(0., 0., 10.); //kG, 1.1 T : TDR

    //latest magnetic field from 4.1.6 sample
    SpdFieldMap1_8 *MagField = new SpdFieldMap1_8("full_map");
    MagField->InitData("field_full1_8.bin");
    SpdRegion* reg = 0;
    reg = MagField->CreateFieldRegion("box");
    reg->SetBoxRegion(-330, 330, -330, 330, -386, 386);
    run->SetField(MagField);
    MagField->Print();

    /* ++++++++++ DEFINE PRIMARY GENERATORS +++++++++++ */
    SpdPrimaryGenerator* primGen = new SpdPrimaryGenerator();

    //-----------PYTHIA 8 GENERATOR-------------------------
    SpdPythia8Generator* P8gen = new SpdPythia8Generator();
    //define beams-------------

    P8gen -> SetBeam(2212, 2212, 27.); // pdg(A), pdg(B), E_cms (CM energy of collision)

    if ( signal )
    {
    	P8gen -> SetParameters("HardQCD:gg2ccbar = on");
    	P8gen -> SetParameters("HardQCD:qqbar2ccbar = on");
        P8gen->SetParameters("PhaseSpace:pTHatMin = 1.");
    	P8gen -> event_filter = new MyFilter();
    }
    else 		
    	P8gen -> SetParameters("SoftQCD:all = on");	// Background


    set_decayer();	  	
	
    P8gen -> SetVerboseLevel(-1);

    primGen -> AddGenerator(P8gen);

    //-----------------------------------
    run->SetGenerator(primGen);

    //smearing
    primGen->SetBeam(0., 0., 0.1, 0.1);
    primGen->SmearGausVertexXY(kTRUE);
    primGen->SetTarget(0., 30.);
    primGen->SmearGausVertexZ(kTRUE);

    primGen->SetVerboseLevel(-10);
    primGen->SetVerbose(0);

    /* >>>>>>>>>>>>>>>>>>>>>>>>>>> INITALIZE RUN <<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    run->Init();

    /* +++++++++++++ CREATE RUN PARAMETERS  +++++++++++ */
    Bool_t MergePars = kFALSE;
    FairParRootFileIo* parOut = new FairParRootFileIo(MergePars);
    if (MergePars) parOut->open(parFile.Data());
    else parOut->open(parFile.Data(),"RECREATE");
    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    rtdb->setOutput(parOut);

    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RUN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    TStopwatch timer;
    timer.Start();
    run->Run(nEvents);
    timer.Stop();

    /* ++++++++++++  SAVE RUN PARAMETERS  +++++++++++++ */
    rtdb->saveOutput();

    /*--------------------------------------------------*/
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    //-------------------------------------------------------------  
    cout << endl << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Output file is             " << outFile << endl;
    cout << "Parameter file is          " << parFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl;
    cout << endl;

    /*--------------------------------------------------*/
    SpdCommonGeoMapper::Instance()->PrintGeometry();

    /*--------------------------------------------------*/
    gApplication->Terminate();
}