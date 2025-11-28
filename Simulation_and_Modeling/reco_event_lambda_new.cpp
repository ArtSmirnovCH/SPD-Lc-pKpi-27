void reco_event_lambda_new(Int_t seed = -1)
{
  	Int_t nEvents{}; 
  
 	TString inFile = Form("run_%d.root", seed);			// Input file (MC data)
    TString parFile = Form("params_%d.root", seed);			// Input file with parameters
    TString outFile = Form("reco_full_%d.root", seed); 			// Output file
  
	gRandom->SetSeed(seed);

    SpdRunAna* Run = new SpdRunAna();
    Run->SetContainerStatic(true);

    FairFileSource* SourceFile = new FairFileSource(inFile);
    Run->SetSource(SourceFile);

    SpdRootFileSink* OutputFile = new SpdRootFileSink(outFile);
    Run->SetSink(OutputFile);

    FairRuntimeDb* rtdb = Run->GetRuntimeDb();
    FairParRootFileIo* parInput = new FairParRootFileIo();
    parInput->open(parFile.Data());
    rtdb->setFirstInput(parInput);

    TStopwatch timer;
    timer.Start();

    //-------------------------------------------------------//
    SpdMCEventMaker* event_maker = new SpdMCEventMaker();
    event_maker->SetVerboseLevel(1);
    Run->AddTask(event_maker);

    //  //--------------------------------------------------------//
    //  SpdItsMCHitProducer* its_hits_producer = new SpdItsMCHitProducer();
    //  its_hits_producer->SetVerboseLevel(1);
    //  Run->AddTask(its_hits_producer);

    SpdDssdMCHitProducer* vd_hits_producer = new SpdDssdMCHitProducer();
    vd_hits_producer->SetVerboseLevel(1);
    Run->AddTask(vd_hits_producer);

    //--------------------------------------------------------//
    SpdTsMCHitProducer* ts_hits_producer = new SpdTsMCHitProducer();
    ts_hits_producer->SetVerboseLevel(1);
    Run->AddTask(ts_hits_producer);

    //--------------------------------------------------------//
    SpdTofMCHitProducer* tof_hits_producer = new SpdTofMCHitProducer();
    tof_hits_producer->SetVerboseLevel(1);
    Run->AddTask(tof_hits_producer);

    //--------------------------------------------------------//
    SpdEcalMCHitProducer* ecal_hits_producer = new SpdEcalMCHitProducer();
    ecal_hits_producer->SetVerboseLevel(1);
    Run->AddTask(ecal_hits_producer);

    //--------------------------------------------------------//
    SpdRsMCHitProducer* rs_hits_producer = new SpdRsMCHitProducer();
    rs_hits_producer->SetVerboseLevel(1);
    Run->AddTask(rs_hits_producer);

    //--------------------------------------------------------//
    SpdRsMCRHitProducer* rs_rhits_producer = new SpdRsMCRHitProducer();  // realistic hits
    rs_rhits_producer->SetVerboseLevel(1);
    Run->AddTask(rs_rhits_producer);

    //--------------------------------------------------------//
    SpdBbcMCHitProducer* bbc_hits_producer = new SpdBbcMCHitProducer();
    bbc_hits_producer->SetVerboseLevel(1);
    Run->AddTask(bbc_hits_producer);

    //--------------------------------------------------------//
    SpdZdcMCHitProducer* zdc_hits_producer = new SpdZdcMCHitProducer();
    zdc_hits_producer->SetVerboseLevel(1);
    Run->AddTask(zdc_hits_producer);

    //  //--------------------------------------------------------//
    //  SpdAegMCHitProducer* aeg_hits_producer = new SpdAegMCHitProducer();
    //  aeg_hits_producer->SetVerboseLevel(1);
    //  Run->AddTask(aeg_hits_producer);

    //--------------------------------------------------------//
    SpdFarichMCHitProducer *farich_hits_producer = new SpdFarichMCHitProducer();
    farich_hits_producer->SetVerboseLevel(1);
    Run->AddTask(farich_hits_producer);

    //--------------------------------------------------------//
    SpdMCTracksFinder* track_finder = new SpdMCTracksFinder();
    track_finder->SetVerboseLevel(1);

    SpdTrackFitterGF* track_fitter = track_finder->Fitter();
    track_fitter->SetVerboseLevel(0);

    track_finder->CheckMinItsHits(true,0);
    track_finder->CheckMinTsHits(true,6);
    track_finder->CheckMinHits(true,7);
    track_finder->CheckMaxPartGeneration(true, 3);
    track_finder->CheckMinPartPt(true,0.1);
    track_finder->CheckMinPartMomentum(true,0.15);

    Run->AddTask(track_finder);

    //--------------------------------------------------------//  
    SpdMCVerticesFitter* mcvtxs_fitter = new SpdMCVerticesFitter();
    mcvtxs_fitter->SetFitSecondaries(false);
    mcvtxs_fitter->SetVerboseLevel(1);
    Run->AddTask(mcvtxs_fitter);

    //--------------------------------------------------------//  
    SpdRCVerticesFinder* rcvtxs_finder = new SpdRCVerticesFinder();
    //rcvtxs_finder->SetMinItsHits(0);//for no SVD case only  
    rcvtxs_finder->SetMinItsHits(2);
    rcvtxs_finder->SetFitSecondaries(false);
    rcvtxs_finder->SetVerboseLevel(1);
    Run->AddTask(rcvtxs_finder);

    //--------------------------------------------------------//
    SpdEcalRCMaker* ecal_rc = new SpdEcalRCMaker();
    ecal_rc->SetVerboseLevel(1);
    Run->AddTask(ecal_rc);

    //--------------------------------------------------------//
    SpdEcalClusterMCInfoMaker* ecal_mc_info = new SpdEcalClusterMCInfoMaker();
    ecal_mc_info->SetVerboseLevel(1);
    Run->AddTask(ecal_mc_info);

    //--------------------------------------------------------//
    SpdRsMCClusterMaker* rs_cl = new SpdRsMCClusterMaker();
    rs_cl->SetVerboseLevel(1);
    Run->AddTask(rs_cl);

    //--------------------------------------------------------//
    SpdMCTsParticleProducer* mcts_part = new SpdMCTsParticleProducer();
    mcts_part->SetVerboseLevel(1);
    Run->AddTask(mcts_part);

    //--------------------------------------------------------//
    SpdMCTofParticleProducer* mctof_part = new SpdMCTofParticleProducer();
    mctof_part->SetVerboseLevel(1);
    Run->AddTask(mctof_part);

    //  //--------------------------------------------------------//  
    //  SpdMCAegParticleProducer* mcaeg_part = new SpdMCAegParticleProducer();
    //  mcaeg_part->SetVerboseLevel(1);
    //  Run->AddTask(mcaeg_part);

    ///-------------------FARICH-------------------------------------//
    // [MC-PRODUCER FOR FARICH-PARTICLES]
    // Output: FARICH-particles
    SpdMCFarichParticleProducer *mcfarich_part = new SpdMCFarichParticleProducer();
    mcfarich_part->SetVerboseLevel(1);
    Run->AddTask(mcfarich_part);

    //========================================================================// 
    cout << "Start ... " << endl;
    Run->Initialize();

    //-------------------------------------//     
    cout << "Run ... " << endl;
    Run->Run(0,nEvents);

    //-------------------------------------//  
    OutputFile->Close();
    timer.Stop();

    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();

    cout << endl << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
    cout << endl;
}
