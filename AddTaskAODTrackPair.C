AliAnalysisTaskAODTrackPair *AddTaskAODTrackPair(
						 UInt_t offlineTriggerMask = AliVEvent::kINT7,
						 string spd_multi_correction_file_path =
						 "/Users/syano_mbp2021/analysis/run2/Glueball/SPDMultiCorrection.root",
						 string dimuon_ds_file_path =
						 "/Users/syano_mbp2021/analysis/run2/Glueball/DownScale_Run2_CTP.root",
						 string k0s_cuts_file_path=
						 "/Users/syano_mbp2021/analysis/run2/Glueball/K0sSelectionCuts.root",
						 double min_vtxz = -10, double max_vtxz = 10, int min_vtx_cont = 1,
						 double min_pair_rap = -0.5,
						 double max_pair_rap = 0.5, 
						 double min_pair_pt = 0.50,
						 double max_pair_pt = 999., 
						 double alpha = 0.2,
						 double min_pangle = 0.998,
						 double max_pangle = 0.998,
						 double min_v0Dca = 0.1,
						 double max_v0Dca = 0.1,
						 double min_trackDca = 1.0,
						 double max_trackDca = 1.0,
						 double min_dlength = 5.0,
						 double max_dlength = 100.,
						 double min_lifetime = 20.,
						 double max_lifetime = 20.,
						 string period = "LHC18c",
						 string multi_method = "SPDTracklets", bool onPURej = true,
						 bool onLBcut = true, bool onMuEtaCut = true, bool onMuThetaAbsCut = true,
						 bool onMuMatchAptCut = true, bool onMuMatchLptCut = true,
						 bool onMuMatchHptCut = false, bool onMuChi2Cut = true,
						 bool onMuPdcaCut = true, bool isMC = false, bool isSelectEvt = true,
						 int paircuttype = 1, double min_pairtrackptcut = 0.5,
						 bool onMixingAnalysis = false, bool isMidTrackAnalysis = false,
						 bool isPrimTrackAnalysis = false, bool isV0TrackAnalysis = true,
						 double min_track_pt = 0.0, double max_track_pt = 999.,
						 double min_track_eta = -0.8, double max_track_eta = 0.8,
						 double min_track_p = 0.3, double max_track_p = 2.5,
						 double min_pion_sigma_tpc = -5, double max_pion_sigma_tpc = 5,
						 double min_pion_sigma_tof = -5, double max_pion_sigma_tof = 5,
						 double min_kaon_sigma_tpc = -2, double max_kaon_sigma_tpc = 2,
						 double min_kaon_sigma_tof = -2, double max_kaon_sigma_tof = 2,
						 double min_proton_sigma_tpc = -2, double max_proton_sigma_tpc = 2,
						 double min_proton_sigma_tof = -2, double max_proton_sigma_tof = 2,
						 double findable = 0.8, string dcaxy = "0.0105+0.035/pow(x,1.1)",
						 double dcaz = 2.0, double chi2tpc = 4., double chi2its = 36.,
						 int nclusttpc = 70, int nclustits = 1, double pair_opangle = 0.98,
						 int pid1 = 211, int pid2 = 211,
						 int trackdepth = 100,
						 int poolsize = 100,
						 double readypoolfraction = 0.1,
						 bool onpoolVtx = true,
						 bool onpoolCent = true,
						 bool onpoolPsi = true,
						 bool isManual = true,
						 bool decay_length_cut=true,
						 bool pointangle_cut=true,
						 bool chi2_cut=true,
						 bool pairDCA_cut=true,
						 bool track_distance_cut=true,
						 bool dcaxy_cut=true,
						 bool lifetime_cut=true,
						 bool armenteros_cut=true) {
    
  TFile *input3 = TFile::Open(k0s_cuts_file_path.c_str());
  
  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setMC(isMC);
  utils->setEvtSelection(isSelectEvt);
  utils->setVertexCut(min_vtxz, max_vtxz, min_vtx_cont);
  utils->setPairRapidityCut(min_pair_rap, max_pair_rap);
  utils->setV0SelectCuts(alpha,
			 min_pangle, 
			 max_pangle, 
			 min_v0Dca,
			 max_v0Dca,
			 min_trackDca,
			 max_trackDca,
			 min_dlength,
                         max_dlength,
			 min_lifetime,
			 max_lifetime);
  utils->setPileupRejectionCut(onPURej);
  utils->setPairKinematicCut(paircuttype, min_pairtrackptcut);
  utils->setTrackKinematicCut(min_track_pt, max_track_pt, min_track_eta, max_track_eta, min_track_p, max_track_p);
  utils->setPionSelectSigmaTPC(min_pion_sigma_tpc, max_pion_sigma_tpc);
  utils->setKaonSelectSigmaTPC(min_kaon_sigma_tpc, max_kaon_sigma_tpc);
  utils->setProtonSelectSigmaTPC(min_proton_sigma_tpc, max_proton_sigma_tpc);
  utils->setPionSelectSigmaTOF(min_pion_sigma_tof, max_pion_sigma_tof);
  utils->setKaonSelectSigmaTOF(min_kaon_sigma_tof, max_kaon_sigma_tof);
  utils->setProtonSelectSigmaTOF(min_proton_sigma_tof, max_proton_sigma_tof);
  utils->setTrackQualities(findable, dcaxy.c_str(), dcaz, chi2tpc, chi2its, nclusttpc, nclustits);
  utils->setPairTargetPIDs(pid1, pid2);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAODMuonEventSelection",
            "No analysis manager to connect to");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAODMuonEventSelection",
            "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskAODTrackPair *task = new AliAnalysisTaskAODTrackPair("dimuon");
  if (isSelectEvt) {
    task->SelectCollisionCandidates(offlineTriggerMask);
  }

  task->setUtils(utils);
  task->setMC(isMC);
  task->setCentralityMethod(multi_method);
  task->setEvtMixingTrackDepth(trackdepth);
  task->setEvtMixingPoolSize(poolsize);
  task->setEvtMixingReadyFraction(readypoolfraction);
  task->setEvtMixingPoolVtxZ(onpoolVtx);
  task->setEvtMixingPoolCent(onpoolCent);
  task->setEvtMixingPoolPsi(onpoolPsi);
  task->setMixingAnalysis(onMixingAnalysis);
  task->setVtxZRange(min_vtxz,max_vtxz);
  task->setManualV0Analysis(isManual);  
  task->setK0sCuts(input3,decay_length_cut,pointangle_cut,chi2_cut,pairDCA_cut,track_distance_cut,dcaxy_cut,lifetime_cut,armenteros_cut);  
  task->setK0sRapRange(min_pair_rap,max_pair_rap);
  task->setK0sPtRange(min_pair_pt,max_pair_pt);
  task->setK0sDaughterTrackPRange(min_track_p,max_track_p);
  task->setK0sDaughterTrackPtRange(min_track_pt,max_track_pt);
  task->setK0sDaughterTrackEtaRange(min_track_eta,max_track_eta);
  task->setK0sDaughterTrackChi2TPC(chi2tpc);
  task->setK0sDaughterTrackChi2ITS(chi2its);
  task->setK0sDaughterTrackNClustTPC(nclusttpc);
  task->setK0sDaughterTrackNClustSPD(nclustits);
  task->setK0sDaughterTrackFindableTPC(findable);
  mgr->AddTask(task);

  cout << "min_vtxz=" << min_vtxz << endl;
  cout << "max_vtxz=" << max_vtxz << endl;
  cout << "min_vtx_cont=" << min_vtx_cont << endl;
  cout << "min_pair_rap=" << min_pair_rap << endl;
  cout << "max_pair_rap=" << max_pair_rap << endl;
  cout << "alpha=" << alpha << endl;
  cout << "min_pangle=" << min_pangle << endl;
  cout << "max_pangle=" << max_pangle << endl;
  cout << "min_v0Dca=" << min_v0Dca << endl;
  cout << "max_v0Dca=" << max_v0Dca << endl;
  cout << "min_trackDca=" << min_trackDca << endl;
  cout << "max_trackDca=" << max_trackDca << endl;
  cout << "min_dlength=" << min_dlength << endl;
  cout << "max_dlength=" << max_dlength << endl;
  cout << "min_lifetime=" << min_lifetime << endl;
  cout << "max_lifetime=" << max_lifetime << endl;
  cout << "period=" << period << endl;
  cout << "multi_method=" << multi_method << endl;
  cout << "onPURej=" << onPURej << endl;
  cout << "onLBcut=" << onLBcut << endl;
  cout << "onMuEtaCut=" << onMuEtaCut << endl;
  cout << "onMuThetaAbsCut=" << onMuThetaAbsCut << endl;
  cout << "onMuMatchAptCut=" << onMuMatchAptCut << endl;
  cout << "onMuMatchLptCut=" << onMuMatchLptCut << endl;
  cout << "onMuMatchHptCut=" << onMuMatchHptCut << endl;
  cout << "onMuChi2Cut=" << onMuChi2Cut << endl;
  cout << "onMuPdcaCut=" << onMuPdcaCut << endl;
  cout << "isMC=" << isMC << endl;
  cout << "isSelectEvt=" << isSelectEvt << endl;
  cout << "paircuttype=" << paircuttype << endl;
  cout << "min_pairtrackptcut=" << min_pairtrackptcut << endl;
  cout << "onMixingAnalysis=" << onMixingAnalysis << endl;
  cout << "isMidTrackAnalysis=" << isMidTrackAnalysis << endl;
  cout << "isPrimTrackAnalysis=" << isPrimTrackAnalysis << endl;
  cout << "isV0TrackAnalysis=" << isV0TrackAnalysis << endl;
  cout << "min_track_pt=" << min_track_pt << endl;
  cout << "max_track_pt=" << max_track_pt << endl;
  cout << "min_track_eta=" << min_track_eta << endl;
  cout << "max_track_eta=" << max_track_eta << endl;
  cout << "min_track_p=" << min_track_p << endl;
  cout << "max_track_p=" << max_track_p << endl;
  cout << "min_pion_sigma_tpc=" << min_pion_sigma_tpc << endl;
  cout << "max_pion_sigma_tpc=" << max_pion_sigma_tpc << endl;
  cout << "min_pion_sigma_tpc=" << min_pion_sigma_tof << endl;
  cout << "max_pion_sigma_tpc=" << max_pion_sigma_tof << endl;
  cout << "min_kaon_sigma_tpc=" << min_kaon_sigma_tpc << endl;
  cout << "max_kaon_sigma_tpc=" << max_kaon_sigma_tpc << endl;
  cout << "min_kaon_sigma_tpc=" << min_kaon_sigma_tof << endl;
  cout << "max_kaon_sigma_tpc=" << max_kaon_sigma_tof << endl;
  cout << "min_proton_sigma_tpc=" << min_proton_sigma_tpc << endl;
  cout << "max_proton_sigma_tpc=" << max_proton_sigma_tpc << endl;
  cout << "min_proton_sigma_tpc=" << min_proton_sigma_tof << endl;
  cout << "max_proton_sigma_tpc=" << max_proton_sigma_tof << endl;
  cout << "findable=" << findable << endl;
  cout << "dcaxy=" << dcaxy.c_str() << endl;
  cout << "dcaz=" << dcaz << endl;
  cout << "chi2tpc=" << chi2tpc << endl;
  cout << "chi2its=" << chi2its << endl;
  cout << "nclusttpc=" << nclusttpc << endl;
  cout << "nclustits=" << nclustits << endl;
  cout << "pair_opangle=" << pair_opangle << endl;
  cout << "pid1=" << pid1 << endl;
  cout << "pid2=" << pid2 << endl;
  cout << "trackdepth=" << trackdepth << endl;
  cout << "poolsize=" << poolsize << endl;
  cout << "readypoolfraction=" << readypoolfraction << endl;
  cout << "onpoolVtx=" << onpoolVtx << endl;
  cout << "onpoolCent=" << onpoolCent << endl;
  cout << "onpoolPsi=" << onpoolPsi << endl;
  cout << "isManual=" << isManual << endl;
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("dimuon", TList::Class(),
                           AliAnalysisManager::kOutputContainer, "Dimuon.root");

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}
