AliAnalysisTaskAODTrackPairMC* AddTaskAODTrackPairMC(uint offlineTriggerMask = AliVEvent::kINT7,
						     float min_vtxz =-10,
						     float max_vtxz = 10,
						     int min_vtx_cont = 1,
						     float min_pair_rap = -4.0,
						     float max_pair_rap = -2.5,
						     string period = "LHC18c",
						     string multi_method="SPDTracklets",
						     bool onPURej = true,
						     bool onLBcut = true,
						     bool onMuEtaCut = true,
						     bool onMuThetaAbsCut = true,
						     bool onMuMatchAptCut = true,
						     bool onMuMatchLptCut = true,
						     bool onMuMatchHptCut = false,
						     bool onMuChi2Cut = true,
						     bool onMuPdcaCut = true,
						     bool isMC=false,
						     bool isSelectEvt=true,
						     int paircuttype=1,
						     double min_pairtrackptcut=0.5,
						     bool onMixingAnalysis=false,
						     bool isMidMuonAnalysis=false
						     )

{
    
  AliMuonTrackCuts* fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetIsMC(isMC);
  fMuonTrackCuts->SetAllowDefaultParams(true);
  int selectionMask = 0;
  if(onMuEtaCut) selectionMask |= AliMuonTrackCuts::kMuEta;
  if(onMuThetaAbsCut) selectionMask |=AliMuonTrackCuts::kMuThetaAbs;
  if(onMuMatchAptCut) selectionMask |=AliMuonTrackCuts::kMuMatchApt;
  if(onMuMatchLptCut) selectionMask |=AliMuonTrackCuts::kMuMatchLpt;
  if(onMuMatchHptCut) selectionMask |=AliMuonTrackCuts::kMuMatchHpt;
  if(onMuPdcaCut) selectionMask |=AliMuonTrackCuts::kMuPdca;
  if(onMuChi2Cut) selectionMask |=AliMuonTrackCuts::kMuTrackChiSquare;    
  fMuonTrackCuts->SetFilterMask(selectionMask);
  
  TFile* input = TFile::Open("./DownScale_Run2_CTP.root");
  TFile* input2 = TFile::Open("./SPDMultiCorrection.root");

  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setMC(isMC);
  utils->setEvtSelection(isSelectEvt);
  utils->setDownScalingHist(input);
  utils->setSPDTrkCorrHist(input2,period);
  utils->setVertexCut(min_vtxz,max_vtxz,min_vtx_cont);
  utils->setPairRapidityCut(min_pair_rap,max_pair_rap);
  utils->setV0SelectCuts(alpha,pangle,v0Dca,trackDca,min_dlength,max_dlength);
  utils->setPileupRejectionCut(onPURej);
  utils->setLocalBoardCut(onLBcut);
  utils->setMultiEstimateMethod(multi_method);
  utils->setMuonTrackCut(fMuonTrackCuts);
  utils->setPairKinematicCut(paircuttype,min_pairtrackptcut);
  utils->setMidMuonAna(isMidMuonAnalysis);  
  utils->setPeriod(period);
  if (isMidMuonAnalysis) {
    utils->setMidTrackKinematicRange(0.05,999,min_track_eta,max_track_eta);
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAODMuonEventSelection", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAODMuonEventSelection", "This task requires an input event handler");
    return NULL;
  }
    
  AliAnalysisTaskAODTrackPairMC* task = new AliAnalysisTaskAODTrackPairMC("dimuon");
  if(isSelectEvt){
    task->SelectCollisionCandidates(offlineTriggerMask);
  }
  task->setUtils(utils);
  task->setEvtMixingTrackDepth(100);
  task->setEvtMixingPoolSize(100);
  task->setEvtMixingReadyFraction(0.1);
  task->setEvtMixingPoolVtxZ(true);
  task->setEvtMixingPoolCent(true);
  task->setEvtMixingPoolPsi(true);
  task->setMixingAnalysis(onMixingAnalysis);
  task->setMC(isMC);
  task->setMidMuonAna(isMidMuonAnalysis);
  mgr->AddTask(task);
  
  cout<<"min_vtxz="<< min_vtxz <<endl;
  cout<<"max_vtxz="<< max_vtxz <<endl;
  cout<<"min_vtx_cont="<<  min_vtx_cont <<endl;
  cout<<"min_pair_rap="<< min_pair_rap <<endl;
  cout<<"max_pair_rap="<< max_pair_rap <<endl;
  cout<<"min_track_eta="<< min_track_eta <<endl;
  cout<<"max_track_eta="<< max_track_eta <<endl;
  cout<<"alpha="<< alpha <<endl;
  cout<<"pangle="<< pangle <<endl;
  cout<<"v0Dca="<< v0Dca <<endl;
  cout<<"trackDca="<< trackDca <<endl;
  cout<<"min_dlength="<< min_dlength <<endl;
  cout<<"max_dlength="<< max_dlength <<endl;
  cout<<"period="<< period <<endl;
  cout<<"multi_method="<< multi_method <<endl;
  cout<<"onPURej="<< onPURej <<endl;
  cout<<"onLBcut="<< onLBcut <<endl;
  cout<<"onMuEtaCut="<< onMuEtaCut <<endl;
  cout<<"onMuThetaAbsCut="<< onMuThetaAbsCut <<endl;
  cout<<"onMuMatchAptCut="<< onMuMatchAptCut <<endl;
  cout<<"onMuMatchLptCut="<< onMuMatchLptCut <<endl;
  cout<<"onMuMatchHptCut="<< onMuMatchHptCut <<endl;
  cout<<"onMuChi2Cut="<< onMuChi2Cut <<endl;
  cout<<"onMuPdcaCut="<< onMuPdcaCut <<endl;
  cout<<"isMC="<< isMC <<endl;
  cout<<"isSelectEvt="<< isSelectEvt <<endl;
  cout<<"paircuttype="<< paircuttype <<endl;
  cout<<"min_pairtrackptcut="<< min_pairtrackptcut <<endl;
  cout<<"onMixingAnalysis="<< onMixingAnalysis <<endl;
  cout<<"isMidMuonAnalysis="<< isMidMuonAnalysis <<endl;
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("dimuon",TList::Class(),AliAnalysisManager::kOutputContainer,"Dimuon.root");
  
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput1);

  return task;

}
