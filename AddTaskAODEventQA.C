AliAnalysisTaskAODEventStudy* AddTaskAODEventQA(UInt_t offlineTriggerMask = AliVEvent::kAny,
						float min_vtxz =-10.,
						float max_vtxz = 10.,
						float min_pair_rap = -4.0,
						float max_pair_rap = -2.5,
						string multi_method="SPDTracklets",
						bool onPURej = true,
						bool onLBcut = true,
						bool onMuEtaCut = true,
						bool onMuThetaAbsCut = true,
						bool onMuMatchAptCut = true,
						bool onMuMatchLptCut = true,
						bool onMuMatchHptCut = true,
						bool onMuChi2Cut = true,
						bool onMuPdcaCut = true,
						bool isMC=false)
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

  AliAnalysisTaskAODTrackPairUtils *utils = new AliAnalysisTaskAODTrackPairUtils();
  utils->setDownScalingHist(input);
  utils->setVertexCut(min_vtxz,max_vtxz);
  utils->setPairRapidityCut(min_pair_rap,max_pair_rap);
  utils->setPileupRejectionCut(onPURej);
  utils->setLocalBoardCut(onLBcut);
  utils->setMultiEstimateMethod(multi_method);
  utils->setMuonTrackCut(fMuonTrackCuts);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAODMuonEventSelection", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAODMuonEventSelection", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskAODEventStudy* task = new AliAnalysisTaskAODEventStudy("runQA");
  task->SelectCollisionCandidates(offlineTriggerMask);
  task->setUtils(utils);
  task->setMC(isMC);
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("taskQA",TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput1);
  
  return task;

}
