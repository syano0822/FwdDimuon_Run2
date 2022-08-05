#include "AliAnalysisTaskAODTrackPairUtils.h"
#include "AliAnalysisTaskAODEventStudy.h"
#include "AliAnalysisTaskAODTrackPair.h"

AliAnalysisGrid* CreateAlienHandler(TString period, TString run_mode, Bool_t isJDL, TString type,bool onMixingAnalysis);

void runAnalysisMidMuonMC(TString runPeriod = "LHC16k",
			  TString run_mode  = "test",
			  Bool_t isJDL      = true,
			  TString type      = "data",
			  Bool_t  local     = false,
			  bool isMix        = false)
{
  // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
  gInterpreter->ProcessLine(".include ./");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include ./");
#endif
  
  bool isMC = false;  
  bool isEventSelection = true;

  if(type != "data"){
    isMC = true;
    if(type == "LHC20f10a" || type == "LHC20f10b" || type == "LHC20f10c" || type == "Pythia6Perugia2011CCSemiMuonic" ||
       type == "UncorrFlatMuons_GEANT4" ){
      isEventSelection = false;
    }
  }
  
  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

#if !defined (__CINT__) || defined (__CLING__)
  if(isEventSelection) {    
    if(run_mode=="full" || run_mode=="test" ) {
      AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d,%d)",isMC,true)));
    }
  }
#else
  if(isEventSelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);
  }
#endif

#if !defined (__CINT__) || defined (__CLING__)
  if(isEventSelection){
    if(run_mode=="full" || run_mode=="test" ) {
      AliMultSelectionTask *multSelTask=reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C(%d)",false)));
      multSelTask->SetUseDefaultCalib(true);
      if(isMC){
	multSelTask->SetUseDefaultMCCalib(true);
      } else{
	multSelTask->SetUseDefaultCalib(true); 
      }
    }
  }
#else
  if(isEventSelection){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *multSelTask = AddTaskMultSelection(false);
    multSelTask->SetUseDefaultCalib(true);
    if(isMC) {
      multSelTask->SetUseDefaultMCCalib(true);
    } else {
      multSelTask->SetUseDefaultCalib(true);
    }
  }
#endif

#if !defined (__CINT__) || defined (__CLING__) //ROOT6
  AliAnalysisTaskPIDResponse *pidtask=reinterpret_cast<AliAnalysisTaskPIDResponse*>
    (gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C()"));
#else //ROOT5
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidtask=AddTaskPIDResponse();
#endif
 
  TFile* input = TFile::Open("./DownScale_Run2_CTP.root");
  
  unsigned int offlineTriggerMask;
  if(isEventSelection){
    offlineTriggerMask = AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7;
  } else {
    offlineTriggerMask = AliVEvent::kAny;
  }  

  float min_vtxz =-10;
  float max_vtxz = 10;
  int min_vtx_cont = 1;
  float min_pair_rap = -0.8;
  float max_pair_rap =  0.8;
  string multi_method="SPDTracklets";
  bool onPURej = true;
  bool onLBcut = false;
  bool onMuEtaCut = true;
  bool onMuThetaAbsCut = false;
  bool onMuMatchAptCut = false;
  bool onMuMatchLptCut = false;
  bool onMuMatchHptCut = false;
  bool onMuChi2Cut = true;
  bool onMuPdcaCut = false;
  bool isSelectEvt = isEventSelection;
  int paircuttype = 0;
  float min_pairtrackptcut = 0.0;  
  bool onMixingAnalysis = isMix;
  bool isMidMuonAnalysis = true;

#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gInterpreter->LoadMacro("AliAnalysisTaskAODEventStudy.cxx++g");
  gInterpreter->LoadMacro("AliAnalysisTaskAODTrackPair.cxx++g");
  AliAnalysisTaskAODTrackPair *task
    = reinterpret_cast<AliAnalysisTaskAODTrackPair*>(gInterpreter->ExecuteMacro(Form("AddTaskAODTrackPair.C(%u, %f, %f, %d, %f, %f, \"%s\", %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d,%d,%f,%d,%d)",
										     offlineTriggerMask,
										     min_vtxz,
										     max_vtxz,
										     min_vtx_cont,
										     min_pair_rap,
										     max_pair_rap,
										     multi_method.c_str(),
										     onPURej,
										     onLBcut,
										     onMuEtaCut,
										     onMuThetaAbsCut,
										     onMuMatchAptCut,
										     onMuMatchLptCut,
										     onMuMatchHptCut,
										     onMuChi2Cut,
										     onMuPdcaCut,
										     isMC,
										     isSelectEvt,
										     paircuttype,
										     min_pairtrackptcut,
										     onMixingAnalysis,
										     isMidMuonAnalysis
										     )));
#else
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskAODEventStudy.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPair.cxx++g");
  gROOT->LoadMacro("AddTaskAODTrackPair.C");
  AliAnalysisTaskAODTrackPair* task = AddTaskAODTrackPair(offlineTriggerMask,
							  min_vtxz,
							  max_vtxz,
							  min_vtx_cont,
							  min_pair_rap,
							  max_pair_rap,
							  multi_method.c_str(),
							  onPURej,
							  onLBcut,
							  onMuEtaCut,
							  onMuThetaAbsCut,
							  onMuMatchAptCut,
							  onMuMatchLptCut,
							  onMuMatchHptCut,
							  onMuChi2Cut,
							  onMuPdcaCut,
							  isMC,
							  isSelectEvt,
							  paircuttype,
							  min_pairtrackptcut,
							  onMixingAnalysis,
							  isMidMuonAnalysis
							  );
#endif
  
  if(!mgr->InitAnalysis()) return;

  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);
 
  if(local) {
    return;
  } else {
    AliAnalysisGrid *alienHandler = CreateAlienHandler(runPeriod,run_mode,isJDL,type,onMixingAnalysis);
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }
  
  if(!mgr->InitAnalysis()) {
    return;
  }
  
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);

  mgr->StartAnalysis("grid");
  
}

AliAnalysisGrid* CreateAlienHandler(TString runPeriod, TString run_mode, Bool_t isJDL, TString type, bool onMixingAnalysis){
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  plugin->SetMergeViaJDL(isJDL);
  if(run_mode=="terminate")
    plugin->SetRunMode("terminate");
  if(run_mode=="full")
    plugin->SetRunMode("full");
  if(run_mode=="test")
    plugin->SetRunMode("test");
  if(run_mode=="offline")
    plugin->SetRunMode("offline");
  
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20220401-1");
  plugin->SetDefaultOutputs(kTRUE);
  
  if (onMixingAnalysis) {
    plugin->SetGridWorkingDir("PWGLF/AOD/"+runPeriod+"/LowMassMidDimuonMix/"+type);
  } else {
    plugin->SetGridWorkingDir("PWGLF/AOD/"+runPeriod+"/LowMassMidDimuon/"+type);
  }
  plugin->SetGridOutputDir("output");

  plugin->SetSplitMaxInputFileNumber(30);    
  plugin->SetNrunsPerMaster(45);
  
  if(type == "data"){
    plugin->SetRunPrefix("000");
    if (runPeriod.Contains("LHC18c")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18c");
      if (runPeriod.Contains("FAST")) {
	plugin->SetDataPattern("/pass2_FAST/AOD/ AliAOD.root");
      } else {
	plugin->SetDataPattern("/pass2_CENT/AOD/ AliAOD.root");
      }
    }
  } else {
    if (runPeriod.Contains("LHC18c")){
      if (runPeriod.Contains("FAST")) {
	if (runPeriod.Contains("EXTRA")){
	  plugin->SetGridDataDir("/alice/sim/2018/LHC18h1_extra");
	} else {
	  plugin->SetGridDataDir("/alice/sim/2018/LHC18h1_fast");
	}	
	plugin->SetDataPattern("/AOD237/ AliAOD.root");
      } else {
	plugin->SetGridDataDir("/alice/sim/2018/LHC18h1_cent");
	plugin->SetDataPattern("/AOD237/ AliAOD.root");
      }
    }
  }
  
  if(runPeriod.Contains("LHC18c")){
    plugin->AddRunNumber(285958); plugin->AddRunNumber(285957); plugin->AddRunNumber(285946); plugin->AddRunNumber(285917);
    plugin->AddRunNumber(285893); plugin->AddRunNumber(285892); plugin->AddRunNumber(285869); plugin->AddRunNumber(285851);
    plugin->AddRunNumber(285830); plugin->AddRunNumber(285812); plugin->AddRunNumber(285811); plugin->AddRunNumber(285810);
    plugin->AddRunNumber(285806); plugin->AddRunNumber(285805); plugin->AddRunNumber(285804); plugin->AddRunNumber(285781);
    plugin->AddRunNumber(285778); plugin->AddRunNumber(285777); plugin->AddRunNumber(285756); plugin->AddRunNumber(285755);
    plugin->AddRunNumber(285753); plugin->AddRunNumber(285722); plugin->AddRunNumber(285698); plugin->AddRunNumber(285666);
    plugin->AddRunNumber(285664); plugin->AddRunNumber(285663); plugin->AddRunNumber(285662); plugin->AddRunNumber(285643);
    plugin->AddRunNumber(285642); plugin->AddRunNumber(285641); plugin->AddRunNumber(285640); plugin->AddRunNumber(285639);
    plugin->AddRunNumber(285603); plugin->AddRunNumber(285602); plugin->AddRunNumber(285601); plugin->AddRunNumber(285599);
    plugin->AddRunNumber(285578); plugin->AddRunNumber(285577); plugin->AddRunNumber(285576); plugin->AddRunNumber(285575);
    plugin->AddRunNumber(285557); plugin->AddRunNumber(285550); plugin->AddRunNumber(285545); plugin->AddRunNumber(285516);
    plugin->AddRunNumber(285515); plugin->AddRunNumber(285497); plugin->AddRunNumber(285496); plugin->AddRunNumber(285481);
    plugin->AddRunNumber(285471);
  }

  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  plugin->SetAnalysisSource("AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODEventStudy.cxx AliAnalysisTaskAODTrackPair.cxx");
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libANALYSISalice.so libANALYSISaliceBase.so libCORRFW.so libOADB.so libCore.so libTree.so libGeom.so libVMC.so libPhysics.so "
			    "AliAnalysisTaskAODTrackPairUtils.h AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODEventStudy.h AliAnalysisTaskAODEventStudy.cxx AliAnalysisTaskAODTrackPair.h AliAnalysisTaskAODTrackPair.cxx");
  
  //Set Job
  plugin->SetExecutableCommand("aliroot -b -q");
  plugin->SetAnalysisMacro("analysis_syano_"+type+"_"+runPeriod+".C");
  plugin->SetExecutable("analysis_syano_"+type+"_"+runPeriod+".sh");
  
  if(type == "data") plugin->SetNtestFiles(1);
  else               plugin->SetNtestFiles(1);
  
  
  
  plugin->SetOutputToRunNo();
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName("analysis_syano_MC.jdl");
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");
  
  return plugin;
}
