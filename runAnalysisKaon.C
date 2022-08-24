#include "AliAnalysisTaskAODTrackPairUtils.h"
#include "AliAnalysisTaskAODEventStudy.h"
#include "AliAnalysisTaskAODTrackPair.h"

AliAnalysisGrid* CreateAlienHandler(string period, string run_mode, Bool_t isJDL, string type,bool onMixingAnalysis);

void runAnalysisKaon(string runPeriod = "LHC16k",
		     string run_mode  = "test",
		     Bool_t isJDL      = true,
		     string type      = "data",
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
  
  //AliAnalysisTaskPIDResponse *pidTask=reinterpret_cast<AliAnalysisTaskPIDResponse*>
  //(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C (%d,%d,%d,%d)",isMC,false,true,1)));    
  
  AliAnalysisTaskPIDResponse *pidtask=reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C()"));
  
  TFile* input = TFile::Open("./DownScale_Run2_CTP.root");
  
  unsigned int offlineTriggerMask;
  offlineTriggerMask = AliVEvent::kINT7;

  float min_vtxz =-10;
  float max_vtxz = 10;
  int min_vtx_cont = 1;
  float min_pair_rap = -0.5;
  float max_pair_rap = 0.5;
  float alpha = 0.2;
  float pangle = 0.998;
  float v0Dca = 0.1; 
  float trackDca = 1.0;
  float min_dlength = 5.0;
  float max_dlength = 100.;
  string period = runPeriod;
  string multi_method="SPDTracklets";
  bool onPURej = true;
  bool onLBcut = true;
  bool onMuEtaCut = true;
  bool onMuThetaAbsCut = true;
  bool onMuMatchAptCut = true;
  bool onMuMatchLptCut = true;
  bool onMuMatchHptCut = false;
  bool onMuChi2Cut = true;
  bool onMuPdcaCut = true;
  bool isSelectEvt = isEventSelection;
  int paircuttype = 0;
  double min_pairtrackptcut = 0.0;
  bool onMixingAnalysis = isMix;
  bool isMidMuonAnalysis = true;
  bool isKaonAnalysis = true;
  bool isK0sAnalysis = false;  
  float min_track_pt = 0.0;
  float max_track_pt = 999.;
  float min_track_eta = -0.8;
  float max_track_eta =  0.8;
  float min_track_p = 0.15;
  float max_track_p = 2.50;
  float min_pion_sigma_tpc = -5;
  float max_pion_sigma_tpc =  5;
  float min_pion_sigma_tof = -5;
  float max_pion_sigma_tof =  5;
  float min_kaon_sigma_tpc = -3;
  float max_kaon_sigma_tpc =  3;
  float min_kaon_sigma_tof = -3;
  float max_kaon_sigma_tof =  3;
  float min_proton_sigma_tpc = -2;
  float max_proton_sigma_tpc =  2;
  float min_proton_sigma_tof = -2;
  float max_proton_sigma_tof =  2;
  float findable = 0.8;
  string dcaxy = "0.0105 + 0.035/pow(x,1.1)";
  float dcaz = 2.0;
  float chi2tpc = 4.;
  float chi2its = 36.;
  int nclusttpc = 70;
  int nclustits = 1;
  
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gInterpreter->LoadMacro("AliAnalysisTaskAODEventStudy.cxx++g");
  gInterpreter->LoadMacro("AliAnalysisTaskAODTrackPair.cxx++g");
  AliAnalysisTaskAODTrackPair *task
    = reinterpret_cast<AliAnalysisTaskAODTrackPair*>
    (gInterpreter->ExecuteMacro(Form("AddTaskAODTrackPair.C(%u, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, \"%s\", \"%s\",%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d,%d,%f,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\"%s\",%f,%f,%f,%d,%d)",
				     offlineTriggerMask,
				     min_vtxz,
				     max_vtxz,
				     min_vtx_cont,
				     min_pair_rap,
				     max_pair_rap,
				     alpha,
				     pangle,
				     v0Dca,
				     trackDca,
				     min_dlength,
				     max_dlength,
				     period.c_str(),
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
				     isMidMuonAnalysis,
				     isKaonAnalysis,
				     isK0sAnalysis,
				     min_track_pt,
				     max_track_pt,
				     min_track_eta,
				     max_track_eta,
				     min_track_p,
				     max_track_p,
				     min_pion_sigma_tpc,
				     max_pion_sigma_tpc,
				     min_pion_sigma_tof,
				     max_pion_sigma_tof,
				     min_kaon_sigma_tpc,
				     max_kaon_sigma_tpc,
				     min_kaon_sigma_tof,
				     max_kaon_sigma_tof,
				     min_proton_sigma_tpc,
				     max_proton_sigma_tpc,
				     min_proton_sigma_tof,
				     max_proton_sigma_tof,
				     findable,
				     dcaxy.c_str(),
				     dcaz,
				     chi2tpc,
				     chi2its,
				     nclusttpc,
				     nclustits
				     ))); 
#else
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskAODEventStudy.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPair.cxx++g");
  gROOT->LoadMacro("AddTaskAODTrackPairMC.C");
  AliAnalysisTaskAODTrackPair* task = AddTaskAODTrackPairMC(
							      offlineTriggerMask,
							      min_vtxz,
							      max_vtxz,
							      min_vtx_cont,
							      min_pair_rap,
							      max_pair_rap,
							      period.c_str(),
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

AliAnalysisGrid* CreateAlienHandler(string runPeriod, string run_mode, Bool_t isJDL, string type, bool onMixingAnalysis){
  
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  
  plugin->SetMergeViaJDL(isJDL);

  if(run_mode=="terminate"){
    plugin->SetRunMode("terminate");
  }
  if(run_mode=="full"){
    plugin->SetRunMode("full");
  }
  if(run_mode=="test"){
    plugin->SetRunMode("test");
  }
  if(run_mode=="offline"){
    plugin->SetRunMode("offline");
  }
    
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20220805_ROOT6-1");
  plugin->SetDefaultOutputs(kTRUE);
  
  if (onMixingAnalysis) {    
    plugin->SetGridWorkingDir(Form("GlueBall_KK_20220824_2/AOD/%s/TrackPairMix/%s",runPeriod.c_str(),type.c_str()));
  } else {
    plugin->SetGridWorkingDir(Form("GlueBall_KK_20220824_2/AOD/%s/TrackPair/%s",runPeriod.c_str(),type.c_str()));
  }
  
  plugin->SetGridOutputDir("output");

  plugin->SetSplitMaxInputFileNumber(15);    
  //plugin->SetSplitMaxInputFileNumber(5);
  //plugin->SetSplitMaxInputFileNumber(45);
  plugin->SetNrunsPerMaster();
  
  //LHC17q, 
  if(type == "data"){
    if (runPeriod.find("LHC17q") != std::string::npos) {
      plugin->SetGridDataDir("/alice/data/2017/LHC17q");
      if(runPeriod.find("woSDD") != std::string::npos){
	plugin->SetDataPattern("/pass1_CENT_woSDD/AOD234/ AliAOD.root");
      } else {
	plugin->SetDataPattern("/pass1_FAST/AOD234/ AliAOD.root");
      }
    } else if (runPeriod.find("LHC17p") != std::string::npos) {
      plugin->SetGridDataDir("/alice/data/2017/LHC17p");
      if(runPeriod.find("woSDD") != std::string::npos){
	plugin->SetDataPattern("/pass2_CENT_woSDD/AOD/ AliAOD.root");
      } else if(runPeriod.find("FAST") != std::string::npos){
	plugin->SetDataPattern("/pass2_FAST/AOD/ AliAOD.root");
      }
    } else if (runPeriod.find("LHC16q") != std::string::npos) {
      plugin->SetGridDataDir("/alice/data/2016/LHC16q");
      if(runPeriod.find("woSDD") != std::string::npos){
	plugin->SetDataPattern("/pass2_CENT_woSDD/AOD244/ AliAOD.root");
      } else if(runPeriod.find("FAST") != std::string::npos){
	plugin->SetDataPattern("/pass2_FAST/AOD244/ AliAOD.root");
      }
    } else if (runPeriod.find("LHC16t") != std::string::npos) {
      plugin->SetGridDataDir("/alice/data/2016/LHC16t");
      if(runPeriod.find("woSDD") != std::string::npos){
	plugin->SetDataPattern("/pass2_CENT_woSDD/AOD244/ AliAOD.root");
      } else if(runPeriod.find("FAST") != std::string::npos){
	plugin->SetDataPattern("/pass2_FAST/AOD244/ AliAOD.root");
      }
    }
    plugin->SetRunPrefix("000");    
  }
  
  if(runPeriod.find("LHC17q") != std::string::npos){//p-p 5 TeV
    //plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282365);
    plugin->AddRunNumber(282366);
    plugin->AddRunNumber(282367);
  }
  
  if(runPeriod.find("LHC17p_FAST") != std::string::npos){//p-p 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282343); plugin->AddRunNumber(282342); plugin->AddRunNumber(282341); plugin->AddRunNumber(282340); plugin->AddRunNumber(282314); plugin->AddRunNumber(282313); 
    plugin->AddRunNumber(282312); plugin->AddRunNumber(282309); plugin->AddRunNumber(282307); plugin->AddRunNumber(282306); plugin->AddRunNumber(282305); plugin->AddRunNumber(282304);
    plugin->AddRunNumber(282303); plugin->AddRunNumber(282302); plugin->AddRunNumber(282247); plugin->AddRunNumber(282230); plugin->AddRunNumber(282229); plugin->AddRunNumber(282227); 
    plugin->AddRunNumber(282224); plugin->AddRunNumber(282206); plugin->AddRunNumber(282189); plugin->AddRunNumber(282147); plugin->AddRunNumber(282146); plugin->AddRunNumber(282127); 
    plugin->AddRunNumber(282126); plugin->AddRunNumber(282125); plugin->AddRunNumber(282123); plugin->AddRunNumber(282122); plugin->AddRunNumber(282120); plugin->AddRunNumber(282119); 
    plugin->AddRunNumber(282118); plugin->AddRunNumber(282099); plugin->AddRunNumber(282098); plugin->AddRunNumber(282078); plugin->AddRunNumber(282051); plugin->AddRunNumber(282050);
    plugin->AddRunNumber(282031); plugin->AddRunNumber(282030); plugin->AddRunNumber(282025); plugin->AddRunNumber(282021); plugin->AddRunNumber(282016); plugin->AddRunNumber(282008);
 }

  if(runPeriod.find("LHC17p_woSDD_No1") != std::string::npos){//p-p 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(282343); plugin->AddRunNumber(282342); plugin->AddRunNumber(282341); plugin->AddRunNumber(282340); plugin->AddRunNumber(282314); plugin->AddRunNumber(282313); 
    plugin->AddRunNumber(282312); plugin->AddRunNumber(282309); plugin->AddRunNumber(282307); plugin->AddRunNumber(282306); plugin->AddRunNumber(282305); plugin->AddRunNumber(282304); 
    plugin->AddRunNumber(282303); plugin->AddRunNumber(282302); plugin->AddRunNumber(282247); plugin->AddRunNumber(282230); plugin->AddRunNumber(282229); plugin->AddRunNumber(282227); 
    plugin->AddRunNumber(282224); plugin->AddRunNumber(282206); plugin->AddRunNumber(282189); plugin->AddRunNumber(282147); plugin->AddRunNumber(282146); plugin->AddRunNumber(282127); 
    plugin->AddRunNumber(282126); plugin->AddRunNumber(282125); plugin->AddRunNumber(282123); plugin->AddRunNumber(282122); plugin->AddRunNumber(282120); plugin->AddRunNumber(282119); 
    plugin->AddRunNumber(282118); plugin->AddRunNumber(282098); plugin->AddRunNumber(282078); plugin->AddRunNumber(282051); plugin->AddRunNumber(282050);
    plugin->AddRunNumber(282031); plugin->AddRunNumber(282030); plugin->AddRunNumber(282025); plugin->AddRunNumber(282021); plugin->AddRunNumber(282016); plugin->AddRunNumber(282008);
  }
  if(runPeriod.find("LHC17p_woSDD_No2") != std::string::npos){//p-p 5 TeV
    plugin->AddRunNumber(282099);
  }

  if(runPeriod.find("LHC16q") != std::string::npos){//p-Pb 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(265309); plugin->AddRunNumber(265332); plugin->AddRunNumber(265334); plugin->AddRunNumber(265336); 
    plugin->AddRunNumber(265338); plugin->AddRunNumber(265339); plugin->AddRunNumber(265342); plugin->AddRunNumber(265343); 
    plugin->AddRunNumber(265344); plugin->AddRunNumber(265377); plugin->AddRunNumber(265378); plugin->AddRunNumber(265381);
    plugin->AddRunNumber(265383); plugin->AddRunNumber(265384); plugin->AddRunNumber(265385); plugin->AddRunNumber(265387);
    plugin->AddRunNumber(265388); plugin->AddRunNumber(265419); plugin->AddRunNumber(265420); plugin->AddRunNumber(265421);
    plugin->AddRunNumber(265422); plugin->AddRunNumber(265424); plugin->AddRunNumber(265425); plugin->AddRunNumber(265426);
    plugin->AddRunNumber(265427); plugin->AddRunNumber(265435); plugin->AddRunNumber(265499); plugin->AddRunNumber(265500);
    plugin->AddRunNumber(265501); plugin->AddRunNumber(265521); plugin->AddRunNumber(265525);
  }

  if(runPeriod.find("LHC16t") != std::string::npos){//p-Pb 5 TeV
    plugin->SetNrunsPerMaster(4);
    plugin->AddRunNumber(267166); plugin->AddRunNumber(267165); plugin->AddRunNumber(267164); plugin->AddRunNumber(267163);
  }

  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  plugin->SetAnalysisSource("AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODEventStudy.cxx AliAnalysisTaskAODTrackPair.cxx");
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libAOD.so libANALYSIS.so libANALYSISalice.so libANALYSISaliceBase.so libCORRFW.so libOADB.so libCore.so libTree.so libGeom.so libVMC.so libPhysics.so "
			    "AliAnalysisTaskAODTrackPairUtils.h AliAnalysisTaskAODTrackPairUtils.cxx AliAnalysisTaskAODEventStudy.h AliAnalysisTaskAODEventStudy.cxx AliAnalysisTaskAODTrackPair.h AliAnalysisTaskAODTrackPair.cxx");
  
  //Set Job
  plugin->SetExecutableCommand("aliroot -b -q");
  //Form("analysis_syano_%s_%s.C",type.c_str(),runPeriod.c_str())
  plugin->SetAnalysisMacro(Form("analysis_syano_%s_%s.C",type.c_str(),runPeriod.c_str()));
  plugin->SetExecutable(Form("analysis_syano_%s_%s.sh",type.c_str(),runPeriod.c_str()));
  
  if(type == "data") plugin->SetNtestFiles(1);
  else               plugin->SetNtestFiles(1);
  
  
  
  plugin->SetOutputToRunNo();
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName("analysis_syano_MC.jdl");
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");
  
  return plugin;
}
