#if !defined(__CINT__) || defined(__CLING__)
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/AddTaskPIDResponse.C>
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
//#include <PWGLF/RESONANCES/macros/mini/AddTaskF1500.C>
//#include <./AddTaskF1500.C>
#include <./AddTaskAODEventQA.C>
#endif

AliAnalysisGrid* CreateAlienHandler(TString period, TString run_mode, Bool_t isJDL, TString type);

void runAnalysis(TString runPeriod = "LHC16k",
		 TString run_mode  = "test",
		 Bool_t isJDL      = true,
		 TString type      = "data",
		 Bool_t  local     = false)
{
  
  // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
  gInterpreter->ProcessLine("./");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include ./");
#endif
  
  Bool_t isMC = false;
  
  if(type == "data"){

  }
  else{
    isMC = true;
  }

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  if(local) {
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("AliAOD.root");
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
  }
  else {
    // also specify the include (header) paths on grid        
    AliAnalysisGrid *alienHandler = CreateAlienHandler(runPeriod,run_mode,isJDL,type);
    if (!alienHandler) return;
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
  }
  
  // -----------------------------------------
  //            PHYSICS SELECTION
  // -----------------------------------------
#if !defined (__CINT__) || defined (__CLING__)
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);//isMC,PileupCuts
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  if(run_mode == "full") AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);
#endif

  // -----------------------------------------
  //               MULT SELECTION
  // -----------------------------------------
#if !defined (__CINT__) || defined (__CLING__)
  AliMultSelectionTask *multSelTask=reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C ()"));
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask *multSelTask = AddTaskMultSelection(false);
  multSelTask->SetUseDefaultCalib(true);
  if(isMC) multSelTask->SetUseDefaultMCCalib(kTRUE);
  else     multSelTask->SetUseDefaultCalib(kTRUE);
#endif

  TFile* input = TFile::Open("./DownScale_Run2_CTP.root");
  
  float min_vtxz =-10;
  float max_vtxz = 10;

  float min_pair_rap = -4.0;
  float max_pair_rap = -2.5;

  string multi_method="SPDTracklets";

  bool onPURej = true;
  bool onLBcut = true;

  gROOT->LoadMacro("./AliAnalysisTaskAODTrackPairUtils.cxx++g");
  gROOT->LoadMacro("./AliAnalysisTaskAODEventStudy.cxx++g");
  gROOT->LoadMacro("./AddTaskAODEventQA.C");

  gROOT->LoadMacro("./AliAnalysisTaskAODTrackPair.cxx++g");
  gROOT->LoadMacro("./AddTaskAODTrackPair.C");

  /*
    UInt_t offlineTriggerMask = AliVEvent::kAny,
    float min_vtxz =-10,
    float max_vtxz = 10,
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
  */

  AliAnalysisTaskAODEventStudy* qa = AddTaskAODEventQA(AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kINT7inMUON,-10,10,-4.0,-2.5,"SPDTracklets",1,1,1,1,0,0,0,1,1,isMC);  
  //AliAnalysisTaskAODTrackPair* dimuon = AddTaskAODTrackPair(AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7,-10,10,-4.0,-2.5,"SPDTracklets",1,1,1,1,0,0,0,1,1,isMC);  
  
  // -----------------------------------------
  //               Add Task V0Reader
  // -----------------------------------------  
#if !defined (__CINT__) || defined (__CLING__)
  //AddTaskGlueball();
#else
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddTaskGlueball.C");
  
  //gROOT->LoadMacro("./AddTaskGlueball.C");
  //AliRsnMiniAnalysisTask *task = AddTaskGlueball("gluball",isMC,1,AliPIDResponse::kPP,AliVEvent::kINT7,0,0,0,32,1,AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,1.0,1.0,0.2,2.5,2300,0.0,10.0,20,0.001,10.0,1,1);
  //AddTaskGlueball("glueball",isMC,1,AliPIDResponse::kPP,AliVEvent::kINT7,0,0,0,5,1,AliRsnCutSetDaughterParticle::kTPCpidTOFveto,2.0,2.0,0.2,2.5,2300,0.0,10.0,20,0.001,10.0,1,1);
#endif
  if(!mgr->InitAnalysis()) return;
  
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);

  mgr->StartAnalysis("grid");
  
}

AliAnalysisGrid* CreateAlienHandler(TString runPeriod, TString run_mode, Bool_t isJDL, TString type){
  
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
  plugin->SetAliPhysicsVersion("vAN-20210822-1");
  plugin->SetDefaultOutputs(kTRUE);
  
  plugin->SetGridWorkingDir("PWGLF/AOD/"+runPeriod+"/LowMassDimuon/"+type);
  plugin->SetGridOutputDir("output");
  
  if(type == "data"){
    if (runPeriod.Contains("LHC16h")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16h");    
      plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
    }
    if (runPeriod.Contains("LHC16i")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16i");    
      plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC16j")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16j");
      plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC16k")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16k");
      plugin->SetDataPattern("/pass3/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC16l")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16l");
      plugin->SetDataPattern("/pass3/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC16o")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16o");
      plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC16p")){
      plugin->SetGridDataDir("/alice/data/2016/LHC16p");
      plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
    }

    else if (runPeriod.Contains("LHC17h")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17h");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17i")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17i");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17k")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17k");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17l")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17l");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17m")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17m");
      //plugin->SetDataPattern("/pass1/AOD234/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17o")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17o");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD203/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC17r")){
      plugin->SetGridDataDir("/alice/data/2017/LHC17r");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18c")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18c");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18d")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18d");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18e")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18e");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18f")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18f");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18l")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18l");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18m")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18m");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass1/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18o")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18o");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass2/AOD/ AliAOD.root");
    }
    else if (runPeriod.Contains("LHC18p")){
      plugin->SetGridDataDir("/alice/data/2018/LHC18p");
      //plugin->SetDataPattern("/pass2/AOD/ AliAOD.root");
      plugin->SetDataPattern("/muon_calo_pass2/AOD/ AliAOD.root");
    }
    
    plugin->SetRunPrefix("000");    
    //plugin->SetNrunsPerMaster(1);
    plugin->SetSplitMaxInputFileNumber(30);    
    plugin->SetNrunsPerMaster(15);
  }

  if(runPeriod.Contains("LHC16h")){
    plugin->AddRunNumber(255467); plugin->AddRunNumber(255466); plugin->AddRunNumber(255465); plugin->AddRunNumber(255463); plugin->AddRunNumber(255447); plugin->AddRunNumber(255442); plugin->AddRunNumber(255440); 
    plugin->AddRunNumber(255415); plugin->AddRunNumber(255402); plugin->AddRunNumber(255398); plugin->AddRunNumber(255352); plugin->AddRunNumber(255351); plugin->AddRunNumber(255350); plugin->AddRunNumber(255283); 
    plugin->AddRunNumber(255280); plugin->AddRunNumber(255276); plugin->AddRunNumber(255275); plugin->AddRunNumber(255256); plugin->AddRunNumber(255255); plugin->AddRunNumber(255253); plugin->AddRunNumber(255252);
    plugin->AddRunNumber(255251); plugin->AddRunNumber(255249); plugin->AddRunNumber(255248); plugin->AddRunNumber(255247); plugin->AddRunNumber(255242); plugin->AddRunNumber(255240); plugin->AddRunNumber(255182);
    plugin->AddRunNumber(255180); plugin->AddRunNumber(255177); plugin->AddRunNumber(255176); plugin->AddRunNumber(255173); plugin->AddRunNumber(255171); plugin->AddRunNumber(255167); plugin->AddRunNumber(255162);
    plugin->AddRunNumber(255159); plugin->AddRunNumber(255154); plugin->AddRunNumber(255111); plugin->AddRunNumber(255091); plugin->AddRunNumber(255086); plugin->AddRunNumber(255085); plugin->AddRunNumber(255082);
    plugin->AddRunNumber(255079); plugin->AddRunNumber(255076); plugin->AddRunNumber(255075); plugin->AddRunNumber(255074); plugin->AddRunNumber(255073); plugin->AddRunNumber(255071); plugin->AddRunNumber(255068);
    plugin->AddRunNumber(255042); plugin->AddRunNumber(255010); plugin->AddRunNumber(255009); plugin->AddRunNumber(255008); plugin->AddRunNumber(254984); plugin->AddRunNumber(254983); plugin->AddRunNumber(254654);
    plugin->AddRunNumber(254653); plugin->AddRunNumber(254652); plugin->AddRunNumber(254651); plugin->AddRunNumber(254649); plugin->AddRunNumber(254648); plugin->AddRunNumber(254646); plugin->AddRunNumber(254644);
    plugin->AddRunNumber(254640); plugin->AddRunNumber(254632); plugin->AddRunNumber(254630); plugin->AddRunNumber(254629); plugin->AddRunNumber(254621); plugin->AddRunNumber(254608); plugin->AddRunNumber(254606);
    plugin->AddRunNumber(254604); plugin->AddRunNumber(254419);
  }
  if(runPeriod.Contains("LHC16j")){
    plugin->AddRunNumber(256420); plugin->AddRunNumber(256418); plugin->AddRunNumber(256417); plugin->AddRunNumber(256415); plugin->AddRunNumber(256373); plugin->AddRunNumber(256372); plugin->AddRunNumber(256371);
    plugin->AddRunNumber(256368); plugin->AddRunNumber(256366); plugin->AddRunNumber(256365); plugin->AddRunNumber(256364); plugin->AddRunNumber(256363); plugin->AddRunNumber(256362); plugin->AddRunNumber(256361);
    plugin->AddRunNumber(256356); plugin->AddRunNumber(256311); plugin->AddRunNumber(256307); plugin->AddRunNumber(256302); plugin->AddRunNumber(256298); plugin->AddRunNumber(256297); plugin->AddRunNumber(256295);
    plugin->AddRunNumber(256292); plugin->AddRunNumber(256290); plugin->AddRunNumber(256289); plugin->AddRunNumber(256287); plugin->AddRunNumber(256284); plugin->AddRunNumber(256283); plugin->AddRunNumber(256282);
    plugin->AddRunNumber(256281); plugin->AddRunNumber(256231); plugin->AddRunNumber(256228); plugin->AddRunNumber(256227); plugin->AddRunNumber(256223); plugin->AddRunNumber(256222); plugin->AddRunNumber(256219);
    plugin->AddRunNumber(256215); plugin->AddRunNumber(256213); plugin->AddRunNumber(256212); plugin->AddRunNumber(256210); plugin->AddRunNumber(256204); plugin->AddRunNumber(256169); plugin->AddRunNumber(256161);
    plugin->AddRunNumber(256158); plugin->AddRunNumber(256157); plugin->AddRunNumber(256156); plugin->AddRunNumber(256149); plugin->AddRunNumber(256148); plugin->AddRunNumber(256147); plugin->AddRunNumber(256146);
  }
  if(runPeriod.Contains("LHC16k")){
    plugin->AddRunNumber(258537); plugin->AddRunNumber(258499); plugin->AddRunNumber(258498); plugin->AddRunNumber(258477); plugin->AddRunNumber(258456); plugin->AddRunNumber(258454); plugin->AddRunNumber(258452);
    plugin->AddRunNumber(258426); plugin->AddRunNumber(258399); plugin->AddRunNumber(258393); plugin->AddRunNumber(258391); plugin->AddRunNumber(258388); plugin->AddRunNumber(258387); plugin->AddRunNumber(258359);
    plugin->AddRunNumber(258336); plugin->AddRunNumber(258332); plugin->AddRunNumber(258307); plugin->AddRunNumber(258306); plugin->AddRunNumber(258303); plugin->AddRunNumber(258302); plugin->AddRunNumber(258301);
    plugin->AddRunNumber(258299); plugin->AddRunNumber(258280); plugin->AddRunNumber(258278); plugin->AddRunNumber(258274); plugin->AddRunNumber(258273); plugin->AddRunNumber(258271); plugin->AddRunNumber(258270);
    plugin->AddRunNumber(258258); plugin->AddRunNumber(258257); plugin->AddRunNumber(258256); plugin->AddRunNumber(258204); plugin->AddRunNumber(258203); plugin->AddRunNumber(258202); plugin->AddRunNumber(258197);
    plugin->AddRunNumber(258178); plugin->AddRunNumber(258117); plugin->AddRunNumber(258114); plugin->AddRunNumber(258113); plugin->AddRunNumber(258109); plugin->AddRunNumber(258108); plugin->AddRunNumber(258107);
    plugin->AddRunNumber(258063); plugin->AddRunNumber(258062); plugin->AddRunNumber(258060); plugin->AddRunNumber(258059); plugin->AddRunNumber(258049); plugin->AddRunNumber(258048); plugin->AddRunNumber(258045);
    plugin->AddRunNumber(258042); plugin->AddRunNumber(258041); plugin->AddRunNumber(258039); plugin->AddRunNumber(258019); plugin->AddRunNumber(258017); plugin->AddRunNumber(258014); plugin->AddRunNumber(258012);
    plugin->AddRunNumber(258008); plugin->AddRunNumber(257989); plugin->AddRunNumber(257986); plugin->AddRunNumber(257979); plugin->AddRunNumber(257963); plugin->AddRunNumber(257960); plugin->AddRunNumber(257958);
    plugin->AddRunNumber(257957); plugin->AddRunNumber(257939); plugin->AddRunNumber(257937); plugin->AddRunNumber(257936); plugin->AddRunNumber(257932); plugin->AddRunNumber(257912); plugin->AddRunNumber(257901);
    plugin->AddRunNumber(257893); plugin->AddRunNumber(257892); plugin->AddRunNumber(257737); plugin->AddRunNumber(257735); plugin->AddRunNumber(257734); plugin->AddRunNumber(257733); plugin->AddRunNumber(257727);
    plugin->AddRunNumber(257725); plugin->AddRunNumber(257724); plugin->AddRunNumber(257697); plugin->AddRunNumber(257694); plugin->AddRunNumber(257688); plugin->AddRunNumber(257687); plugin->AddRunNumber(257685);
    plugin->AddRunNumber(257684); plugin->AddRunNumber(257682); plugin->AddRunNumber(257644); plugin->AddRunNumber(257642); plugin->AddRunNumber(257636); plugin->AddRunNumber(257635); plugin->AddRunNumber(257632);
    plugin->AddRunNumber(257630); plugin->AddRunNumber(257606); plugin->AddRunNumber(257605); plugin->AddRunNumber(257604); plugin->AddRunNumber(257601); plugin->AddRunNumber(257595); plugin->AddRunNumber(257594);
    plugin->AddRunNumber(257592); plugin->AddRunNumber(257590); plugin->AddRunNumber(257588); plugin->AddRunNumber(257587); plugin->AddRunNumber(257566); plugin->AddRunNumber(257565); plugin->AddRunNumber(257564);
    plugin->AddRunNumber(257563); plugin->AddRunNumber(257562); plugin->AddRunNumber(257561); plugin->AddRunNumber(257560); plugin->AddRunNumber(257541); plugin->AddRunNumber(257540); plugin->AddRunNumber(257531);
    plugin->AddRunNumber(257530); plugin->AddRunNumber(257492); plugin->AddRunNumber(257491); plugin->AddRunNumber(257490); plugin->AddRunNumber(257488); plugin->AddRunNumber(257487); plugin->AddRunNumber(257474);
    plugin->AddRunNumber(257468); plugin->AddRunNumber(257457); plugin->AddRunNumber(257433); plugin->AddRunNumber(257364); plugin->AddRunNumber(257358); plugin->AddRunNumber(257330); plugin->AddRunNumber(257322);
    plugin->AddRunNumber(257320); plugin->AddRunNumber(257318); plugin->AddRunNumber(257260); plugin->AddRunNumber(257224); plugin->AddRunNumber(257095); plugin->AddRunNumber(257092); plugin->AddRunNumber(257086);
    plugin->AddRunNumber(257084); plugin->AddRunNumber(257083); plugin->AddRunNumber(257082); plugin->AddRunNumber(257080); plugin->AddRunNumber(257077); plugin->AddRunNumber(257071); plugin->AddRunNumber(257026);
    plugin->AddRunNumber(257021); plugin->AddRunNumber(257012); plugin->AddRunNumber(257011); plugin->AddRunNumber(256944); plugin->AddRunNumber(256942); plugin->AddRunNumber(256941); plugin->AddRunNumber(256697);
    plugin->AddRunNumber(256695); plugin->AddRunNumber(256694); plugin->AddRunNumber(256691); plugin->AddRunNumber(256684); plugin->AddRunNumber(256681); plugin->AddRunNumber(256677); plugin->AddRunNumber(256676);
    plugin->AddRunNumber(256658); plugin->AddRunNumber(256620); plugin->AddRunNumber(256619); plugin->AddRunNumber(256591); plugin->AddRunNumber(256567); plugin->AddRunNumber(256565); plugin->AddRunNumber(256564);
    plugin->AddRunNumber(256561); plugin->AddRunNumber(256560); plugin->AddRunNumber(256557); plugin->AddRunNumber(256556); plugin->AddRunNumber(256554); plugin->AddRunNumber(256552); plugin->AddRunNumber(256512);
    plugin->AddRunNumber(256510); plugin->AddRunNumber(256506); plugin->AddRunNumber(256504);
  }
  if(runPeriod.Contains("LHC16o")){
    plugin->AddRunNumber(264035); plugin->AddRunNumber(264033); plugin->AddRunNumber(263985); plugin->AddRunNumber(263984); plugin->AddRunNumber(263981); plugin->AddRunNumber(263979); plugin->AddRunNumber(263978);
    plugin->AddRunNumber(263977); plugin->AddRunNumber(263923); plugin->AddRunNumber(263920); plugin->AddRunNumber(263917); plugin->AddRunNumber(263916); plugin->AddRunNumber(263905); plugin->AddRunNumber(263866);
    plugin->AddRunNumber(263863); plugin->AddRunNumber(263861); plugin->AddRunNumber(263830); plugin->AddRunNumber(263829); plugin->AddRunNumber(263824); plugin->AddRunNumber(263823); plugin->AddRunNumber(263813);
    plugin->AddRunNumber(263810); plugin->AddRunNumber(263803); plugin->AddRunNumber(263793); plugin->AddRunNumber(263792); plugin->AddRunNumber(263790); plugin->AddRunNumber(263787); plugin->AddRunNumber(263786);
    plugin->AddRunNumber(263785); plugin->AddRunNumber(263784); plugin->AddRunNumber(263744); plugin->AddRunNumber(263743); plugin->AddRunNumber(263741); plugin->AddRunNumber(263739); plugin->AddRunNumber(263738);
    plugin->AddRunNumber(263737); plugin->AddRunNumber(263691); plugin->AddRunNumber(263690); plugin->AddRunNumber(263689); plugin->AddRunNumber(263682); plugin->AddRunNumber(263662); plugin->AddRunNumber(263657);
    plugin->AddRunNumber(263654); plugin->AddRunNumber(263653); plugin->AddRunNumber(263652); plugin->AddRunNumber(263647); plugin->AddRunNumber(263529); plugin->AddRunNumber(263497); plugin->AddRunNumber(263496);
    plugin->AddRunNumber(263490); plugin->AddRunNumber(263487); plugin->AddRunNumber(263332); plugin->AddRunNumber(262858); plugin->AddRunNumber(262855); plugin->AddRunNumber(262853); plugin->AddRunNumber(262849);
    plugin->AddRunNumber(262847); plugin->AddRunNumber(262844); plugin->AddRunNumber(262842); plugin->AddRunNumber(262841); plugin->AddRunNumber(262778); plugin->AddRunNumber(262777); plugin->AddRunNumber(262776);
    plugin->AddRunNumber(262768); plugin->AddRunNumber(262760); plugin->AddRunNumber(262727); plugin->AddRunNumber(262725); plugin->AddRunNumber(262723); plugin->AddRunNumber(262719); plugin->AddRunNumber(262717);
    plugin->AddRunNumber(262713); plugin->AddRunNumber(262705); plugin->AddRunNumber(262635); plugin->AddRunNumber(262632); plugin->AddRunNumber(262628); plugin->AddRunNumber(262594); plugin->AddRunNumber(262593);
    plugin->AddRunNumber(262583); plugin->AddRunNumber(262578); plugin->AddRunNumber(262574); plugin->AddRunNumber(262572); plugin->AddRunNumber(262571); plugin->AddRunNumber(262570); plugin->AddRunNumber(262569);
    plugin->AddRunNumber(262568); plugin->AddRunNumber(262567); plugin->AddRunNumber(262563); plugin->AddRunNumber(262537); plugin->AddRunNumber(262533); plugin->AddRunNumber(262532); plugin->AddRunNumber(262528);
    plugin->AddRunNumber(262492); plugin->AddRunNumber(262487); plugin->AddRunNumber(262451); plugin->AddRunNumber(262430); plugin->AddRunNumber(262428); plugin->AddRunNumber(262424); plugin->AddRunNumber(262423);
    plugin->AddRunNumber(262422); plugin->AddRunNumber(262419); plugin->AddRunNumber(262418);
  }
  if(runPeriod.Contains("LHC16p")){
    plugin->AddRunNumber(264347); plugin->AddRunNumber(264346); plugin->AddRunNumber(264345); plugin->AddRunNumber(264341); plugin->AddRunNumber(264336); plugin->AddRunNumber(264312); plugin->AddRunNumber(264305);
    plugin->AddRunNumber(264281); plugin->AddRunNumber(264279); plugin->AddRunNumber(264277); plugin->AddRunNumber(264273); plugin->AddRunNumber(264267); plugin->AddRunNumber(264266); plugin->AddRunNumber(264265);
    plugin->AddRunNumber(264264); plugin->AddRunNumber(264262); plugin->AddRunNumber(264261); plugin->AddRunNumber(264260); plugin->AddRunNumber(264259); plugin->AddRunNumber(264238); plugin->AddRunNumber(264233);
    plugin->AddRunNumber(264232); plugin->AddRunNumber(264198); plugin->AddRunNumber(264197); plugin->AddRunNumber(264194); plugin->AddRunNumber(264188); plugin->AddRunNumber(264168); plugin->AddRunNumber(264164);
    plugin->AddRunNumber(264138); plugin->AddRunNumber(264137); plugin->AddRunNumber(264129); plugin->AddRunNumber(264110); plugin->AddRunNumber(264109); plugin->AddRunNumber(264086); plugin->AddRunNumber(264085);
    plugin->AddRunNumber(264082); plugin->AddRunNumber(264078); plugin->AddRunNumber(264076);
  }
  
  if(runPeriod.Contains("LHC17h")){
    plugin->AddRunNumber(273103); plugin->AddRunNumber(273101); plugin->AddRunNumber(273100); plugin->AddRunNumber(273099); plugin->AddRunNumber(273077); plugin->AddRunNumber(273010); plugin->AddRunNumber(273009);
    plugin->AddRunNumber(272985); plugin->AddRunNumber(272983); plugin->AddRunNumber(272976); plugin->AddRunNumber(272949); plugin->AddRunNumber(272947); plugin->AddRunNumber(272939); plugin->AddRunNumber(272935);
    plugin->AddRunNumber(272934); plugin->AddRunNumber(272933); plugin->AddRunNumber(272932); plugin->AddRunNumber(272905); plugin->AddRunNumber(272903); plugin->AddRunNumber(272880); plugin->AddRunNumber(272873);
    plugin->AddRunNumber(272871); plugin->AddRunNumber(272870); plugin->AddRunNumber(272836); plugin->AddRunNumber(272835); plugin->AddRunNumber(272834); plugin->AddRunNumber(272833); plugin->AddRunNumber(272829);
    plugin->AddRunNumber(272828); plugin->AddRunNumber(272784); plugin->AddRunNumber(272783); plugin->AddRunNumber(272782); plugin->AddRunNumber(272762); plugin->AddRunNumber(272760); plugin->AddRunNumber(272749);
    plugin->AddRunNumber(272747); plugin->AddRunNumber(272746); plugin->AddRunNumber(272692); plugin->AddRunNumber(272691); plugin->AddRunNumber(272620); plugin->AddRunNumber(272619); plugin->AddRunNumber(272608);
    plugin->AddRunNumber(272607); plugin->AddRunNumber(272585); plugin->AddRunNumber(272577); plugin->AddRunNumber(272575); plugin->AddRunNumber(272574); plugin->AddRunNumber(272521); plugin->AddRunNumber(272469);
    plugin->AddRunNumber(272468); plugin->AddRunNumber(272466); plugin->AddRunNumber(272463); plugin->AddRunNumber(272462); plugin->AddRunNumber(272461); plugin->AddRunNumber(272414); plugin->AddRunNumber(272413);
    plugin->AddRunNumber(272411); plugin->AddRunNumber(272400); plugin->AddRunNumber(272394); plugin->AddRunNumber(272360); plugin->AddRunNumber(272359); plugin->AddRunNumber(272335); plugin->AddRunNumber(272194);
    plugin->AddRunNumber(272156); plugin->AddRunNumber(272155); plugin->AddRunNumber(272154); plugin->AddRunNumber(272153); plugin->AddRunNumber(272152); plugin->AddRunNumber(272151); plugin->AddRunNumber(272123);
    plugin->AddRunNumber(272101); plugin->AddRunNumber(272100); plugin->AddRunNumber(272076); plugin->AddRunNumber(272075); plugin->AddRunNumber(272042); plugin->AddRunNumber(272041); plugin->AddRunNumber(272040);
    plugin->AddRunNumber(272039); plugin->AddRunNumber(272038); plugin->AddRunNumber(272036); plugin->AddRunNumber(272034); plugin->AddRunNumber(272030); plugin->AddRunNumber(272029); plugin->AddRunNumber(272025);
    plugin->AddRunNumber(272020); plugin->AddRunNumber(271970); plugin->AddRunNumber(271969); plugin->AddRunNumber(271962); plugin->AddRunNumber(271955); plugin->AddRunNumber(271953); plugin->AddRunNumber(271946);
    plugin->AddRunNumber(271925); plugin->AddRunNumber(271921); plugin->AddRunNumber(271915); plugin->AddRunNumber(271912); plugin->AddRunNumber(271886); plugin->AddRunNumber(271879); plugin->AddRunNumber(271878);
    plugin->AddRunNumber(271874); plugin->AddRunNumber(271873); plugin->AddRunNumber(271871); plugin->AddRunNumber(271870); plugin->AddRunNumber(271868);
  }
  if(runPeriod.Contains("LHC17i")){
    plugin->AddRunNumber(274442); plugin->AddRunNumber(274390); plugin->AddRunNumber(274387); plugin->AddRunNumber(274385); plugin->AddRunNumber(274364); plugin->AddRunNumber(274363); plugin->AddRunNumber(274360);
    plugin->AddRunNumber(274357); plugin->AddRunNumber(274355); plugin->AddRunNumber(274329); plugin->AddRunNumber(274283); plugin->AddRunNumber(274281); plugin->AddRunNumber(274280); plugin->AddRunNumber(274278);
    plugin->AddRunNumber(274276); plugin->AddRunNumber(274271); plugin->AddRunNumber(274270); plugin->AddRunNumber(274269); plugin->AddRunNumber(274268); plugin->AddRunNumber(274266); plugin->AddRunNumber(274264);
    plugin->AddRunNumber(274263); plugin->AddRunNumber(274259); plugin->AddRunNumber(274232); plugin->AddRunNumber(274212); plugin->AddRunNumber(274148); plugin->AddRunNumber(274147); plugin->AddRunNumber(274125);
    plugin->AddRunNumber(274094); plugin->AddRunNumber(274092); plugin->AddRunNumber(274064); plugin->AddRunNumber(274063); plugin->AddRunNumber(274058); plugin->AddRunNumber(273986); plugin->AddRunNumber(273985);
    plugin->AddRunNumber(273946); plugin->AddRunNumber(273942); plugin->AddRunNumber(273918); plugin->AddRunNumber(273889); plugin->AddRunNumber(273887); plugin->AddRunNumber(273886); plugin->AddRunNumber(273885);
    plugin->AddRunNumber(273825); plugin->AddRunNumber(273824); plugin->AddRunNumber(273719); plugin->AddRunNumber(273711); plugin->AddRunNumber(273709); plugin->AddRunNumber(273695); plugin->AddRunNumber(273690);
    plugin->AddRunNumber(273689); plugin->AddRunNumber(273687); plugin->AddRunNumber(273654); plugin->AddRunNumber(273653); plugin->AddRunNumber(273593); plugin->AddRunNumber(273592); plugin->AddRunNumber(273591);
  }
  if(runPeriod.Contains("LHC17k")){
    plugin->AddRunNumber(276508); plugin->AddRunNumber(276507); plugin->AddRunNumber(276506); plugin->AddRunNumber(276500); plugin->AddRunNumber(276462); plugin->AddRunNumber(276461); plugin->AddRunNumber(276439);
    plugin->AddRunNumber(276438); plugin->AddRunNumber(276437); plugin->AddRunNumber(276435); plugin->AddRunNumber(276434); plugin->AddRunNumber(276432); plugin->AddRunNumber(276429); plugin->AddRunNumber(276351);
    plugin->AddRunNumber(276348); plugin->AddRunNumber(276302); plugin->AddRunNumber(276297); plugin->AddRunNumber(276294); plugin->AddRunNumber(276292); plugin->AddRunNumber(276291); plugin->AddRunNumber(276290);
    plugin->AddRunNumber(276259); plugin->AddRunNumber(276230); plugin->AddRunNumber(276205); plugin->AddRunNumber(276178); plugin->AddRunNumber(276177); plugin->AddRunNumber(276170); plugin->AddRunNumber(276169);
    plugin->AddRunNumber(276166); plugin->AddRunNumber(276145); plugin->AddRunNumber(276141); plugin->AddRunNumber(276140); plugin->AddRunNumber(276108); plugin->AddRunNumber(276105); plugin->AddRunNumber(276104);
    plugin->AddRunNumber(276102); plugin->AddRunNumber(276099); plugin->AddRunNumber(276098); plugin->AddRunNumber(275664); plugin->AddRunNumber(275661); plugin->AddRunNumber(275657); plugin->AddRunNumber(275650);
    plugin->AddRunNumber(275648); plugin->AddRunNumber(275624); plugin->AddRunNumber(275559); plugin->AddRunNumber(275558); plugin->AddRunNumber(275515); plugin->AddRunNumber(275472); plugin->AddRunNumber(275471);
    plugin->AddRunNumber(275467); plugin->AddRunNumber(275459); plugin->AddRunNumber(275457); plugin->AddRunNumber(275453); plugin->AddRunNumber(275452); plugin->AddRunNumber(275448); plugin->AddRunNumber(275406);
    plugin->AddRunNumber(275404); plugin->AddRunNumber(275401); plugin->AddRunNumber(275369); plugin->AddRunNumber(275361); plugin->AddRunNumber(275360); plugin->AddRunNumber(275357); plugin->AddRunNumber(275332);
    plugin->AddRunNumber(275328); plugin->AddRunNumber(275283); plugin->AddRunNumber(275247); plugin->AddRunNumber(275246); plugin->AddRunNumber(275245); plugin->AddRunNumber(275188); plugin->AddRunNumber(275177);
    plugin->AddRunNumber(275175); plugin->AddRunNumber(275174); plugin->AddRunNumber(275173); plugin->AddRunNumber(275151); plugin->AddRunNumber(275150); plugin->AddRunNumber(275149); plugin->AddRunNumber(275076);
    plugin->AddRunNumber(275075); plugin->AddRunNumber(275073); plugin->AddRunNumber(275070); plugin->AddRunNumber(275068); plugin->AddRunNumber(275067); plugin->AddRunNumber(274979); plugin->AddRunNumber(274978);
    plugin->AddRunNumber(274886); plugin->AddRunNumber(274884); plugin->AddRunNumber(274883); plugin->AddRunNumber(274882); plugin->AddRunNumber(274822); plugin->AddRunNumber(274817); plugin->AddRunNumber(274815);
    plugin->AddRunNumber(274811); plugin->AddRunNumber(274807); plugin->AddRunNumber(274806); plugin->AddRunNumber(274803); plugin->AddRunNumber(274802); plugin->AddRunNumber(274801); plugin->AddRunNumber(274743);
    plugin->AddRunNumber(274736); plugin->AddRunNumber(274708);
  }
  if(runPeriod.Contains("LHC17l")){
    plugin->AddRunNumber(278216); plugin->AddRunNumber(278215); plugin->AddRunNumber(278191); plugin->AddRunNumber(278189); plugin->AddRunNumber(278167); plugin->AddRunNumber(278166); plugin->AddRunNumber(278165);
    plugin->AddRunNumber(278164); plugin->AddRunNumber(278163); plugin->AddRunNumber(278130); plugin->AddRunNumber(278127); plugin->AddRunNumber(278126); plugin->AddRunNumber(278123); plugin->AddRunNumber(278122);
    plugin->AddRunNumber(278121); plugin->AddRunNumber(277996); plugin->AddRunNumber(277991); plugin->AddRunNumber(277989); plugin->AddRunNumber(277988); plugin->AddRunNumber(277987); plugin->AddRunNumber(277952);
    plugin->AddRunNumber(277930); plugin->AddRunNumber(277907); plugin->AddRunNumber(277904); plugin->AddRunNumber(277903); plugin->AddRunNumber(277901); plugin->AddRunNumber(277900); plugin->AddRunNumber(277899);
    plugin->AddRunNumber(277898); plugin->AddRunNumber(277897); plugin->AddRunNumber(277876); plugin->AddRunNumber(277870); plugin->AddRunNumber(277848); plugin->AddRunNumber(277847); plugin->AddRunNumber(277842);
    plugin->AddRunNumber(277841); plugin->AddRunNumber(277836); plugin->AddRunNumber(277834); plugin->AddRunNumber(277801); plugin->AddRunNumber(277800); plugin->AddRunNumber(277799); plugin->AddRunNumber(277795);
    plugin->AddRunNumber(277794); plugin->AddRunNumber(277749); plugin->AddRunNumber(277747); plugin->AddRunNumber(277746); plugin->AddRunNumber(277725); plugin->AddRunNumber(277577); plugin->AddRunNumber(277576);
    plugin->AddRunNumber(277575); plugin->AddRunNumber(277574); plugin->AddRunNumber(277537); plugin->AddRunNumber(277536); plugin->AddRunNumber(277531); plugin->AddRunNumber(277530); plugin->AddRunNumber(277479);
    plugin->AddRunNumber(277478); plugin->AddRunNumber(277476); plugin->AddRunNumber(277473); plugin->AddRunNumber(277472); plugin->AddRunNumber(277470); plugin->AddRunNumber(277418); plugin->AddRunNumber(277417);
    plugin->AddRunNumber(277389); plugin->AddRunNumber(277386); plugin->AddRunNumber(277384); plugin->AddRunNumber(277383); plugin->AddRunNumber(277360); plugin->AddRunNumber(277314); plugin->AddRunNumber(277312);
    plugin->AddRunNumber(277310); plugin->AddRunNumber(277293); plugin->AddRunNumber(277262); plugin->AddRunNumber(277256); plugin->AddRunNumber(277197); plugin->AddRunNumber(277196); plugin->AddRunNumber(277194);
    plugin->AddRunNumber(277193); plugin->AddRunNumber(277189); plugin->AddRunNumber(277188); plugin->AddRunNumber(277184); plugin->AddRunNumber(277183); plugin->AddRunNumber(277182); plugin->AddRunNumber(277181);
    plugin->AddRunNumber(277180); plugin->AddRunNumber(277155); plugin->AddRunNumber(277121); plugin->AddRunNumber(277117); plugin->AddRunNumber(277091); plugin->AddRunNumber(277087); plugin->AddRunNumber(277082);
    plugin->AddRunNumber(277079); plugin->AddRunNumber(277076); plugin->AddRunNumber(277073); plugin->AddRunNumber(277037); plugin->AddRunNumber(277017); plugin->AddRunNumber(277016); plugin->AddRunNumber(277015);
    plugin->AddRunNumber(276972); plugin->AddRunNumber(276971); plugin->AddRunNumber(276970); plugin->AddRunNumber(276969); plugin->AddRunNumber(276920); plugin->AddRunNumber(276917); plugin->AddRunNumber(276916);
    plugin->AddRunNumber(276762); plugin->AddRunNumber(276675); plugin->AddRunNumber(276674); plugin->AddRunNumber(276672); plugin->AddRunNumber(276671); plugin->AddRunNumber(276670); plugin->AddRunNumber(276669);
    plugin->AddRunNumber(276644); plugin->AddRunNumber(276608); plugin->AddRunNumber(276557); plugin->AddRunNumber(276553); plugin->AddRunNumber(276552); plugin->AddRunNumber(276551);

  }
  if(runPeriod.Contains("LHC17m")){
    plugin->AddRunNumber(280140); plugin->AddRunNumber(280135); plugin->AddRunNumber(280134); plugin->AddRunNumber(280131); plugin->AddRunNumber(280126); plugin->AddRunNumber(280118); plugin->AddRunNumber(280114);
    plugin->AddRunNumber(280111); plugin->AddRunNumber(280108); plugin->AddRunNumber(280066); plugin->AddRunNumber(280052); plugin->AddRunNumber(280051); plugin->AddRunNumber(280049); plugin->AddRunNumber(279955);
    plugin->AddRunNumber(279954); plugin->AddRunNumber(279952); plugin->AddRunNumber(279893); plugin->AddRunNumber(279890); plugin->AddRunNumber(279886); plugin->AddRunNumber(279884); plugin->AddRunNumber(279880);
    plugin->AddRunNumber(279879); plugin->AddRunNumber(279855); plugin->AddRunNumber(279854); plugin->AddRunNumber(279853); plugin->AddRunNumber(279830); plugin->AddRunNumber(279827); plugin->AddRunNumber(279826);
    plugin->AddRunNumber(279773); plugin->AddRunNumber(279749); plugin->AddRunNumber(279747); plugin->AddRunNumber(279719); plugin->AddRunNumber(279718); plugin->AddRunNumber(279715); plugin->AddRunNumber(279689);
    plugin->AddRunNumber(279688); plugin->AddRunNumber(279684); plugin->AddRunNumber(279683); plugin->AddRunNumber(279682); plugin->AddRunNumber(279679); plugin->AddRunNumber(279677); plugin->AddRunNumber(279676);
    plugin->AddRunNumber(279642); plugin->AddRunNumber(279641); plugin->AddRunNumber(279600); plugin->AddRunNumber(279598); plugin->AddRunNumber(279597); plugin->AddRunNumber(279583); plugin->AddRunNumber(279565);
    plugin->AddRunNumber(279564); plugin->AddRunNumber(279563); plugin->AddRunNumber(279562); plugin->AddRunNumber(279561); plugin->AddRunNumber(279560); plugin->AddRunNumber(279559); plugin->AddRunNumber(279488);
    plugin->AddRunNumber(279487); plugin->AddRunNumber(279483); plugin->AddRunNumber(279441); plugin->AddRunNumber(279439); plugin->AddRunNumber(279435); plugin->AddRunNumber(279410); plugin->AddRunNumber(279391);
    plugin->AddRunNumber(279355); plugin->AddRunNumber(279354); plugin->AddRunNumber(279349); plugin->AddRunNumber(279348); plugin->AddRunNumber(279344); plugin->AddRunNumber(279342); plugin->AddRunNumber(279312);
    plugin->AddRunNumber(279310); plugin->AddRunNumber(279309); plugin->AddRunNumber(279274); plugin->AddRunNumber(279273); plugin->AddRunNumber(279270); plugin->AddRunNumber(279268); plugin->AddRunNumber(279267);
    plugin->AddRunNumber(279265); plugin->AddRunNumber(279264); plugin->AddRunNumber(279242); plugin->AddRunNumber(279238); plugin->AddRunNumber(279235); plugin->AddRunNumber(279234); plugin->AddRunNumber(279208);
    plugin->AddRunNumber(279207); plugin->AddRunNumber(279201); plugin->AddRunNumber(279199); plugin->AddRunNumber(279157); plugin->AddRunNumber(279155); plugin->AddRunNumber(279130); plugin->AddRunNumber(279125);
    plugin->AddRunNumber(279123); plugin->AddRunNumber(279122); plugin->AddRunNumber(279117); plugin->AddRunNumber(279106); plugin->AddRunNumber(279075); plugin->AddRunNumber(279074); plugin->AddRunNumber(279073);
    plugin->AddRunNumber(279068); plugin->AddRunNumber(279044); plugin->AddRunNumber(279043); plugin->AddRunNumber(279041); plugin->AddRunNumber(279038); plugin->AddRunNumber(279037); plugin->AddRunNumber(279036);
    plugin->AddRunNumber(279008); plugin->AddRunNumber(279007); plugin->AddRunNumber(279005); plugin->AddRunNumber(278999); plugin->AddRunNumber(278964); plugin->AddRunNumber(278963); plugin->AddRunNumber(278959);
    plugin->AddRunNumber(278941); plugin->AddRunNumber(278939); plugin->AddRunNumber(278936); plugin->AddRunNumber(278915); plugin->AddRunNumber(278914);
  }
  if(runPeriod.Contains("LHC17o")){
    plugin->AddRunNumber(281961); plugin->AddRunNumber(281956); plugin->AddRunNumber(281953); plugin->AddRunNumber(281946); plugin->AddRunNumber(281940); plugin->AddRunNumber(281939); plugin->AddRunNumber(281931);
    plugin->AddRunNumber(281928); plugin->AddRunNumber(281918); plugin->AddRunNumber(281916); plugin->AddRunNumber(281915); plugin->AddRunNumber(281894); plugin->AddRunNumber(281893); plugin->AddRunNumber(281892);
    plugin->AddRunNumber(281755); plugin->AddRunNumber(281754); plugin->AddRunNumber(281753); plugin->AddRunNumber(281751); plugin->AddRunNumber(281750); plugin->AddRunNumber(281741); plugin->AddRunNumber(281713);
    plugin->AddRunNumber(281709); plugin->AddRunNumber(281707); plugin->AddRunNumber(281706); plugin->AddRunNumber(281705); plugin->AddRunNumber(281672); plugin->AddRunNumber(281667); plugin->AddRunNumber(281664);
    plugin->AddRunNumber(281658); plugin->AddRunNumber(281655); plugin->AddRunNumber(281654); plugin->AddRunNumber(281651); plugin->AddRunNumber(281645); plugin->AddRunNumber(281642); plugin->AddRunNumber(281640);
    plugin->AddRunNumber(281635); plugin->AddRunNumber(281634); plugin->AddRunNumber(281633); plugin->AddRunNumber(281592); plugin->AddRunNumber(281583); plugin->AddRunNumber(281581); plugin->AddRunNumber(281580);
    plugin->AddRunNumber(281574); plugin->AddRunNumber(281569); plugin->AddRunNumber(281568); plugin->AddRunNumber(281563); plugin->AddRunNumber(281562); plugin->AddRunNumber(281557); plugin->AddRunNumber(281511);
    plugin->AddRunNumber(281509); plugin->AddRunNumber(281477); plugin->AddRunNumber(281475); plugin->AddRunNumber(281450); plugin->AddRunNumber(281449); plugin->AddRunNumber(281446); plugin->AddRunNumber(281444);
    plugin->AddRunNumber(281441); plugin->AddRunNumber(281415); plugin->AddRunNumber(281321); plugin->AddRunNumber(281301); plugin->AddRunNumber(281277); plugin->AddRunNumber(281275); plugin->AddRunNumber(281244);
    plugin->AddRunNumber(281243); plugin->AddRunNumber(281242); plugin->AddRunNumber(281241); plugin->AddRunNumber(281240); plugin->AddRunNumber(281213); plugin->AddRunNumber(281212); plugin->AddRunNumber(281191);
    plugin->AddRunNumber(281190); plugin->AddRunNumber(281181); plugin->AddRunNumber(281180); plugin->AddRunNumber(281179); plugin->AddRunNumber(281081); plugin->AddRunNumber(281080); plugin->AddRunNumber(281079);
    plugin->AddRunNumber(281062); plugin->AddRunNumber(281061); plugin->AddRunNumber(281060); plugin->AddRunNumber(281036); plugin->AddRunNumber(281035); plugin->AddRunNumber(281033); plugin->AddRunNumber(281032);
    plugin->AddRunNumber(280998); plugin->AddRunNumber(280997); plugin->AddRunNumber(280996); plugin->AddRunNumber(280994); plugin->AddRunNumber(280990); plugin->AddRunNumber(280947); plugin->AddRunNumber(280943);
    plugin->AddRunNumber(280940); plugin->AddRunNumber(280936); plugin->AddRunNumber(280897); plugin->AddRunNumber(280890); plugin->AddRunNumber(280881); plugin->AddRunNumber(280880); plugin->AddRunNumber(280856);
    plugin->AddRunNumber(280848); plugin->AddRunNumber(280847); plugin->AddRunNumber(280845); plugin->AddRunNumber(280844); plugin->AddRunNumber(280842); plugin->AddRunNumber(280793); plugin->AddRunNumber(280792);
    plugin->AddRunNumber(280786); plugin->AddRunNumber(280768); plugin->AddRunNumber(280767); plugin->AddRunNumber(280766); plugin->AddRunNumber(280765); plugin->AddRunNumber(280764); plugin->AddRunNumber(280763);
    plugin->AddRunNumber(280761); plugin->AddRunNumber(280756); plugin->AddRunNumber(280755); plugin->AddRunNumber(280754); plugin->AddRunNumber(280753); plugin->AddRunNumber(280706); plugin->AddRunNumber(280705);
    plugin->AddRunNumber(280681); plugin->AddRunNumber(280679); plugin->AddRunNumber(280676); plugin->AddRunNumber(280671); plugin->AddRunNumber(280650); plugin->AddRunNumber(280648); plugin->AddRunNumber(280647);
    plugin->AddRunNumber(280645); plugin->AddRunNumber(280639); plugin->AddRunNumber(280637); plugin->AddRunNumber(280634); plugin->AddRunNumber(280613); plugin->AddRunNumber(280583); plugin->AddRunNumber(280581);
    plugin->AddRunNumber(280576); plugin->AddRunNumber(280575); plugin->AddRunNumber(280574); plugin->AddRunNumber(280551); plugin->AddRunNumber(280550); plugin->AddRunNumber(280547); plugin->AddRunNumber(280546);
    plugin->AddRunNumber(280519); plugin->AddRunNumber(280518); plugin->AddRunNumber(280448); plugin->AddRunNumber(280447); plugin->AddRunNumber(280446); plugin->AddRunNumber(280445); plugin->AddRunNumber(280443);
    plugin->AddRunNumber(280419); plugin->AddRunNumber(280418); plugin->AddRunNumber(280415); plugin->AddRunNumber(280413); plugin->AddRunNumber(280412); plugin->AddRunNumber(280406); plugin->AddRunNumber(280405);
    plugin->AddRunNumber(280403); plugin->AddRunNumber(280375); plugin->AddRunNumber(280374); plugin->AddRunNumber(280352); plugin->AddRunNumber(280351); plugin->AddRunNumber(280350); plugin->AddRunNumber(280349);
    plugin->AddRunNumber(280348); plugin->AddRunNumber(280312); plugin->AddRunNumber(280310); plugin->AddRunNumber(280290); plugin->AddRunNumber(280286); plugin->AddRunNumber(280285); plugin->AddRunNumber(280284);
    plugin->AddRunNumber(280283); plugin->AddRunNumber(280282);
  }
  if(runPeriod.Contains("LHC17r")){
    plugin->AddRunNumber(282704); plugin->AddRunNumber(282703); plugin->AddRunNumber(282702); plugin->AddRunNumber(282700); plugin->AddRunNumber(282677); plugin->AddRunNumber(282676); plugin->AddRunNumber(282673);
    plugin->AddRunNumber(282671); plugin->AddRunNumber(282670); plugin->AddRunNumber(282668); plugin->AddRunNumber(282667); plugin->AddRunNumber(282666); plugin->AddRunNumber(282653); plugin->AddRunNumber(282651);
    plugin->AddRunNumber(282629); plugin->AddRunNumber(282622); plugin->AddRunNumber(282620); plugin->AddRunNumber(282618); plugin->AddRunNumber(282615); plugin->AddRunNumber(282609); plugin->AddRunNumber(282608);
    plugin->AddRunNumber(282607); plugin->AddRunNumber(282606); plugin->AddRunNumber(282580); plugin->AddRunNumber(282579); plugin->AddRunNumber(282575); plugin->AddRunNumber(282573); plugin->AddRunNumber(282546);
    plugin->AddRunNumber(282545); plugin->AddRunNumber(282544); plugin->AddRunNumber(282528); plugin->AddRunNumber(282504);
  }

  if(runPeriod.Contains("LHC18c")){
    plugin->AddRunNumber(285958); plugin->AddRunNumber(285957); plugin->AddRunNumber(285946); plugin->AddRunNumber(285917); plugin->AddRunNumber(285893); plugin->AddRunNumber(285892); plugin->AddRunNumber(285869);
    plugin->AddRunNumber(285851); plugin->AddRunNumber(285830); plugin->AddRunNumber(285812); plugin->AddRunNumber(285811); plugin->AddRunNumber(285810); plugin->AddRunNumber(285806); plugin->AddRunNumber(285805);
    plugin->AddRunNumber(285804); plugin->AddRunNumber(285781); plugin->AddRunNumber(285778); plugin->AddRunNumber(285777); plugin->AddRunNumber(285756); plugin->AddRunNumber(285755); plugin->AddRunNumber(285754);
    plugin->AddRunNumber(285753); plugin->AddRunNumber(285752); plugin->AddRunNumber(285751); plugin->AddRunNumber(285722); plugin->AddRunNumber(285698); plugin->AddRunNumber(285697); plugin->AddRunNumber(285664);
    plugin->AddRunNumber(285663); plugin->AddRunNumber(285662); plugin->AddRunNumber(285659); plugin->AddRunNumber(285643); plugin->AddRunNumber(285642); plugin->AddRunNumber(285641); plugin->AddRunNumber(285640);
    plugin->AddRunNumber(285639); plugin->AddRunNumber(285603); plugin->AddRunNumber(285602); plugin->AddRunNumber(285601); plugin->AddRunNumber(285599); plugin->AddRunNumber(285578); plugin->AddRunNumber(285577);
    plugin->AddRunNumber(285576); plugin->AddRunNumber(285575); plugin->AddRunNumber(285557); plugin->AddRunNumber(285515); plugin->AddRunNumber(285497); plugin->AddRunNumber(285496);
  }
  if(runPeriod.Contains("LHC18d")){
    plugin->AddRunNumber(286350); plugin->AddRunNumber(286349); plugin->AddRunNumber(286348); plugin->AddRunNumber(286345); plugin->AddRunNumber(286340); plugin->AddRunNumber(286337); plugin->AddRunNumber(286336); plugin->AddRunNumber(286314); plugin->AddRunNumber(286313); plugin->AddRunNumber(286312); plugin->AddRunNumber(286311); plugin->AddRunNumber(286310); plugin->AddRunNumber(286309); plugin->AddRunNumber(286308); plugin->AddRunNumber(286289); plugin->AddRunNumber(286288); plugin->AddRunNumber(286287); plugin->AddRunNumber(286284); plugin->AddRunNumber(286282); plugin->AddRunNumber(286261); plugin->AddRunNumber(286258); plugin->AddRunNumber(286257); plugin->AddRunNumber(286254); plugin->AddRunNumber(286230); plugin->AddRunNumber(286229); plugin->AddRunNumber(286203); plugin->AddRunNumber(286202); plugin->AddRunNumber(286201); plugin->AddRunNumber(286199); plugin->AddRunNumber(286198); plugin->AddRunNumber(286159); plugin->AddRunNumber(286130); plugin->AddRunNumber(286129); plugin->AddRunNumber(286127); plugin->AddRunNumber(286124); plugin->AddRunNumber(286064); plugin->AddRunNumber(286028); plugin->AddRunNumber(286027); plugin->AddRunNumber(286026); plugin->AddRunNumber(286025); plugin->AddRunNumber(286018); plugin->AddRunNumber(286014); plugin->AddRunNumber(285980); plugin->AddRunNumber(285979); plugin->AddRunNumber(285978);
  }
  if(runPeriod.Contains("LHC18e")){
    plugin->AddRunNumber(286937); plugin->AddRunNumber(286936); plugin->AddRunNumber(286933); plugin->AddRunNumber(286932); plugin->AddRunNumber(286931); plugin->AddRunNumber(286930); plugin->AddRunNumber(286911); plugin->AddRunNumber(286910); plugin->AddRunNumber(286908); plugin->AddRunNumber(286907); plugin->AddRunNumber(286877); plugin->AddRunNumber(286876); plugin->AddRunNumber(286874); plugin->AddRunNumber(286852); plugin->AddRunNumber(286850); plugin->AddRunNumber(286848); plugin->AddRunNumber(286846); plugin->AddRunNumber(286810); plugin->AddRunNumber(286809); plugin->AddRunNumber(286805); plugin->AddRunNumber(286801); plugin->AddRunNumber(286799); plugin->AddRunNumber(286731); plugin->AddRunNumber(286695); plugin->AddRunNumber(286661); plugin->AddRunNumber(286653); plugin->AddRunNumber(286633); plugin->AddRunNumber(286594); plugin->AddRunNumber(286592); plugin->AddRunNumber(286591); plugin->AddRunNumber(286569); plugin->AddRunNumber(286568); plugin->AddRunNumber(286567); plugin->AddRunNumber(286566); plugin->AddRunNumber(286509); plugin->AddRunNumber(286508); plugin->AddRunNumber(286502); plugin->AddRunNumber(286501); plugin->AddRunNumber(286455); plugin->AddRunNumber(286454); plugin->AddRunNumber(286428); plugin->AddRunNumber(286427); plugin->AddRunNumber(286426); plugin->AddRunNumber(286380);
  }
  if(runPeriod.Contains("LHC18f")){
    plugin->AddRunNumber(287977); plugin->AddRunNumber(287975); plugin->AddRunNumber(287941); plugin->AddRunNumber(287923); plugin->AddRunNumber(287784); plugin->AddRunNumber(287783); plugin->AddRunNumber(287658); plugin->AddRunNumber(287657); plugin->AddRunNumber(287656); plugin->AddRunNumber(287654); plugin->AddRunNumber(287578); plugin->AddRunNumber(287576); plugin->AddRunNumber(287575); plugin->AddRunNumber(287573); plugin->AddRunNumber(287524); plugin->AddRunNumber(287521); plugin->AddRunNumber(287520); plugin->AddRunNumber(287518); plugin->AddRunNumber(287517); plugin->AddRunNumber(287516); plugin->AddRunNumber(287513); plugin->AddRunNumber(287486); plugin->AddRunNumber(287484); plugin->AddRunNumber(287481); plugin->AddRunNumber(287480); plugin->AddRunNumber(287451); plugin->AddRunNumber(287413); plugin->AddRunNumber(287389); plugin->AddRunNumber(287388); plugin->AddRunNumber(287387); plugin->AddRunNumber(287385); plugin->AddRunNumber(287381); plugin->AddRunNumber(287380); plugin->AddRunNumber(287360); plugin->AddRunNumber(287358); plugin->AddRunNumber(287356); plugin->AddRunNumber(287355); plugin->AddRunNumber(287353); plugin->AddRunNumber(287349); plugin->AddRunNumber(287347); plugin->AddRunNumber(287346); plugin->AddRunNumber(287344); plugin->AddRunNumber(287343); plugin->AddRunNumber(287325); plugin->AddRunNumber(287324); plugin->AddRunNumber(287323); plugin->AddRunNumber(287283); plugin->AddRunNumber(287254); plugin->AddRunNumber(287251); plugin->AddRunNumber(287250); plugin->AddRunNumber(287249); plugin->AddRunNumber(287248); plugin->AddRunNumber(287209); plugin->AddRunNumber(287208); plugin->AddRunNumber(287204); plugin->AddRunNumber(287203); plugin->AddRunNumber(287202); plugin->AddRunNumber(287201); plugin->AddRunNumber(287155); plugin->AddRunNumber(287137); plugin->AddRunNumber(287077); plugin->AddRunNumber(287072); plugin->AddRunNumber(287071); plugin->AddRunNumber(287066); plugin->AddRunNumber(287064); plugin->AddRunNumber(287063); plugin->AddRunNumber(287021); plugin->AddRunNumber(287000);
  }
  if(runPeriod.Contains("LHC18l")){
    plugin->AddRunNumber(289971); plugin->AddRunNumber(289966); plugin->AddRunNumber(289943); plugin->AddRunNumber(289941); plugin->AddRunNumber(289940); plugin->AddRunNumber(289935); plugin->AddRunNumber(289931); plugin->AddRunNumber(289928); plugin->AddRunNumber(289888); plugin->AddRunNumber(289884); plugin->AddRunNumber(289880); plugin->AddRunNumber(289857); plugin->AddRunNumber(289856); plugin->AddRunNumber(289855); plugin->AddRunNumber(289852); plugin->AddRunNumber(289849); plugin->AddRunNumber(289830); plugin->AddRunNumber(289816); plugin->AddRunNumber(289815); plugin->AddRunNumber(289814); plugin->AddRunNumber(289811); plugin->AddRunNumber(289808); plugin->AddRunNumber(289775); plugin->AddRunNumber(289757); plugin->AddRunNumber(289731); plugin->AddRunNumber(289729); plugin->AddRunNumber(289724); plugin->AddRunNumber(289723); plugin->AddRunNumber(289721); plugin->AddRunNumber(289666); plugin->AddRunNumber(289664); plugin->AddRunNumber(289660); plugin->AddRunNumber(289659); plugin->AddRunNumber(289658); plugin->AddRunNumber(289657); plugin->AddRunNumber(289654); plugin->AddRunNumber(289632); plugin->AddRunNumber(289626); plugin->AddRunNumber(289625); plugin->AddRunNumber(289582); plugin->AddRunNumber(289581); plugin->AddRunNumber(289579); plugin->AddRunNumber(289577); plugin->AddRunNumber(289576); plugin->AddRunNumber(289574); plugin->AddRunNumber(289547); plugin->AddRunNumber(289494); plugin->AddRunNumber(289493); plugin->AddRunNumber(289468); plugin->AddRunNumber(289466); plugin->AddRunNumber(289465); plugin->AddRunNumber(289463); plugin->AddRunNumber(289462); plugin->AddRunNumber(289444); plugin->AddRunNumber(289426); plugin->AddRunNumber(289373); plugin->AddRunNumber(289370); plugin->AddRunNumber(289369); plugin->AddRunNumber(289368); plugin->AddRunNumber(289367); plugin->AddRunNumber(289366); plugin->AddRunNumber(289365); plugin->AddRunNumber(289363); plugin->AddRunNumber(289356); plugin->AddRunNumber(289355); plugin->AddRunNumber(289354); plugin->AddRunNumber(289353); plugin->AddRunNumber(289309); plugin->AddRunNumber(289308); plugin->AddRunNumber(289306); plugin->AddRunNumber(289303); plugin->AddRunNumber(289300); plugin->AddRunNumber(289280); plugin->AddRunNumber(289278); plugin->AddRunNumber(289277); plugin->AddRunNumber(289276); plugin->AddRunNumber(289275); plugin->AddRunNumber(289254); plugin->AddRunNumber(289253); plugin->AddRunNumber(289249); plugin->AddRunNumber(289247); plugin->AddRunNumber(289243); plugin->AddRunNumber(289242); plugin->AddRunNumber(289241); plugin->AddRunNumber(289240);
  }
  if(runPeriod.Contains("LHC18m")){
    plugin->AddRunNumber(292397); plugin->AddRunNumber(292298); plugin->AddRunNumber(292274); plugin->AddRunNumber(292273); plugin->AddRunNumber(292270); plugin->AddRunNumber(292269); plugin->AddRunNumber(292265); plugin->AddRunNumber(292242); plugin->AddRunNumber(292241); plugin->AddRunNumber(292240); plugin->AddRunNumber(292192); plugin->AddRunNumber(292168); plugin->AddRunNumber(292167); plugin->AddRunNumber(292166); plugin->AddRunNumber(292164); plugin->AddRunNumber(292163); plugin->AddRunNumber(292162); plugin->AddRunNumber(292161); plugin->AddRunNumber(292160); plugin->AddRunNumber(292140); plugin->AddRunNumber(292115); plugin->AddRunNumber(292114); plugin->AddRunNumber(292109); plugin->AddRunNumber(292108); plugin->AddRunNumber(292107); plugin->AddRunNumber(292106); plugin->AddRunNumber(292081); plugin->AddRunNumber(292080); plugin->AddRunNumber(292077); plugin->AddRunNumber(292075); plugin->AddRunNumber(292062); plugin->AddRunNumber(292061); plugin->AddRunNumber(292060); plugin->AddRunNumber(292040); plugin->AddRunNumber(292012); plugin->AddRunNumber(291982); plugin->AddRunNumber(291976); plugin->AddRunNumber(291953); plugin->AddRunNumber(291948); plugin->AddRunNumber(291945); plugin->AddRunNumber(291944); plugin->AddRunNumber(291943); plugin->AddRunNumber(291942); plugin->AddRunNumber(291803); plugin->AddRunNumber(291796); plugin->AddRunNumber(291795); plugin->AddRunNumber(291769); plugin->AddRunNumber(291760); plugin->AddRunNumber(291756); plugin->AddRunNumber(291755); plugin->AddRunNumber(291729); plugin->AddRunNumber(291706); plugin->AddRunNumber(291698); plugin->AddRunNumber(291697); plugin->AddRunNumber(291694); plugin->AddRunNumber(291692); plugin->AddRunNumber(291690); plugin->AddRunNumber(291665); plugin->AddRunNumber(291661); plugin->AddRunNumber(291657); plugin->AddRunNumber(291626); plugin->AddRunNumber(291625); plugin->AddRunNumber(291624); plugin->AddRunNumber(291622); plugin->AddRunNumber(291618); plugin->AddRunNumber(291615); plugin->AddRunNumber(291614); plugin->AddRunNumber(291590); plugin->AddRunNumber(291485); plugin->AddRunNumber(291484); plugin->AddRunNumber(291482); plugin->AddRunNumber(291481); plugin->AddRunNumber(291457); plugin->AddRunNumber(291456); plugin->AddRunNumber(291453); plugin->AddRunNumber(291451); plugin->AddRunNumber(291447); plugin->AddRunNumber(291446); plugin->AddRunNumber(291420); plugin->AddRunNumber(291419); plugin->AddRunNumber(291417); plugin->AddRunNumber(291416); plugin->AddRunNumber(291402); plugin->AddRunNumber(291400); plugin->AddRunNumber(291399); plugin->AddRunNumber(291397); plugin->AddRunNumber(291375); plugin->AddRunNumber(291373); plugin->AddRunNumber(291363); plugin->AddRunNumber(291362); plugin->AddRunNumber(291361); plugin->AddRunNumber(291360); plugin->AddRunNumber(291286); plugin->AddRunNumber(291285); plugin->AddRunNumber(291284); plugin->AddRunNumber(291283); plugin->AddRunNumber(291282); plugin->AddRunNumber(291265); plugin->AddRunNumber(291263); plugin->AddRunNumber(291110); plugin->AddRunNumber(291100); plugin->AddRunNumber(291066); plugin->AddRunNumber(291065); plugin->AddRunNumber(291041); plugin->AddRunNumber(291037); plugin->AddRunNumber(291035); plugin->AddRunNumber(291006); plugin->AddRunNumber(291005); plugin->AddRunNumber(291004); plugin->AddRunNumber(291003); plugin->AddRunNumber(291002); plugin->AddRunNumber(290980); plugin->AddRunNumber(290979); plugin->AddRunNumber(290976); plugin->AddRunNumber(290975); plugin->AddRunNumber(290948); plugin->AddRunNumber(290944); plugin->AddRunNumber(290943); plugin->AddRunNumber(290935); plugin->AddRunNumber(290932); plugin->AddRunNumber(290895); plugin->AddRunNumber(290894); plugin->AddRunNumber(290892); plugin->AddRunNumber(290862); plugin->AddRunNumber(290860); plugin->AddRunNumber(290853); plugin->AddRunNumber(290848); plugin->AddRunNumber(290790); plugin->AddRunNumber(290787); plugin->AddRunNumber(290776); plugin->AddRunNumber(290774); plugin->AddRunNumber(290769); plugin->AddRunNumber(290766); plugin->AddRunNumber(290764); plugin->AddRunNumber(290742); plugin->AddRunNumber(290721); plugin->AddRunNumber(290699); plugin->AddRunNumber(290696); plugin->AddRunNumber(290692); plugin->AddRunNumber(290687); plugin->AddRunNumber(290665); plugin->AddRunNumber(290660); plugin->AddRunNumber(290658); plugin->AddRunNumber(290645); plugin->AddRunNumber(290632); plugin->AddRunNumber(290627); plugin->AddRunNumber(290615); plugin->AddRunNumber(290614); plugin->AddRunNumber(290613); plugin->AddRunNumber(290612); plugin->AddRunNumber(290590); plugin->AddRunNumber(290553); plugin->AddRunNumber(290550); plugin->AddRunNumber(290549); plugin->AddRunNumber(290544); plugin->AddRunNumber(290540); plugin->AddRunNumber(290539); plugin->AddRunNumber(290538); plugin->AddRunNumber(290501); plugin->AddRunNumber(290499); plugin->AddRunNumber(290469); plugin->AddRunNumber(290467); plugin->AddRunNumber(290459); plugin->AddRunNumber(290458); plugin->AddRunNumber(290456); plugin->AddRunNumber(290428); plugin->AddRunNumber(290427); plugin->AddRunNumber(290425); plugin->AddRunNumber(290423); plugin->AddRunNumber(290421); plugin->AddRunNumber(290420); plugin->AddRunNumber(290418); plugin->AddRunNumber(290411); plugin->AddRunNumber(290404); plugin->AddRunNumber(290401); plugin->AddRunNumber(290375); plugin->AddRunNumber(290374); plugin->AddRunNumber(290350); plugin->AddRunNumber(290327); plugin->AddRunNumber(290324); plugin->AddRunNumber(290323); plugin->AddRunNumber(290300); plugin->AddRunNumber(290297); plugin->AddRunNumber(290293); plugin->AddRunNumber(290254); plugin->AddRunNumber(290223); plugin->AddRunNumber(290222);
  }
  if(runPeriod.Contains("LHC18o")){
    plugin->AddRunNumber(293898); plugin->AddRunNumber(293896); plugin->AddRunNumber(293893); plugin->AddRunNumber(293891); plugin->AddRunNumber(293886); plugin->AddRunNumber(293856); plugin->AddRunNumber(293831); plugin->AddRunNumber(293830); plugin->AddRunNumber(293829); plugin->AddRunNumber(293809); plugin->AddRunNumber(293807); plugin->AddRunNumber(293806); plugin->AddRunNumber(293805); plugin->AddRunNumber(293802); plugin->AddRunNumber(293799); plugin->AddRunNumber(293776); plugin->AddRunNumber(293774); plugin->AddRunNumber(293773); plugin->AddRunNumber(293741); plugin->AddRunNumber(293740); plugin->AddRunNumber(293698); plugin->AddRunNumber(293696); plugin->AddRunNumber(293695); plugin->AddRunNumber(293692); plugin->AddRunNumber(293691); plugin->AddRunNumber(293588); plugin->AddRunNumber(293587); plugin->AddRunNumber(293497); plugin->AddRunNumber(293496); plugin->AddRunNumber(293494); plugin->AddRunNumber(293475); plugin->AddRunNumber(293474); plugin->AddRunNumber(293424); plugin->AddRunNumber(293413); plugin->AddRunNumber(293392); plugin->AddRunNumber(293391); plugin->AddRunNumber(293388); plugin->AddRunNumber(293386); plugin->AddRunNumber(293368);
  }
  if(runPeriod.Contains("LHC18p")){
    plugin->AddRunNumber(294925); plugin->AddRunNumber(294916); plugin->AddRunNumber(294884); plugin->AddRunNumber(294883); plugin->AddRunNumber(294880); plugin->AddRunNumber(294877); plugin->AddRunNumber(294875); plugin->AddRunNumber(294852); plugin->AddRunNumber(294818); plugin->AddRunNumber(294817); plugin->AddRunNumber(294816); plugin->AddRunNumber(294815); plugin->AddRunNumber(294813); plugin->AddRunNumber(294809); plugin->AddRunNumber(294775); plugin->AddRunNumber(294774); plugin->AddRunNumber(294772); plugin->AddRunNumber(294769); plugin->AddRunNumber(294749); plugin->AddRunNumber(294747); plugin->AddRunNumber(294743); plugin->AddRunNumber(294742); plugin->AddRunNumber(294741); plugin->AddRunNumber(294722); plugin->AddRunNumber(294721); plugin->AddRunNumber(294718); plugin->AddRunNumber(294716); plugin->AddRunNumber(294715); plugin->AddRunNumber(294710); plugin->AddRunNumber(294703); plugin->AddRunNumber(294653); plugin->AddRunNumber(294636); plugin->AddRunNumber(294634); plugin->AddRunNumber(294633); plugin->AddRunNumber(294632); plugin->AddRunNumber(294593); plugin->AddRunNumber(294591); plugin->AddRunNumber(294590); plugin->AddRunNumber(294588); plugin->AddRunNumber(294587); plugin->AddRunNumber(294586); plugin->AddRunNumber(294563); plugin->AddRunNumber(294558); plugin->AddRunNumber(294556); plugin->AddRunNumber(294553); plugin->AddRunNumber(294531); plugin->AddRunNumber(294530); plugin->AddRunNumber(294529); plugin->AddRunNumber(294527); plugin->AddRunNumber(294526); plugin->AddRunNumber(294525); plugin->AddRunNumber(294524); plugin->AddRunNumber(294503); plugin->AddRunNumber(294502); plugin->AddRunNumber(294310); plugin->AddRunNumber(294308); plugin->AddRunNumber(294307); plugin->AddRunNumber(294305); plugin->AddRunNumber(294242); plugin->AddRunNumber(294241); plugin->AddRunNumber(294212); plugin->AddRunNumber(294210); plugin->AddRunNumber(294208); plugin->AddRunNumber(294205); plugin->AddRunNumber(294201); plugin->AddRunNumber(294200); plugin->AddRunNumber(294199); plugin->AddRunNumber(294156); plugin->AddRunNumber(294155); plugin->AddRunNumber(294154); plugin->AddRunNumber(294152); plugin->AddRunNumber(294131); plugin->AddRunNumber(294128); plugin->AddRunNumber(294013); plugin->AddRunNumber(294012); plugin->AddRunNumber(294011); plugin->AddRunNumber(294010); plugin->AddRunNumber(294009);
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
  else               plugin->SetNtestFiles(10);
  
  
  
  plugin->SetOutputToRunNo();
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName("analysis_syano_MC.jdl");
  plugin->SetPrice(1);      
  plugin->SetSplitMode("se");
  
  return plugin;
}
