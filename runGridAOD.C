void runGridAOD(TString period="LHC17q_woSDD",TString mode = "full",Bool_t isJDL=true,TString type="data",Bool_t isMC=false,Bool_t isAOD=true){
  
  Bool_t isMC = false;
  if(type != "data")
    isMC = true;
  
  LoadLib();

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(period,mode,isJDL,isAOD,type);  
  if (!alienHandler) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Task_syano");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
    
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  UInt_t offlineTriggerMask    = AliVEvent::kINT7;
  
  gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,true);
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* taskPID = AddTaskPIDResponse(isMC,false,true,1); 
    
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * MultiTask = AddTaskMultSelection(kFALSE);
  MultiTask->SetUseDefaultCalib(kTRUE);

  if(isMC) MultiTask->SetUseDefaultMCCalib(kTRUE);
  else     MultiTask->SetUseDefaultCalib(kTRUE);
  
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPairUtils.cxx++g");
  
  gROOT->LoadMacro("AliAnalysisTaskAODTrackPairEventCounter.cxx++g");
  gROOT->LoadMacro("AddTaskAODTrackPairEventCounter.C");
  AliAnalysisTaskAODTrackPairEventCounter* taskEvent = AddTaskAODTrackPairEventCounter("EventCounter",period,offlineTriggerMask);
  taskEvent->SetMC(isMC);

  gROOT->LoadMacro("AliAnalysisTaskAODTrackPair.cxx++g");
  gROOT->LoadMacro("AddTaskAODTrackPair.C");
  AliAnalysisTaskAODTrackPair* taskPair = AddTaskAODTrackPair("TrackPair",period,offlineTriggerMask);
  taskPair->SetMC(isMC);
  
  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
};

// Loading libraries
void LoadLib() {
  printf("Loading libraries...");
  // loading libraries                                                                                                                                                                                  
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY");
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libOADB");
  gSystem->Load("libCORRFW");
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  printf("done\n");


}

