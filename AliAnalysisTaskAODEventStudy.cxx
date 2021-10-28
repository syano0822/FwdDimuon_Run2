#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TRefArray.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolMuon.h"
#include "AliEventPoolManager.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"


#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODDimuon.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

#include "AliEventplane.h"
#include "AliMultSelection.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMuonTrackCuts.h"
#include "AliMultSelectionTask.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskAODEventStudy.h"
#include "THnSparse.h"

#include "AliTHn.h"

#include "iostream"

#include "AliAnalysisTaskAODTrackPairUtils.h"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODEventStudy)       


//________________________________________________________________________
AliAnalysisTaskAODEventStudy::AliAnalysisTaskAODEventStudy() : AliAnalysisTaskSE(),
  fEvent(NULL),
  fUtils(NULL),

  fRunNumber(999999),
  
  fNTrackMatchAll(0),
  fNTrackMatchApt(0),
  fNTrackMatchLpt(0),
  fNTrackMatchHpt(0),
  
  fIsMC(false),
  fIsCINT7(false),
  fIsCMSL7(false),
  fIsCMSH7(false),
  fIsCMUL7(false),
  fIsCMLL7(false),  
  fIs0MSL(false),
  fIs0MSH(false),
  fIs0MUL(false),
  fIs0MLL(false),
  
  trig4NormalizationMap(),

  fHistEvent(),
  fHistEvent4RunQA(),
  
  fHistNTrackMatchAll(),
  fHistNTrackMatchApt(),
  fHistNTrackMatchLpt(),
  fHistNTrackMatchHpt(),
  fHistMeanPtMatchAll(),
  fHistMeanPtMatchApt(),
  fHistMeanPtMatchLpt(),
  fHistMeanPtMatchHpt(),
  
  fHistCentrality(),
  fHistCentVtxZ(),  
  fHistEtaPt(),
  
  fHistSPDTrkEta05Centrality(),
  fHistSPDTrkEta10Centrality(),
  fHistSPDTrkEta15Centrality(),
  fHistSPDTrkEta20Centrality(),
  fHistSPDTrkEtaAllCentrality(),
  fHistSPDTrkSPDClust(),

  fHistChV0ACentrality(),
  fHistChV0CCentrality(),
  fHistChV0MCentrality(),

  fHistCorrV0MCentrality(),
  fHistCorrV0ACentrality(),
  fHistCorrV0CCentrality(),
  fHistCorrSPDTrkCentrality(),

//fHistCorrChV0AChV0C(),
//fHistCorrCentSPDTrkV0M(),
//fHistCorrCentSPDTrkV0A(),
//fHistCorrCentSPDTrkV0C(),
//fHistCorrCentV0AV0C(),

  fHistPtMatchAll(),
  fHistPtMatchApt(),
  fHistPtMatchLpt(),
  fHistPtMatchHpt(),
  
  fHistRabs(NULL),
  fHistTrackChi2(),
  fHistTrigChi2(),
  fHistChVtxEta(NULL)
{
  
}

//________________________________________________________________________
AliAnalysisTaskAODEventStudy::AliAnalysisTaskAODEventStudy(const char* name) : 
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fUtils(NULL),
  fIsMC(false),
  
  fRunNumber(999999),
  
  fNTrackMatchAll(0),
  fNTrackMatchApt(0),
  fNTrackMatchLpt(0),
  fNTrackMatchHpt(0),

  fIsCINT7(false),
  fIsCMSL7(false),
  fIsCMSH7(false),
  fIsCMUL7(false),
  fIsCMLL7(false),  
  fIs0MSL(false),
  fIs0MSH(false),
  fIs0MUL(false),
  fIs0MLL(false),

  trig4NormalizationMap(),

  fHistEvent(),
  fHistEvent4RunQA(),

  fHistNTrackMatchAll(),
  fHistNTrackMatchApt(),
  fHistNTrackMatchLpt(),
  fHistNTrackMatchHpt(),
  fHistMeanPtMatchAll(),
  fHistMeanPtMatchApt(),
  fHistMeanPtMatchLpt(),
  fHistMeanPtMatchHpt(),
  
  fHistCentrality(),
  fHistCentVtxZ(),  
  fHistEtaPt(),

  fHistSPDTrkEta05Centrality(),
  fHistSPDTrkEta10Centrality(),
  fHistSPDTrkEta15Centrality(),
  fHistSPDTrkEta20Centrality(),
  fHistSPDTrkEtaAllCentrality(),
  fHistSPDTrkSPDClust(),

  fHistChV0ACentrality(),
  fHistChV0CCentrality(),
  fHistChV0MCentrality(),

  fHistCorrV0MCentrality(),
  fHistCorrV0ACentrality(),
  fHistCorrV0CCentrality(),
  fHistCorrSPDTrkCentrality(),

//fHistCorrChV0AChV0C(),
//fHistCorrCentSPDTrkV0M(),
//fHistCorrCentSPDTrkV0A(),
//fHistCorrCentSPDTrkV0C(),
//fHistCorrCentV0AV0C(),

  fHistPtMatchAll(),
  fHistPtMatchApt(),
  fHistPtMatchLpt(),
  fHistPtMatchHpt(),
  
  fHistRabs(NULL),
  fHistTrackChi2(),
  fHistTrigChi2(),
  fHistChVtxEta(NULL)
{  
    
  trig4NormalizationMap.emplace(trig4Normalization::CINT7,"CINT7");
  trig4NormalizationMap.emplace(trig4Normalization::CMLL7,"CMLL7");
  trig4NormalizationMap.emplace(trig4Normalization::CMUL7,"CMUL7");
  trig4NormalizationMap.emplace(trig4Normalization::CMSL7,"CMSL7");
  trig4NormalizationMap.emplace(trig4Normalization::CMUL7orCMLL7,"CMUL7orCMLL7");
  trig4NormalizationMap.emplace(trig4Normalization::CMUL7andCMLL7,"CMUL7andCMLL7");
  trig4NormalizationMap.emplace(trig4Normalization::CINT7and0MULor0MLL,"CINT7and0MULor0MLL");
  trig4NormalizationMap.emplace(trig4Normalization::CINT7and0MSL,"CINT7and0MSL");
  trig4NormalizationMap.emplace(trig4Normalization::CMSL7and0MULor0MLL,"CMSL7and0MULor0MLL");
  trig4NormalizationMap.emplace(trig4Normalization::CMSL7and0MLL,"CMSL7and0MLL");
  trig4NormalizationMap.emplace(trig4Normalization::CMSL7and0MUL,"CMSL7and0MUL");  

  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1,  TList::Class());
}

AliAnalysisTaskAODEventStudy::~AliAnalysisTaskAODEventStudy() 
{
  
}

void AliAnalysisTaskAODEventStudy::UserCreateOutputObjects()
{
  
  TList* fOutputList  = new TList();
  fOutputList->SetOwner(true);

  for (auto itrig = trig4NormalizationMap.begin(); itrig != trig4NormalizationMap.end(); ++itrig){
    fHistEvent.push_back(new TH2F(Form("fHistEvent_%s",itrig->second.c_str()),"",2500,0,2500,100,0,100));
    fHistEvent4RunQA.push_back(new TH1F(Form("fHistEvent4RunQA_%s",itrig->second.c_str()),"",2500,0,2500));    
    fOutputList->Add(fHistEvent.back());
    fOutputList->Add(fHistEvent4RunQA.back());
  }
  
  string trigName[]={"CINT7","CMSL7","CMUL7","CMLL7","CMUL7orCMLL7","CMUL7andCMLL7"};

  for(uint itrig=0; itrig<sizeof(trigName)/sizeof(string);++itrig){
    
    fHistNTrackMatchAll[itrig] = new TH1F(Form("fHistNTrackMatchAll_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistNTrackMatchApt[itrig] = new TH1F(Form("fHistNTrackMatchApt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistNTrackMatchLpt[itrig] = new TH1F(Form("fHistNTrackMatchLpt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistNTrackMatchHpt[itrig] = new TH1F(Form("fHistNTrackMatchHpt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistMeanPtMatchAll[itrig] = new TH1F(Form("fHistMeanPtMatchAll_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistMeanPtMatchApt[itrig] = new TH1F(Form("fHistMeanPtMatchApt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistMeanPtMatchLpt[itrig] = new TH1F(Form("fHistMeanPtMatchLpt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    fHistMeanPtMatchHpt[itrig] = new TH1F(Form("fHistMeanPtMatchHpt_%s",trigName[itrig].c_str()),"",2500,0,2500);
    
    fOutputList->Add(fHistNTrackMatchAll[itrig]);
    fOutputList->Add(fHistNTrackMatchApt[itrig]);
    fOutputList->Add(fHistNTrackMatchLpt[itrig]);
    fOutputList->Add(fHistNTrackMatchHpt[itrig]);
    fOutputList->Add(fHistMeanPtMatchAll[itrig]);
    fOutputList->Add(fHistMeanPtMatchApt[itrig]);
    fOutputList->Add(fHistMeanPtMatchLpt[itrig]);
    fOutputList->Add(fHistMeanPtMatchHpt[itrig]);

    fHistPtMatchAll[itrig] = new TH1F(Form("fHistPtMatchAll_%s",trigName[itrig].c_str()),"",500,0,50);
    fHistPtMatchApt[itrig] = new TH1F(Form("fHistPtMatchApt_%s",trigName[itrig].c_str()),"",500,0,50);
    fHistPtMatchLpt[itrig] = new TH1F(Form("fHistPtMatchLpt_%s",trigName[itrig].c_str()),"",500,0,50);
    fHistPtMatchHpt[itrig] = new TH1F(Form("fHistPtMatchHpt_%s",trigName[itrig].c_str()),"",500,0,50);
    fOutputList->Add(fHistPtMatchAll[itrig]);
    fOutputList->Add(fHistPtMatchApt[itrig]);
    fOutputList->Add(fHistPtMatchLpt[itrig]);
    fOutputList->Add(fHistPtMatchHpt[itrig]);
    
    fHistCentrality[itrig] = new TH1F(Form("fHistCentrality_%s",trigName[itrig].c_str()),"",100,0,100);
    fHistCentVtxZ[itrig] = new TH2F(Form("fHistCentVtxZ_%s",trigName[itrig].c_str()),"",60,-30,30,100,0,100);  
    fHistEtaPt[itrig] = new TH2F(Form("fHistEtaPt_%s",trigName[itrig].c_str()),"",25,2.0,4.5,500,0,50);  
    fOutputList->Add(fHistCentrality[itrig]);
    fOutputList->Add(fHistCentVtxZ[itrig]);
    fOutputList->Add(fHistEtaPt[itrig]);  

    fHistSPDTrkEta05Centrality[itrig] = new TH2F(Form("fHistSPDTrkEta05Centrality_%s",trigName[itrig].c_str()),"",200,0,200,100,0,100);
    fHistSPDTrkEta10Centrality[itrig] = new TH2F(Form("fHistSPDTrkEta10Centrality_%s",trigName[itrig].c_str()),"",200,0,200,100,0,100);
    fHistSPDTrkEta15Centrality[itrig] = new TH2F(Form("fHistSPDTrkEta15Centrality_%s",trigName[itrig].c_str()),"",200,0,200,100,0,100);
    fHistSPDTrkEta20Centrality[itrig] = new TH2F(Form("fHistSPDTrkEta20Centrality_%s",trigName[itrig].c_str()),"",200,0,200,100,0,100);
    fHistSPDTrkEtaAllCentrality[itrig] = new TH2F(Form("fHistSPDTrkEtaAllCentrality_%s",trigName[itrig].c_str()),"",200,0,200,100,0,100);
    fOutputList->Add(fHistSPDTrkEta05Centrality[itrig]);
    fOutputList->Add(fHistSPDTrkEta10Centrality[itrig]);
    fOutputList->Add(fHistSPDTrkEta15Centrality[itrig]);
    fOutputList->Add(fHistSPDTrkEta20Centrality[itrig]);
    fOutputList->Add(fHistSPDTrkEtaAllCentrality[itrig]);

    fHistChV0ACentrality[itrig] = new TH2F(Form("fHistChV0ACentrality_%s",trigName[itrig].c_str()),"",50,0,500,100,0,100);
    fHistChV0CCentrality[itrig] = new TH2F(Form("fHistChV0CCentrality_%s",trigName[itrig].c_str()),"",50,0,500,100,0,100);
    fHistChV0MCentrality[itrig] = new TH2F(Form("fHistChV0MCentrality_%s",trigName[itrig].c_str()),"",100,0,1000,100,0,100);
    fOutputList->Add(fHistChV0ACentrality[itrig]);
    fOutputList->Add(fHistChV0CCentrality[itrig]);
    fOutputList->Add(fHistChV0MCentrality[itrig]);
    
    fHistCorrV0MCentrality[itrig] = new TH2F(Form("fHistCorrV0MCentrality_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrV0ACentrality[itrig] = new TH2F(Form("fHistCorrV0ACentrality_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrV0CCentrality[itrig] = new TH2F(Form("fHistCorrV0CCentrality_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrSPDTrkCentrality[itrig] = new TH2F(Form("fHistCorrSPDTrkCentrality_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);

    fOutputList->Add(fHistCorrV0MCentrality[itrig]);
    fOutputList->Add(fHistCorrV0ACentrality[itrig]);
    fOutputList->Add(fHistCorrV0CCentrality[itrig]);
    fOutputList->Add(fHistCorrSPDTrkCentrality[itrig]);

    /*
    fHistCorrChV0AChV0C[itrig] = new TH2F(Form("fHistCorrChV0AChV0C_%s",trigName[itrig].c_str()),"",50,0,500,50,0,500);
    fHistCorrCentSPDTrkV0M[itrig] = new TH2F(Form("fHistCorrCentSPDTrkV0M_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrCentSPDTrkV0A[itrig] = new TH2F(Form("fHistCorrCentSPDTrkV0A_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrCentSPDTrkV0C[itrig] = new TH2F(Form("fHistCorrCentSPDTrkV0C_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fHistCorrCentV0AV0C[itrig] = new TH2F(Form("fHistCorrCentV0AV0C_%s",trigName[itrig].c_str()),"",100,0,100,100,0,100);
    fOutputList->Add(fHistCorrChV0AChV0C[itrig]);
    fOutputList->Add(fHistCorrCentSPDTrkV0M[itrig]);
    fOutputList->Add(fHistCorrCentSPDTrkV0A[itrig]);
    fOutputList->Add(fHistCorrCentSPDTrkV0C[itrig]);
    fOutputList->Add(fHistCorrCentV0AV0C[itrig]);
    */
    fHistSPDTrkSPDClust[itrig] = new TH2F(Form("fHistSPDTrkSPDClust_%s",trigName[itrig].c_str()),"",150,0,150,400,0,400);
    fOutputList->Add(fHistSPDTrkSPDClust[itrig]);
  }
  
  string matchName[]={"All","Apt","Lpt","Hpt"};

  for(uint imatch=0; imatch<sizeof(matchName)/sizeof(string); ++imatch){
    fHistTrackChi2[imatch] = new TH2F(Form("fHistTrackChi2_%s",matchName[imatch].c_str()),"",100,0,50,100,0,10);
    fHistTrigChi2[imatch] = new TH2F(Form("fHistTrigChi2_%s",matchName[imatch].c_str()),"",100,0,50,100,0,10);
    fOutputList->Add(fHistTrackChi2[imatch]);
    fOutputList->Add(fHistTrigChi2[imatch]);
  }

  fHistRabs = new TH1F("fHistRabs","",100,0,100);
  fOutputList->Add(fHistRabs);  
  
  fHistChVtxEta = new TH2F("fHistChVtxEta","",100,-2.5,2.5,60,-30,30);  
  fOutputList->Add(fHistChVtxEta);

  PostData(1,fOutputList); 
}

void AliAnalysisTaskAODEventStudy::UserExec(Option_t *){

  if(!Initialize()) return;  
  if(!fUtils->isAcceptEvent()) return;
  if(fabs(fUtils->getVtxZ()>10.))cout<<fUtils->getVtxZ()<<endl;
  RunQA();
  MultiplicityQA();
  TrackQA();
}

bool AliAnalysisTaskAODEventStudy::Initialize()
{
  
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    
  if( !fUtils->setEvent(fEvent,fInputHandler) )
    return false;
  
  fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7, fIs0MSL, fIs0MSH, fIs0MUL, fIs0MLL);
  
  if( fRunNumber != fEvent->GetRunNumber() ){
    fRunNumber = fUtils->getRunnumber();
    AliMuonTrackCuts* trackCut = fUtils->getMuonTrackCuts();
    trackCut->SetRun(fInputHandler);
  }
  
  fNTrackMatchAll=0;
  fNTrackMatchApt=0;
  fNTrackMatchLpt=0;
  fNTrackMatchHpt=0;

  return true;

}

void AliAnalysisTaskAODEventStudy::RunQA()
{
  if(fIsCINT7){
    fHistEvent[trig4Normalization::CINT7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CINT7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }  
  
  if(fIsCMSL7){
    fHistEvent[trig4Normalization::CMSL7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMSL7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }
  
  if(fIsCMUL7){
    fHistEvent[trig4Normalization::CMUL7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMUL7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }
  
  if(fIsCMLL7){
    fHistEvent[trig4Normalization::CMLL7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMLL7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }  
  
  if(fIsCMUL7 || fIsCMLL7){
    fHistEvent[trig4Normalization::CMUL7orCMLL7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass()); 
    fHistEvent4RunQA[trig4Normalization::CMUL7orCMLL7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }
  
  if(fIsCMUL7 && fIsCMLL7){
    fHistEvent[trig4Normalization::CMUL7andCMLL7] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMUL7andCMLL7] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }
  if( fIsCINT7 ){
    if( fIs0MUL || fIs0MLL ){
      fHistEvent[trig4Normalization::CINT7and0MULor0MLL] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
      fHistEvent4RunQA[trig4Normalization::CINT7and0MULor0MLL] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
    }
  }  
  if( fIsCINT7 && fIs0MSL ){
    fHistEvent[trig4Normalization::CINT7and0MSL] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CINT7and0MSL] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }  
  if( fIsCMSL7 ){
    if( fIs0MUL || fIs0MLL ){
      fHistEvent[trig4Normalization::CMSL7and0MULor0MLL] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
      fHistEvent4RunQA[trig4Normalization::CMSL7and0MULor0MLL] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
    }
  }
  if( fIsCMSL7 && fIs0MLL ){
    fHistEvent[trig4Normalization::CMSL7and0MLL] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMSL7and0MLL] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }
  if( fIsCMSL7 && fIs0MUL ){
    fHistEvent[trig4Normalization::CMSL7and0MUL] -> Fill(fUtils->getRunnumberIndex(),fUtils->getCentClass());
    fHistEvent4RunQA[trig4Normalization::CMSL7and0MUL] -> Fill(fUtils->getRunnumberIndex(),1./fUtils->getDS());
  }

}

void AliAnalysisTaskAODEventStudy::FillMultiplicityQA(int index)
{
  fHistCentrality[index]->Fill(fUtils->getCentClass());
  fHistCentVtxZ[index]->Fill(fUtils->getVtxZ(),fUtils->getCentClass());
  fHistSPDTrkEta05Centrality[index]->Fill(fUtils->getNSPDTrkInfo(0),fUtils->getCentClass());
  fHistSPDTrkEta10Centrality[index]->Fill(fUtils->getNSPDTrkInfo(1),fUtils->getCentClass());
  fHistSPDTrkEta15Centrality[index]->Fill(fUtils->getNSPDTrkInfo(2),fUtils->getCentClass());
  fHistSPDTrkEta20Centrality[index]->Fill(fUtils->getNSPDTrkInfo(3),fUtils->getCentClass());
  fHistSPDTrkEtaAllCentrality[index]->Fill(fUtils->getNSPDTrkInfo(4),fUtils->getCentClass());
  fHistChV0ACentrality[index]->Fill(fUtils->getVzeroInfo(0),fUtils->getCentClass());
  fHistChV0CCentrality[index]->Fill(fUtils->getVzeroInfo(1),fUtils->getCentClass());
  fHistChV0MCentrality[index]->Fill(fUtils->getVzeroInfo(2),fUtils->getCentClass());
  fHistCorrSPDTrkCentrality[index]->Fill(fUtils->getCentClass(0),fUtils->getCentClass());
  fHistCorrV0MCentrality[index]->Fill(fUtils->getCentClass(1),fUtils->getCentClass());
  fHistCorrV0ACentrality[index]->Fill(fUtils->getCentClass(2),fUtils->getCentClass());
  fHistCorrV0CCentrality[index]->Fill(fUtils->getCentClass(3),fUtils->getCentClass());
}

void AliAnalysisTaskAODEventStudy::MultiplicityQA()
{
  
  if( fIsCINT7 ){
    FillMultiplicityQA(0);
  }
  if( fIsCMSL7 ){ 
    FillMultiplicityQA(1);
  }
  if( fIsCMUL7 ){ 
    FillMultiplicityQA(2);
  }
  if( fIsCMLL7 ){ 
    FillMultiplicityQA(3);
  }
  if( fIsCMUL7 || fIsCMLL7 ){ 
    FillMultiplicityQA(4);
  }
  if( fIsCMUL7 && fIsCMLL7 ){ 
    FillMultiplicityQA(5);
  }

  if(!(fIsCMSL7 || fIsCMUL7 || fIsCMLL7)) return;
  
  AliAODTracklets *tracklet = (AliAODTracklets*)fEvent->GetTracklets(); 

  if(!tracklet) return;

  if(tracklet){
    for(Int_t iTrk=0; iTrk < tracklet->GetNumberOfTracklets(); ++iTrk){
      Double_t eta = -1*TMath::Log(TMath::Tan(tracklet->GetTheta(iTrk)/2.));
      Double_t phi = tracklet->GetPhi(iTrk);
      fHistChVtxEta->Fill(eta,fUtils->getVtxZ());
    }
  }
  
}

void AliAnalysisTaskAODEventStudy::FillTrackQA(int index,AliAODTrack *track)
{
  fHistPtMatchAll[index]->Fill(track->Pt());
  fHistNTrackMatchAll[index]->Fill(fUtils->getRunnumberIndex());
  fHistMeanPtMatchAll[index]->Fill(fUtils->getRunnumberIndex(),track->Pt());
  if(track->GetMatchTrigger()>0){
    fHistPtMatchApt[index]->Fill(track->Pt());	
    fHistNTrackMatchApt[index]->Fill(fUtils->getRunnumberIndex());
    fHistMeanPtMatchApt[index]->Fill(fUtils->getRunnumberIndex(),track->Pt());
  }
  if(track->GetMatchTrigger()>1){
    fHistPtMatchLpt[index]->Fill(track->Pt());
    fHistEtaPt[index]->Fill(fabs(track->Eta()),track->Pt());
    fHistNTrackMatchLpt[index]->Fill(fUtils->getRunnumberIndex());
    fHistMeanPtMatchLpt[index]->Fill(fUtils->getRunnumberIndex(),track->Pt());
  }
  if(track->GetMatchTrigger()>2){
    fHistPtMatchHpt[index]->Fill(track->Pt());
    fHistNTrackMatchHpt[index]->Fill(fUtils->getRunnumberIndex());
    fHistMeanPtMatchHpt[index]->Fill(fUtils->getRunnumberIndex(),track->Pt());
  }
}

void AliAnalysisTaskAODEventStudy::TrackQA()
{
  
  int nTrack = fEvent->GetNumberOfTracks();
  
  for(int iTrack=0; iTrack<nTrack; ++iTrack){
    
    AliAODTrack* track = (AliAODTrack*)fEvent->GetTrack(iTrack);
    
    if(!fUtils->isAcceptMuonTrack(track)) continue;
 
    if( fIsCINT7 ){ 
      FillTrackQA(0,track);
    }

    if( fIsCMSL7 ){ 
      FillTrackQA(1,track);
    }

    if( fIsCMUL7 ){ 
      FillTrackQA(2,track);
    }

    if( fIsCMLL7 ){ 
      FillTrackQA(3,track);
    }

    if( fIsCMUL7 || fIsCMLL7 ){ 
      FillTrackQA(4,track);
    }

    if( fIsCMUL7 && fIsCMLL7 ){ 
      FillTrackQA(5,track);
    }
        
    if( fIsCMSL7 ||  fIsCMUL7 || fIsCMLL7 ){
      fHistTrackChi2[0]->Fill(track->Pt(),track->Chi2perNDF());
      fHistTrigChi2[0]->Fill(track->Pt(),track->GetChi2MatchTrigger());
      if(track->GetMatchTrigger()>0){
	fHistTrackChi2[1]->Fill(track->Pt(),track->Chi2perNDF());
	fHistTrigChi2[1]->Fill(track->Pt(),track->GetChi2MatchTrigger());
      }
      if(track->GetMatchTrigger()>1){
	fHistTrackChi2[2]->Fill(track->Pt(),track->Chi2perNDF());
	fHistTrigChi2[2]->Fill(track->Pt(),track->GetChi2MatchTrigger());
      }
      if(track->GetMatchTrigger()>2){
	fHistTrackChi2[3]->Fill(track->Pt(),track->Chi2perNDF());
	fHistTrigChi2[3]->Fill(track->Pt(),track->GetChi2MatchTrigger());
      }

      if(track->GetMatchTrigger()>1){
	fHistRabs->Fill(track->GetRAtAbsorberEnd());
      }
    }

  }


  //isAcceptMuonTrack;

}
