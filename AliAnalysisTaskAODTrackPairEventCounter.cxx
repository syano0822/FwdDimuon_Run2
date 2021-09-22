#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TRefArray.h"
#include "TLorentzVector.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolMuon.h"
#include "AliEventPoolManager.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliPIDCombined.h"
#include "AliAODDimuon.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliAODPid.h"
#include "AliTOFHeader.h"

#include "AliMultSelection.h"
#include "AliAnalysisMuonUtility.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskAODTrackPairEventCounter.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPairEventCounter)       

//________________________________________________________________________
AliAnalysisTaskAODTrackPairEventCounter::AliAnalysisTaskAODTrackPairEventCounter() : AliAnalysisTaskSE(),
  fIsMC(NULL),
  fOutputList(NULL),
  fPIDResponse(NULL),
  fEvent(NULL),
  fHandler(NULL),
  fMCTrackArray(NULL),
  fMultSelection(NULL),
  fCent(0.),
  fVtxZ(0.),
  fPsi(0.),
  fHistCentrality(NULL),
  fHistPsi(NULL),
  fHistAllPt(NULL),
  fHistPionPt(NULL),
  fHistKaonPt(NULL)
{
  
  
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPairEventCounter::AliAnalysisTaskAODTrackPairEventCounter(const char *name) : AliAnalysisTaskSE(name),
												     fIsMC(NULL),
												     fOutputList(NULL),
												     fPIDResponse(NULL),
												     fEvent(NULL),
												     fHandler(NULL),
												     fMCTrackArray(NULL),
												     fMultSelection(NULL),
												     fCent(0.),
												     fVtxZ(0.),
												     fPsi(0.),
												     fHistCentrality(NULL),
												     fHistPsi(NULL),
												     fHistAllPt(NULL),
												     fHistPionPt(NULL),
  
  fHistKaonPt(NULL)

{ 
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPairEventCounter::~AliAnalysisTaskAODTrackPairEventCounter() 
{
  
}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPairEventCounter::UserCreateOutputObjects()
{
  
  AliAnalysisManager   *man          = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if(inputHandler->GetPIDResponse()) fPIDResponse=inputHandler->GetPIDResponse();
  
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);
  
  fHistCentrality = new TH1F("fHistCentrality","fHistCentrality",100,0,100);  
  fOutputList->Add(fHistCentrality);

  fHistPsi = new TH1F("fHistPsi","fHistPsi",100,-TMath::Pi(),TMath::Pi());  
  fOutputList->Add(fHistPsi);
  
  fHistAllPt = new TH1F("fHistAllPt","",200,0,20);
  fOutputList->Add(fHistAllPt);

  fHistPionPt = new TH1F("fHistPionPt","",200,0,20);
  fOutputList->Add(fHistPionPt);

  fHistKaonPt = new TH1F("fHistKaonPt","",200,0,20);
  fOutputList->Add(fHistKaonPt);

  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPairEventCounter::UserExec(Option_t *){

  if(!InitEvent()) return;
  
  if(!AliAnalysisTaskAODTrackPairUtils::IsGoodEvent(fEvent,fMultSelection)) return;

  fHistCentrality->Fill(fCent);
  fHistPsi->Fill(fPsi);
  
  SingleTrackAnalysis();

  return;
}

Bool_t AliAnalysisTaskAODTrackPairEventCounter::InitEvent()
{

  if(!fPIDResponse) return false;

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fEvent) return false;;
  
  fMultSelection = (AliMultSelection * )fEvent->FindListObject("MultSelection");
  if(!fMultSelection){
    return false;
  }

  fMCTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(fIsMC && !fMCTrackArray){
    return false;
  }
  
  fEvent->SetUniqueID(0);
  
  AliAODHeader* header = (AliAODHeader*)fEvent->GetHeader();
  AliEventplane* eventplane = (AliEventplane*)fEvent->GetEventplane();
  AliCentrality* cent = fEvent->GetCentrality();
  AliAODVertex * vtx  = fEvent->GetPrimaryVertex();

  fVtxZ = vtx->GetZ();
  fPsi = eventplane->GetEventplane("V0",fEvent,2);
  fCent = fMultSelection->GetMultiplicityPercentile("V0M", false);
  
  cout<<fCent<<"        "<<fEvent->GetFiredTriggerClasses()<<endl;

  return true;
}

Bool_t AliAnalysisTaskAODTrackPairEventCounter::SingleTrackAnalysis()
{

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);
    
    fHistAllPt->Fill(track1->Pt());

    if(!AliAnalysisTaskAODTrackPairUtils::IsGoodTrack(track1)) continue;
    
    if(AliAnalysisTaskAODTrackPairUtils::IsGoodPionTrack(track1,fPIDResponse)) fHistPionPt->Fill(track1->Pt());
    if(AliAnalysisTaskAODTrackPairUtils::IsGoodKaonTrack(track1,fPIDResponse)) fHistKaonPt->Fill(track1->Pt());
    
  }
  
  AliAnalysisTaskAODTrackPairUtils::FilterGoodTracksInEvent(fEvent,fPIDResponse);

  return true;
}
