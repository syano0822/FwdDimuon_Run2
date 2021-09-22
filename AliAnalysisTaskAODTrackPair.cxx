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
#include "AliAnalysisTaskAODTrackPair.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)       

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair() : 
  AliAnalysisTaskSE(),
  fEvent(NULL),
  fPoolMuonTrackMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),  
  fRunNumber(-99999),
  fTrackDepth(1000),
  fPoolSize(1),
  fReadyFraction(0.1),
  fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7),
  fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  fIsCINT7(false),
  fIsCMSL7(false),
  fIsCMSH7(false),
  fIsCMUL7(false),
  fIsCMLL7(false),
  fOutputList(NULL),
  fEventCounter(NULL),
  fSparseULSDimuon(NULL),
  fSparseLSppDimuon(NULL),
  fSparseLSmmDimuon(NULL),
  fSparseMixULSDimuon(NULL),
  fSparseMixLSppDimuon(NULL),
  fSparseMixLSmmDimuon(NULL)
{ 
  
}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char* name) : 
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fPoolMuonTrackMgr(NULL),
  fUtils(NULL),
  fMCTrackArray(NULL),
  fIsMC(false),
  fRunNumber(-99999),
  fTrackDepth(1000),
  fPoolSize(1),
  fReadyFraction(0.1),
  fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7),
  fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  fIsCINT7(false),
  fIsCMSL7(false),
  fIsCMSH7(false),
  fIsCMUL7(false),
  fIsCMLL7(false),
  fOutputList(NULL),
  fEventCounter(NULL),
  fSparseULSDimuon(NULL),
  fSparseLSppDimuon(NULL),
  fSparseLSmmDimuon(NULL),
  fSparseMixULSDimuon(NULL),
  fSparseMixLSppDimuon(NULL),
  fSparseMixLSmmDimuon(NULL)
{ 
  
  
  double fCentBins[] = {-1,10,20,30,40,50,60,70,80,90,101};
  double fVtxBins[] = {-50,-10.5,-6,-2,0,2,6,10.5,50};
  double fPsiBins[] = {-10,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,10};
  
  int fNCentBins = sizeof(fCentBins)/sizeof(double)-1;
  int fNVtxZBins = sizeof(fVtxBins)/sizeof(double)-1;
  int fNPsiBins = sizeof(fPsiBins)/sizeof(double)-1;

  fPoolMuonTrackMgr = new AliEventPoolManager(fPoolSize,fTrackDepth,
					      fNCentBins,(double*)fCentBins,
					      fNVtxZBins,(double*)fVtxBins,
					      fNPsiBins,(double*)fPsiBins);  
  fPoolMuonTrackMgr->SetTargetValues(fTrackDepth,(double)fReadyFraction,fPoolSize);  
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPair::~AliAnalysisTaskAODTrackPair() 
{
  
}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPair::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  double bins_cent_hist[]={0.,1.,2.,3.,4.,5.,10.,20.,30.,50.,70.,100.};
  int binnum_cent_hist = sizeof(bins_cent_hist)/sizeof(double) - 1;
  
  double bins_event_hist[]={0,1,2,3,4,5,6,7,8,9,10};
  int binnum_event_hist = sizeof(bins_event_hist)/sizeof(double) - 1;

  std::string event_label[]={"CMUL7","CMLL7","CMUL7orCMLL7","CMUL7andCMLL7","CMUL7withDS","CMLL7withDS","CMUL7orCMLL7withDS","CMUL7andCMLL7withDS"};
  fEventCounter = new TH2F("fEventCounter","",binnum_event_hist,bins_event_hist,binnum_cent_hist,bins_cent_hist);
  for(int iname=0; iname<sizeof(event_label)/sizeof(std::string); ++iname)
    fEventCounter->GetXaxis()->SetBinLabel(iname+1,event_label[iname].c_str());
  fOutputList->Add(fEventCounter);

  const Int_t binnum_Cent            = 4;
  Int_t bins_Cent[binnum_Cent]       = {800,6,300,binnum_cent_hist};
  double minbins_Cent[binnum_Cent] = {0,2.5,0,0};
  double maxbins_Cent[binnum_Cent] = {4,4.0,30,100};
  TString  namebins_Cent[]        = {"mass","rap","pt","centrality"};
  
  fSparseULSDimuon = new THnSparseF("fSparseULSDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseLSppDimuon = new THnSparseF("fSparseLSppDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseLSmmDimuon = new THnSparseF("fSparseLSmmDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseULSDimuon->SetBinEdges(3,bins_cent_hist);
  fSparseLSppDimuon->SetBinEdges(3,bins_cent_hist);
  fSparseLSmmDimuon->SetBinEdges(3,bins_cent_hist);
  fOutputList->Add(fSparseULSDimuon);
  fOutputList->Add(fSparseLSppDimuon);
  fOutputList->Add(fSparseLSmmDimuon);
  
  fSparseMixULSDimuon = new THnSparseF("fSparseMixULSDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseMixLSppDimuon = new THnSparseF("fSparseMixLSppDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseMixLSmmDimuon = new THnSparseF("fSparseMixLSmmDimuon","",binnum_Cent,bins_Cent,minbins_Cent,maxbins_Cent);
  fSparseULSDimuon->SetBinEdges(3,bins_cent_hist);
  fSparseLSppDimuon->SetBinEdges(3,bins_cent_hist);
  fSparseLSmmDimuon->SetBinEdges(3,bins_cent_hist);
  fOutputList->Add(fSparseMixULSDimuon);
  fOutputList->Add(fSparseMixLSppDimuon);
  fOutputList->Add(fSparseMixLSmmDimuon);
  
  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{  
  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;
  
  MuonPairAnalysis();
  MuonPairAnalysisEveMixing();
}

bool AliAnalysisTaskAODTrackPair::Initialize()
{

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());

  if( !fUtils->setEvent(fEvent,fInputHandler) )
    return false;

  fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);

  if( fRunNumber != fEvent->GetRunNumber() ){
    fRunNumber = fUtils->getRunnumber();
    AliMuonTrackCuts* trackCut = fUtils->getMuonTrackCuts();
    trackCut->SetRun(fInputHandler);
  }
  
  return true;
}


AliEventPool* AliAnalysisTaskAODTrackPair::setPool()
{
  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi();

  AliEventPool* pool = (AliEventPool*)fPoolMuonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);
  
  return pool;
}

bool AliAnalysisTaskAODTrackPair::MuonPairAnalysisEveMixing(){
  
  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForMixing) ) return false;
  
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  AliEventPool* pool = setPool();  
  
  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);        
    if(!fUtils->isAcceptMuonTrack(track1)) continue;
    
    if (pool->IsReady()){
      
      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){
	
	TObjArray* poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);	
	
	for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){

	  AliAODTrack* __track2__ = (AliAODTrack*)poolTracks->At(iTrack2);
	  
	  AliAODTrack* track2 = (AliAODTrack*)__track2__->Clone();

	  if(!fUtils->isAcceptMuonTrack(track2)) continue;
	  
	  AliAODDimuon* dimuon = new AliAODDimuon();
	  dimuon->SetMuons(track1,track2);

	  if(!fUtils->isAcceptDimuon(dimuon)) continue;
	  
	  double fill[]={dimuon->M(),dimuon->Pt(),fabs(dimuon->Y()),fUtils->getCentClass()};
	  
	  if(dimuon->Charge() == 0){
	    fSparseMixULSDimuon->Fill(fill,(double)1./fUtils->getDS());
	  }
	  else if(dimuon->Charge() > 0){
	    fSparseMixLSppDimuon->Fill(fill,(double)1./fUtils->getDS());
	  }
	  else{
	    fSparseMixLSmmDimuon->Fill(fill,(double)1./fUtils->getDS());
	  }
	  
	  delete track2;
	  delete dimuon;
	  
	}//end of loop track2
      }// end of loop iMixEvt
    }//poolPion->IsReady()    
    
    fTrackArray->Add(track1);
   
  }//end of loop track1

  TObjArray* fTrackArrayClone = (TObjArray*)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if(fTrackArrayClone->GetEntriesFast()>0){
    pool->UpdatePool(fTrackArrayClone);
  }
  
}

bool AliAnalysisTaskAODTrackPair::MuonPairAnalysis()
{ 
  
  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForSame) ) return false;

  if(fIsCMUL7){
    fEventCounter->Fill(1,fUtils->getCentClass());
    fEventCounter->Fill(5,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMLL7){
    fEventCounter->Fill(2,fUtils->getCentClass());
    fEventCounter->Fill(6,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 || fIsCMLL7){
    fEventCounter->Fill(3,fUtils->getCentClass());
    fEventCounter->Fill(7,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 && fIsCMLL7){
    fEventCounter->Fill(4,fUtils->getCentClass());
    fEventCounter->Fill(8,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }


  //cout<<fEvent->GetFiredTriggerClasses()<<endl;

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){
   
    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);
    
    if(!fUtils->isAcceptMuonTrack(track1)) continue;

    for(Int_t iTrack2=iTrack1+1; iTrack2<nTrack; ++iTrack2){
   
      AliAODTrack *track2 = (AliAODTrack*)fEvent->GetTrack(iTrack2);
      
      if(!fUtils->isAcceptMuonTrack(track2)) continue;
      
      AliAODDimuon* dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1,track2);

      if(!fUtils->isAcceptDimuon(dimuon)) continue;
      
      double fill[]={dimuon->M(),dimuon->Pt(),fabs(dimuon->Y()),fUtils->getCentClass()};
      
      if(dimuon->Charge() == 0){
	fSparseULSDimuon->Fill(fill,(double)1./fUtils->getDS());
      }
      else if(dimuon->Charge() > 0){
	fSparseLSppDimuon->Fill(fill,(double)1./fUtils->getDS());
      }
      else{
	fSparseLSmmDimuon->Fill(fill,(double)1./fUtils->getDS());
      }

      delete dimuon;

    }//end of loop track2
  }//end of loop track1 
}
