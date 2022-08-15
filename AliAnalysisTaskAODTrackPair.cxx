#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TRefArray.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliEventPoolMuon.h"
#include "AliEventPoolManager.h"
#include "AliGenEventHeader.h"

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
  fIsMC(false),
  fIsMidMuonAna(false),
  fIsMixingAnalysis(false),
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

  fHistTrackPairPtBalance(NULL),
  fHistTrackPairLocalBoardPair(NULL),

  fHistTrackEta(NULL),
  fHistTrackThetaAbs(NULL),
  fHistTrackTriggerMatch(NULL),
  fHistTrackPDCA(NULL),
  fHistTrackChiSquare(NULL),
  fHistTriggerChiSquare(NULL),

  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),

  fTreeULSDimuon(NULL),
  fTreeLSppDimuon(NULL),
  fTreeLSmmDimuon(NULL),

  fTreeMixULSDimuon(NULL),
  fTreeMixLSppDimuon(NULL),
  fTreeMixLSmmDimuon(NULL),

  RecDimuonPt(0.),
  RecDimuonRap(0.),
  RecDimuonMass(0.),
  RecDimuonCent(0.),
  RecDimuonDS(0.),

  fTreeMidMuon(NULL),
  fTrackPt(0.),
  fTrackP(0.),
  fTrackTheta(0.),
  fTrackPhi(0.),
  fTrackLength(0.),
  fTrackBeta(0.),
  fTrackTrackChi2perNDF(0.),
  fTrackTrackITSNcls(0.),
  fTrackTrackTPCNcls(0.),
  fTrackTrackTOFNcls(0.),
  fTrackTrackTPCChi2(0.),
  fTrackTrackITSChi2(0.),
  fTrackTPCCrossedRows(0.),
  fTrackTPCFindableNcls(0.),
  fTrackTOFBCTime(0.),
  fTrackTOFKinkIndex(0.),
  fTrackDCAxy(0.),
  fTrackDCAz(0.),
  fTrackTPCsigmaMuon(0.),
  fTrackTOFsigmaMuon(0.),
  fTrackTrueMuonLabel(false)
{

}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char* name) :
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fPoolMuonTrackMgr(NULL),
  fUtils(NULL),
  fIsMC(false),
  fIsMidMuonAna(false),
  fIsMixingAnalysis(false),
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

  fHistTrackPairPtBalance(NULL),
  fHistTrackPairLocalBoardPair(NULL),

  fHistTrackEta(NULL),
  fHistTrackThetaAbs(NULL),
  fHistTrackTriggerMatch(NULL),
  fHistTrackPDCA(NULL),
  fHistTrackChiSquare(NULL),
  fHistTriggerChiSquare(NULL),

  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),

  fTreeULSDimuon(NULL),
  fTreeLSppDimuon(NULL),
  fTreeLSmmDimuon(NULL),

  fTreeMixULSDimuon(NULL),
  fTreeMixLSppDimuon(NULL),
  fTreeMixLSmmDimuon(NULL),

  RecDimuonPt(0.),
  RecDimuonRap(0.),
  RecDimuonMass(0.),
  RecDimuonCent(0.),
  RecDimuonDS(0.),

  fTreeMidMuon(NULL),
  fTrackPt(0.),
  fTrackP(0.),
  fTrackTheta(0.),
  fTrackPhi(0.),
  fTrackLength(0.),
  fTrackBeta(0.),
  fTrackTrackChi2perNDF(0.),
  fTrackTrackITSNcls(0.),
  fTrackTrackTPCNcls(0.),
  fTrackTrackTOFNcls(0.),
  fTrackTrackTPCChi2(0.),
  fTrackTrackITSChi2(0.),
  fTrackTPCCrossedRows(0.),
  fTrackTPCFindableNcls(0.),
  fTrackTOFBCTime(0.),
  fTrackTOFKinkIndex(0.),
  fTrackDCAxy(0.),
  fTrackDCAz(0.),
  fTrackTPCsigmaMuon(0.),
  fTrackTOFsigmaMuon(0.),
  fTrackTrueMuonLabel(false)
{

  double fCentBins[] = {-1,9,15,21,26,34,42,51,61,99999};
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

  double bins_event_hist[]={0,1,2,3,4,5,6,7,8,9,10};
  int binnum_event_hist = sizeof(bins_event_hist)/sizeof(double) - 1;

  std::string event_label[]=
    {"CMUL7","CMLL7","CMUL7orCMLL7","CMUL7andCMLL7","CMUL7withDS","CMLL7withDS","CMUL7orCMLL7withDS","CMUL7andCMLL7withDS"};
  fEventCounter = new TH2F("fEventCounter","",11,0,11,200,0,200);
  for(unsigned int iname=0; iname<sizeof(event_label)/sizeof(std::string); ++iname) {
    fEventCounter->GetXaxis()->SetBinLabel(iname+1,event_label[iname].c_str());
  }

  fOutputList->Add(fEventCounter);

  if ( !fIsMixingAnalysis ){
    fTreeULSDimuon = new TTree("fTreeULSDimuon","");
    fTreeULSDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeULSDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeULSDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeULSDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeULSDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeULSDimuon);

    fTreeLSppDimuon = new TTree("fTreeLSppDimuon","");
    fTreeLSppDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeLSppDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeLSppDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeLSppDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeLSppDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeLSppDimuon);

    fTreeLSmmDimuon = new TTree("fTreeLSmmDimuon","");
    fTreeLSmmDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeLSmmDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeLSmmDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeLSmmDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeLSmmDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeLSmmDimuon);
  } else {
    fTreeMixULSDimuon = new TTree("fTreeMixULSDimuon","");
    fTreeMixULSDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeMixULSDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeMixULSDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeMixULSDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeMixULSDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeMixULSDimuon);

    fTreeMixLSppDimuon = new TTree("fTreeMixLSppDimuon","");
    fTreeMixLSppDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeMixLSppDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeMixLSppDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeMixLSppDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeMixLSppDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeMixLSppDimuon);

    fTreeMixLSmmDimuon = new TTree("fTreeMixLSmmDimuon","");
    fTreeMixLSmmDimuon->Branch("RecDimuonPt",&RecDimuonPt,"RecDimuonPt/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonRap",&RecDimuonRap,"RecDimuonRap/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonMass",&RecDimuonMass,"RecDimuonMass/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonCent",&RecDimuonCent,"RecDimuonCent/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonDS",&RecDimuonDS,"RecDimuonDS/F");
    fOutputList->Add(fTreeMixLSmmDimuon);
  }
  
  fHistTrackPairPtBalance = new TH2F("fHistTrackPairPtBalance","",50,0,5,50,0,5);
  fHistTrackPairLocalBoardPair = new TH2F("fHistTrackPairLocalBoardPair","",240,0,240,240,0,240);
  fOutputList->Add(fHistTrackPairPtBalance);
  fOutputList->Add(fHistTrackPairLocalBoardPair);

  fHistEventVtxZ = new TH1F("fHistEventVtxZ","",60,-30,30);
  fHistEventCent = new TH1F("fHistEventCent","",100,0,100);
  fHistEventMulti = new TH1F("fHistEventMulti","",200,0,200);
  fHistEventVtxCont = new TH1F("fHistEventVtxCont","",100,0,100);
  fOutputList->Add(fHistEventVtxZ);
  fOutputList->Add(fHistEventCent);
  fOutputList->Add(fHistEventMulti);
  fOutputList->Add(fHistEventVtxCont);

  if (!fIsMidMuonAna){
    fHistTrackEta = new TH2F("fHistTrackEta","",20,0,10,25,-4.5,-2.0);
    fHistTrackThetaAbs = new TH2F("fHistTrackThetaAbs","",20,0,10,60,0,15);
    fHistTrackTriggerMatch = new TH2F("fHistTrackTriggerMatch","",20,0,10,5,0,5);
    fHistTrackPDCA = new TH2F("fHistTrackPDCA","",20,0,10,200,0,20);
    fHistTrackChiSquare = new TH2F("fHistTrackChiSquare","",20,0,10,100,0,10);
    fHistTriggerChiSquare = new TH2F("fHistTriggerChiSquare","",20,0,10,100,0,10);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTrackThetaAbs);
    fOutputList->Add(fHistTrackTriggerMatch);
    fOutputList->Add(fHistTrackPDCA);
    fOutputList->Add(fHistTrackChiSquare);
    fOutputList->Add(fHistTriggerChiSquare);
  }
  
  if (fIsMidMuonAna) {
    fTreeMidMuon = new TTree("fTreeMidMuon","Tree for machine leraning for PID");
    fTreeMidMuon->Branch("fTrackPt", &fTrackPt, "fTrackPt/F");
    fTreeMidMuon->Branch("fTrackP", &fTrackP, "fTrackP/F");
    fTreeMidMuon->Branch("fTrackTheta", &fTrackTheta, "fTrackTheta/F");
    fTreeMidMuon->Branch("fTrackPhi", &fTrackPhi, "fTrackPhi/F");
    fTreeMidMuon->Branch("fTrackLength", &fTrackLength, "fTrackLength/F");
    fTreeMidMuon->Branch("fTrackBeta", &fTrackBeta, "fTrackBeta/F");
    fTreeMidMuon->Branch("fTrackTrackChi2perNDF", &fTrackTrackChi2perNDF, "fTrackTrackChi2perNDF/F");
    fTreeMidMuon->Branch("fTrackTrackITSNcls", &fTrackTrackITSNcls, "fTrackTrackITSNcls/F");
    fTreeMidMuon->Branch("fTrackTrackTPCNcls", &fTrackTrackTPCNcls, "fTrackTrackTPCNcls/F");
    fTreeMidMuon->Branch("fTrackTrackTOFNcls", &fTrackTrackTOFNcls, "fTrackTrackTOFNcls/F");
    fTreeMidMuon->Branch("fTrackTrackTPCChi2", &fTrackTrackTPCChi2, "fTrackTrackTPCChi2/F");
    fTreeMidMuon->Branch("fTrackTrackITSChi2", &fTrackTrackITSChi2, "fTrackTrackITSChi2/F");
    fTreeMidMuon->Branch("fTrackTPCCrossedRows", &fTrackTPCCrossedRows, "fTrackTPCCrossedRows/F");
    fTreeMidMuon->Branch("fTrackTPCFindableNcls", &fTrackTPCFindableNcls, "fTrackTPCFindableNcls/F");
    fTreeMidMuon->Branch("fTrackTOFBCTime", &fTrackTOFBCTime, "fTrackTOFBCTime/F");
    fTreeMidMuon->Branch("fTrackTOFKinkIndex", &fTrackTOFKinkIndex, "fTrackTOFKinkIndex/F");
    fTreeMidMuon->Branch("fTrackDCAxy", &fTrackDCAxy, "fTrackDCAxy/F");
    fTreeMidMuon->Branch("fTrackDCAz", &fTrackDCAz, "fTrackDCAz/F");
    fTreeMidMuon->Branch("fTrackTPCsigmaMuon", &fTrackTPCsigmaMuon, "fTrackTPCsigmaMuon/F");
    fTreeMidMuon->Branch("fTrackTOFsigmaMuon", &fTrackTOFsigmaMuon, "fTrackTOFsigmaMuon/F");
    fTreeMidMuon->Branch("fTrackTrueMuonLabel", &fTrackTrueMuonLabel, "fTrackTrueMuonLabel/B");
    fOutputList->Add(fTreeMidMuon);
  }
    
  PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{

  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;

  EventQA();
  
  if(!fIsMidMuonAna) {
    if ( !fIsMixingAnalysis ) {
      FwdMuonPairAnalysis();
    } else {
      FwdMuonPairAnalysisEveMixing();
    }
  } else {
    if ( !fIsMixingAnalysis ) {
      MidMuonPairAnalysis();      
    } else {
      //MidMuonPairAnalysisEveMixing();
    }
  }

}

bool AliAnalysisTaskAODTrackPair::Initialize() {
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if( !fUtils->setEvent(fEvent,fInputHandler) )
    return false;
  if( fRunNumber != fEvent->GetRunNumber() ){
    fRunNumber = fUtils->getRunnumber();
    AliMuonTrackCuts* trackCut = fUtils->getMuonTrackCuts();
    trackCut->SetRun(fInputHandler);
  }
  fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);
  return true;
}

bool AliAnalysisTaskAODTrackPair::EventQA() {
  fHistEventVtxZ->Fill(fUtils->getVtxZ());
  fHistEventCent->Fill(fUtils->getCentClass());
  fHistEventMulti->Fill(fUtils->getNCorrSPDTrkInfo(1));
  fHistEventVtxCont->Fill(fUtils->getVtxCont());
  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonTrackQA(AliAODTrack* track){
  fHistTrackEta->Fill(track->Pt(),track->Eta());
  fHistTrackThetaAbs->Fill(track->Pt(),AliAnalysisMuonUtility::GetThetaAbsDeg(track));
  fHistTrackTriggerMatch->Fill(track->Pt(),AliAnalysisMuonUtility::GetMatchTrigger(track));
  fHistTrackChiSquare->Fill(track->Pt(),AliAnalysisMuonUtility::GetChi2perNDFtracker(track));
  fHistTriggerChiSquare->Fill(track->Pt(),AliAnalysisMuonUtility::GetChi2MatchTrigger(track));
  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairQA(AliAODDimuon* dimuon){

  AliAODTrack* track1 = dynamic_cast<AliAODTrack*>(dimuon->GetMu(0));
  AliAODTrack* track2 = dynamic_cast<AliAODTrack*>(dimuon->GetMu(1));

  int triggerLB1 = AliAnalysisMuonUtility::GetLoCircuit(track1);
  int triggerLB2 = AliAnalysisMuonUtility::GetLoCircuit(track2);
  
  float pt_max = track1->Pt();
  float pt_min = track2->Pt();

  if ( track1->Pt() > track2->Pt() ) {
    pt_max = track1->Pt();
    pt_min = track2->Pt();
  } else {
    pt_max = track2->Pt();
    pt_min = track1->Pt();
  }

  fHistTrackPairPtBalance->Fill(pt_max,pt_min);
  fHistTrackPairLocalBoardPair->Fill(triggerLB1,triggerLB2);

  return true;
}


bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysisEveMixing(){

  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForMixing) ) return false;

  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ){
    poolVtxZ=fUtils->getVtxZ();
  }
  if(onEvtMixingPoolCent){
    //poolCent=fUtils->getCentClass();
    poolCent=fUtils->getNCorrSPDTrkInfo(1);
  }
  if(onEvtMixingPoolPsi){
    poolPsi=fUtils->getPsi();
  }

  AliEventPool* pool = (AliEventPool*)fPoolMuonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);

    if(!fUtils->isAcceptFwdMuonTrack(track1)) continue;

    if (pool->IsReady()){

      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){

	TObjArray* poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);

	for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){

	  AliAODTrack* __track2__ = (AliAODTrack*)poolTracks->At(iTrack2);

	  AliAODTrack* track2 = (AliAODTrack*)__track2__->Clone();

	  if(!fUtils->isAcceptFwdMuonTrack(track2)) continue;

	  AliAODDimuon* dimuon = new AliAODDimuon();
	  dimuon->SetMuons(track1,track2);

	  if(!fUtils->isAcceptFwdDimuon(dimuon)) continue;

	  RecDimuonPt = dimuon->Pt();
	  RecDimuonRap = fabs(dimuon->Y());
	  RecDimuonMass = dimuon->M();
	  //RecDimuonCent = fUtils->getCentClass();
	  RecDimuonCent = fUtils->getNCorrSPDTrkInfo(1);
	  RecDimuonDS = fUtils->getDS();

	  string fFiredTrigName = string(fEvent->GetFiredTriggerClasses());

	  if(dimuon->Charge() == 0) {
	    fTreeMixULSDimuon->Fill();
	  } else if(dimuon->Charge() > 0) {
	    fTreeMixLSppDimuon->Fill();
	  } else {
	    fTreeMixLSmmDimuon->Fill();
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

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysis()
{

  if(!fIsMC && !(fInputHandler->IsEventSelected() & fTriggerMaskForSame) ){
    return false;
  }

  if(fIsCMUL7){
    fEventCounter->Fill(0.,fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(4.,fUtils->getNCorrSPDTrkInfo(1),(double)1./fUtils->getDS());
  }
  if(fIsCMLL7){
    fEventCounter->Fill(1.,fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(5.,fUtils->getNCorrSPDTrkInfo(1),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 || fIsCMLL7){
    fEventCounter->Fill(2.,fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(6.,fUtils->getNCorrSPDTrkInfo(1),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 && fIsCMLL7){
    fEventCounter->Fill(3.,fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(7.,fUtils->getNCorrSPDTrkInfo(1),(double)1./fUtils->getDS());
  }

  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack* track1;
  AliAODTrack* track2;

  AliAODDimuon* dimuon;

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);
    
    if(!fUtils->isAcceptFwdMuonTrack(track1)) continue;

    FwdMuonTrackQA(track1);

    for(Int_t iTrack2=iTrack1+1; iTrack2<nTrack; ++iTrack2){

      track2 = (AliAODTrack*)fEvent->GetTrack(iTrack2);

      if(!fUtils->isAcceptFwdMuonTrack(track2)) continue;

      dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1,track2);

      if(!fUtils->isAcceptFwdDimuon(dimuon)) continue;

      FwdMuonPairQA(dimuon);
      
      RecDimuonPt = dimuon->Pt();
      RecDimuonRap = fabs(dimuon->Y());
      RecDimuonMass = dimuon->M();
      //RecDimuonCent = fUtils->getCentClass();
      RecDimuonCent = fUtils->getNCorrSPDTrkInfo(1);
      RecDimuonDS = fUtils->getDS();

      if(dimuon->Charge() == 0) {
	fTreeULSDimuon->Fill();
      } else if(dimuon->Charge() > 0) {
	fTreeLSppDimuon->Fill();
      } else {
	fTreeLSmmDimuon->Fill();
      }

      delete dimuon;

    }//end of loop track2
  }//end of loop track1
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonTrackQA(AliAODTrack* track){  

  float sigTOF =track->GetTOFsignal();
  float length =track->GetIntegratedLength();
  float beta =(sigTOF>0) ? (double)length/ (2.99792457999999984e-02 * sigTOF) : 0;

  float dca_xy=9999;
  float dca_z=9999;
  track->GetImpactParameters(dca_xy,dca_z);
    
  AliTOFHeader * tofHeader = (AliTOFHeader*)track->GetTOFHeader();

  fTrackPt = track->Pt();
  fTrackP = track->P();
  fTrackTheta = track->Theta();
  fTrackPhi = track->Phi();
  fTrackLength = track->GetIntegratedLength();
  fTrackBeta = beta;
  fTrackTrackChi2perNDF = track->Chi2perNDF();
  fTrackTrackITSNcls = track->GetITSNcls();
  fTrackTrackTPCNcls = track->GetTPCNcls();
  fTrackTrackTOFNcls = tofHeader->GetNumberOfTOFclusters();
  fTrackTrackTPCChi2 = track->GetTPCchi2();
  fTrackTrackITSChi2 = track->GetITSchi2();
  fTrackTPCCrossedRows = track->GetTPCCrossedRows();
  fTrackTPCFindableNcls = track->GetTPCNclsF();
  fTrackTOFBCTime = track->GetTOFBunchCrossing();
  fTrackTOFKinkIndex = track->GetKinkIndex(0);
  fTrackDCAxy = dca_xy;
  fTrackDCAz = dca_z;
  fTrackTPCsigmaMuon = fUtils->getTPCSigma(track,AliPID::kMuon);
  fTrackTOFsigmaMuon = fUtils->getTOFSigma(track,AliPID::kMuon);
  
  fTrackTrueMuonLabel = false;
  
  fTreeMidMuon->Fill();
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonPairQA(AliAODDimuon* dimuon){  
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonPairAnalysis()
{
  cout<<"Processing mid-rapidity muon analysis"<<endl;
  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack* track1;
  AliAODTrack* track2;

  AliAODDimuon* dimuon;

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);

    if(!fUtils->isAcceptMidTrackQuality(track1)){
      continue;
    }
    
    MidMuonTrackQA(track1);

    if(!fUtils->isAcceptTrackKinematics(track1)){
      continue;
    }
    if(!fUtils->isAcceptMidMuonTrack(track1)){
      continue;
    }

    for(Int_t iTrack2=iTrack1+1; iTrack2<nTrack; ++iTrack2){

      track2 = (AliAODTrack*)fEvent->GetTrack(iTrack2);

      if(!fUtils->isAcceptTrackKinematics(track2)){
	continue;
      }
      if(!fUtils->isAcceptMidMuonTrack(track2)) {
	continue;
      }
      
      dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1,track2);
      
      if(!fUtils->isAcceptMidDimuon(dimuon)){
	continue;
      }
      
      MidMuonPairQA(dimuon);
      
      double fill[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fUtils->getCentClass()};
      
      RecDimuonPt = dimuon->Pt();
      RecDimuonRap = fabs(dimuon->Y());
      RecDimuonMass = dimuon->M();
      RecDimuonCent = fUtils->getCentClass();
      RecDimuonDS = 1.;

      if(dimuon->Charge() == 0) {
	fTreeULSDimuon->Fill();
      } else if(dimuon->Charge() > 0) {
	fTreeLSppDimuon->Fill();
      } else {
	fTreeLSmmDimuon->Fill();
      }

      delete dimuon;

    }//end of loop track2
  }//end of loop track1
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonPairAnalysisEveMixing(){
  
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray -> SetOwner();

  float poolCent=0.;
  float poolVtxZ=0.;
  float poolPsi=0.;

  if(onEvtMixingPoolVtxZ) poolVtxZ=fUtils->getVtxZ();
  if(onEvtMixingPoolCent) poolCent=fUtils->getCentClass();
  if(onEvtMixingPoolPsi) poolPsi=fUtils->getPsi();

  AliEventPool* pool = (AliEventPool*)fPoolMuonTrackMgr -> GetEventPool(poolCent,poolVtxZ,poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){

    AliAODTrack *track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);

    if(!fUtils->isAcceptMidMuonTrack(track1)) continue;

    if (pool->IsReady()){

      for (Int_t iMixEvt=0; iMixEvt<pool->GetCurrentNEvents(); iMixEvt++){

	TObjArray* poolTracks = (TObjArray*)pool->GetEvent(iMixEvt);

	for(Int_t iTrack2=0; iTrack2<poolTracks->GetEntriesFast(); ++iTrack2){

	  AliAODTrack* __track2__ = (AliAODTrack*)poolTracks->At(iTrack2);

	  AliAODTrack* track2 = (AliAODTrack*)__track2__->Clone();

	  if(!fUtils->isAcceptMidMuonTrack(track2)) continue;

	  AliAODDimuon* dimuon = new AliAODDimuon();
	  dimuon->SetMuons(track1,track2);

	  if(!fUtils->isAcceptMidDimuon(dimuon)) continue;

	  double fill[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fUtils->getCentClass()};

	  RecDimuonPt = dimuon->Pt();
	  RecDimuonRap = fabs(dimuon->Y());
	  RecDimuonMass = dimuon->M();
	  RecDimuonCent = fUtils->getCentClass();
	  RecDimuonDS = fUtils->getDS();
	  
	  if(dimuon->Charge() == 0) {
	    fTreeMixULSDimuon->Fill();
	  } else if(dimuon->Charge() > 0) {
	    fTreeMixLSppDimuon->Fill();
	  } else {
	    fTreeMixLSmmDimuon->Fill();
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

  return true;
}
