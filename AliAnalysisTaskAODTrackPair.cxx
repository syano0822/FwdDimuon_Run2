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
  fSparseMixLSmmDimuon(NULL),

  fSparseEtaDalitz(NULL),
  fSparseEta2Body(NULL),
  fSparseRho2Body(NULL),
  fSparseOmegaDalitz(NULL),
  fSparseOmega2Body(NULL),
  fSparsePhi2Body(NULL),
  fSparseEtaPrimeDalitz(NULL),

  fSparseEtaDalitzMC(NULL),
  fSparseEta2BodyMC(NULL),
  fSparseRho2BodyMC(NULL),
  fSparseOmegaDalitzMC(NULL),
  fSparseOmega2BodyMC(NULL),
  fSparsePhi2BodyMC(NULL),
  fSparseEtaPrimeDalitzMC(NULL),

  fHistTrackEta(NULL),
  fHistTrackThetaAbs(NULL),
  fHistTrackTriggerMatch(NULL),
  fHistTrackPDCA(NULL),
  fHistTrackChiSquare(NULL),
  fHistTriggerChiSquare(NULL)

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
  fSparseMixLSmmDimuon(NULL),

  fSparseEtaDalitz(NULL),
  fSparseEta2Body(NULL),
  fSparseRho2Body(NULL),
  fSparseOmegaDalitz(NULL),
  fSparseOmega2Body(NULL),
  fSparsePhi2Body(NULL),
  fSparseEtaPrimeDalitz(NULL),
  
  fSparseEtaDalitzMC(NULL),
  fSparseEta2BodyMC(NULL),
  fSparseRho2BodyMC(NULL),
  fSparseOmegaDalitzMC(NULL),
  fSparseOmega2BodyMC(NULL),
  fSparsePhi2BodyMC(NULL),
  fSparseEtaPrimeDalitzMC(NULL),

  fHistTrackEta(NULL),
  fHistTrackThetaAbs(NULL),
  fHistTrackTriggerMatch(NULL),
  fHistTrackPDCA(NULL),
  fHistTrackChiSquare(NULL),
  fHistTriggerChiSquare(NULL)
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
  
  if(fIsMC){
    
    const Int_t binnum_MC            = 5;
    Int_t bins_MC[binnum_MC]       = {800,6,300,12,300};
    double minbins_MC[binnum_MC] = {0,2.5,0,1.5,0};
    double maxbins_MC[binnum_MC] = {4,4.0,30,4.5,30};
    TString  namebins_MC[]        = {"rec_mass","rec_rap","rec_pt","true_rap","true_pt"};
    fSparseEtaDalitz= new THnSparseF("fSparseEtaDalitz","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);
    fSparseEta2Body= new THnSparseF("fSparseEta2Body","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);    
    fSparseOmegaDalitz= new THnSparseF("fSparseOmegaDalitz","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);
    fSparseOmega2Body= new THnSparseF("fSparseOmega2Body","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);
    fSparsePhi2Body= new THnSparseF("fSparsePhi2Body","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);
    fSparseEtaPrimeDalitz= new THnSparseF("fSparseEtaPrimeDalitz","",binnum_MC,bins_MC,minbins_MC,maxbins_MC);    
    fOutputList->Add(fSparseEtaDalitz);
    fOutputList->Add(fSparseEta2Body);
    fOutputList->Add(fSparseOmegaDalitz);
    fOutputList->Add(fSparseOmega2Body);
    fOutputList->Add(fSparsePhi2Body);
    fOutputList->Add(fSparseEtaPrimeDalitz);

    const Int_t binnum_True            = 2;
    Int_t bins_True[binnum_True]       = {12,300};
    double minbins_True[binnum_True] = {1.5,0};
    double maxbins_True[binnum_True] = {4.5,30};
    TString  namebins_True[]        = {"true_rap","true_pt"};
    fSparseEtaDalitzMC= new THnSparseF("fSparseEtaDalitzMC","",binnum_True,bins_True,minbins_True,maxbins_True);
    fSparseEta2BodyMC= new THnSparseF("fSparseEta2BodyMC","",binnum_True,bins_True,minbins_True,maxbins_True);
    fSparseOmegaDalitzMC= new THnSparseF("fSparseOmegaDalitzMC","",binnum_True,bins_True,minbins_True,maxbins_True);
    fSparseOmega2BodyMC= new THnSparseF("fSparseOmega2BodyMC","",binnum_True,bins_True,minbins_True,maxbins_True);
    fSparsePhi2BodyMC= new THnSparseF("fSparsePhi2BodyMC","",binnum_True,bins_True,minbins_True,maxbins_True);
    fSparseEtaPrimeDalitzMC= new THnSparseF("fSparseEtaPrimeDalitzMC","",binnum_True,bins_True,minbins_True,maxbins_True);    
    fOutputList->Add(fSparseEtaDalitzMC);
    fOutputList->Add(fSparseEta2BodyMC);    
    fOutputList->Add(fSparseOmegaDalitzMC);
    fOutputList->Add(fSparseOmega2BodyMC);
    fOutputList->Add(fSparsePhi2BodyMC);
    fOutputList->Add(fSparseEtaPrimeDalitzMC);
    
    const Int_t binnum_MC_Rho            = 6;
    Int_t bins_MC_Rho[binnum_MC_Rho]       = {800,6,300,12,300,800};
    double minbins_MC_Rho[binnum_MC_Rho] = {0,2.5,0,1.5,0,0};
    double maxbins_MC_Rho[binnum_MC_Rho] = {4,4.0,30,4.5,30,4};
    TString  namebins_MC_Rho[]        = {"rec_mass","rec_rap","rec_pt","true_rap","true_pt","true_mass"};
    fSparseRho2Body= new THnSparseF("fSparseRho2Body","",binnum_MC_Rho,bins_MC_Rho,minbins_MC_Rho,maxbins_MC_Rho);
    fOutputList->Add(fSparseRho2Body);

    const Int_t binnum_True_Rho            = 3;
    Int_t bins_True_Rho[binnum_True_Rho]       = {12,300,800};
    double minbins_True_Rho[binnum_True_Rho] = {1.5,0,0};
    double maxbins_True_Rho[binnum_True_Rho] = {4.5,30,4};
    TString  namebins_True_Rho[]        = {"true_rap","true_pt","true_mass"};
    fSparseRho2BodyMC= new THnSparseF("fSparseRho2BodyMC","",binnum_True_Rho,bins_True_Rho,minbins_True_Rho,maxbins_True_Rho);
    fOutputList->Add(fSparseRho2BodyMC);
  }

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

  PostData(1, fOutputList);    
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *)
{ 

  if(!Initialize()) return;
  if(!fUtils->isAcceptEvent()) return;

  MuonPairAnalysis();
  MuonPairAnalysisEveMixing();
  if(fIsMC) ProcessMC();
}

bool AliAnalysisTaskAODTrackPair::ProcessMC(){
  
  AliAODMCParticle *particle1;

  for(Int_t iTrack1=0; iTrack1<fMCTrackArray->GetSize(); ++iTrack1){
    
    particle1 = (AliAODMCParticle*)fMCTrackArray->At(iTrack1);

    if(!particle1) continue;
    if(!TDatabasePDG::Instance()->GetParticle(particle1->GetPdgCode())) continue;
    
    int pdgcode = particle1->GetPdgCode();
    
    double fill[]={fabs(particle1->Y()),particle1->Pt()};
    
    if(pdgcode == fUtils->fPdgCodeEta){
      if(fUtils->isDalitzProd()){
	fSparseEtaDalitzMC->Fill(fill);
      }
      if(fUtils->is2BodyProd()){
	fSparseEta2BodyMC->Fill(fill);
      }
    }
    else if(pdgcode == fUtils->fPdgCodeRho){
      double fill_rho[]={fabs(particle1->Y()),particle1->Pt(),particle1->GetCalcMass()};
      if(fUtils->is2BodyProd()){
	fSparseRho2BodyMC->Fill(fill_rho);
      }
    }
    else if(pdgcode == fUtils->fPdgCodeOmega){
      if(fUtils->isDalitzProd()){
	fSparseOmegaDalitzMC->Fill(fill);
      }
      if(fUtils->is2BodyProd()){
	fSparseOmega2BodyMC->Fill(fill);
      }	  
    }
    else if(pdgcode == fUtils->fPdgCodePhi){
      if(fUtils->is2BodyProd()){
	fSparsePhi2BodyMC->Fill(fill);
      }	  
    }
    else if(pdgcode == fUtils->fPdgCodeEtaPrime){
      if(fUtils->isDalitzProd()){
	fSparseEtaPrimeDalitzMC->Fill(fill);
      }
    }
  
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::Initialize()
{

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());

  if( !fUtils->setEvent(fEvent,fInputHandler) )
    return false;

  
  if( fRunNumber != fEvent->GetRunNumber() ){
    fRunNumber = fUtils->getRunnumber();
    AliMuonTrackCuts* trackCut = fUtils->getMuonTrackCuts();
    trackCut->SetRun(fInputHandler);
  }

  if(fIsMC){
    AliAODMCHeader* mcHeader = (AliAODMCHeader*) fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());  
    if(!mcHeader) return false;  
    
    TList* headerList = mcHeader->GetCocktailHeaders();
    if(!headerList) return false;
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)headerList->At(i) ;
      TString name = eventHeader2->GetName();

      if(name.Contains("EtaDalitz")){
	fUtils->setDalitzProd(true);
      }
      else if(name.Contains("EtaDirect")){
	fUtils->set2BodyProd(true);
      }
      else if(name.Contains("OmegaDalitz")){
	fUtils->setDalitzProd(true);
      }
      else if(name.Contains("OmegaDirect")){
	fUtils->set2BodyProd(true);
      }
      else if(name.Contains("EtaPrimeDalitz")){
	fUtils->setDalitzProd(true);
      }
      else if(name.Contains("PhiDirect")){
	fUtils->set2BodyProd(true);
      }
      else if(name.Contains("RhoDirect")){
	fUtils->set2BodyProd(true);
      }
      else{
	fUtils->setDalitzProd(false);
	fUtils->set2BodyProd(false);
      }      
    }
    
    Double_t impPar = mcHeader->GetImpactParameter();    
    
    fMCTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCTrackArray) return false;
    fUtils->setMCArray(fMCTrackArray);
  }

  fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);
  
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::MuonTrackQA(AliAODTrack* track){
  
  fHistTrackEta->Fill(track->Pt(),track->Eta());
  fHistTrackThetaAbs->Fill(track->Pt(),AliAnalysisMuonUtility::GetThetaAbsDeg(track));
  fHistTrackTriggerMatch->Fill(track->Pt(),AliAnalysisMuonUtility::GetMatchTrigger(track));
  //fHistTrackPDCA->Fill();
  fHistTrackChiSquare->Fill(track->Pt(),AliAnalysisMuonUtility::GetChi2perNDFtracker(track));
  fHistTriggerChiSquare->Fill(track->Pt(),AliAnalysisMuonUtility::GetChi2MatchTrigger(track));
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::MuonPairAnalysisEveMixing(){
  
  if( !(fInputHandler->IsEventSelected() & fTriggerMaskForMixing) ) return false;
  
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

	  double fill[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fUtils->getCentClass()};
	  
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
  
  if(!fIsMC && !(fInputHandler->IsEventSelected() & fTriggerMaskForSame) ) return false;

  if(fIsCMUL7){
    fEventCounter->Fill(0.,fUtils->getCentClass());
    fEventCounter->Fill(4.,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMLL7){
    fEventCounter->Fill(1.,fUtils->getCentClass());
    fEventCounter->Fill(5.,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 || fIsCMLL7){
    fEventCounter->Fill(2.,fUtils->getCentClass());
    fEventCounter->Fill(6.,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }
  if(fIsCMUL7 && fIsCMLL7){
    fEventCounter->Fill(3.,fUtils->getCentClass());
    fEventCounter->Fill(7.,fUtils->getCentClass(),(double)1./fUtils->getDS());
  }

  Int_t nTrack = fEvent->GetNumberOfTracks();
  
  AliAODMCParticle* particle1;
  AliAODTrack* track1;
  AliAODTrack* track2;
  AliAODDimuon* dimuon;

  for(Int_t iTrack1=0; iTrack1<nTrack; ++iTrack1){
   
    track1 = (AliAODTrack*)fEvent->GetTrack(iTrack1);
    
    if(!fUtils->isAcceptMuonTrack(track1)) continue;

    MuonTrackQA(track1);

    for(Int_t iTrack2=iTrack1+1; iTrack2<nTrack; ++iTrack2){
   
      track2 = (AliAODTrack*)fEvent->GetTrack(iTrack2);
      
      if(!fUtils->isAcceptMuonTrack(track2)) continue;

      dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1,track2);
      
      if(!fUtils->isAcceptDimuon(dimuon)) continue;
     
      double fill[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fUtils->getCentClass()};
      
      string fFiredTrigName = string(fEvent->GetFiredTriggerClasses());
      
      if(dimuon->Charge() == 0){
	fSparseULSDimuon->Fill(fill,(double)1./fUtils->getDS());
      }
      else if(dimuon->Charge() > 0){
	fSparseLSppDimuon->Fill(fill,(double)1./fUtils->getDS());
      }
      else{
	fSparseLSmmDimuon->Fill(fill,(double)1./fUtils->getDS());
      }
      
      if(!fIsMC) continue;

      if(fUtils->isSameMotherPair(track1,track2)){
	int mom_pdg=fUtils->getMotherPdgCode(track1);
	int mom_label=fUtils->getMotherLabel(track1);
	
	particle1 = (AliAODMCParticle*)fMCTrackArray->At(mom_label);

	double fill_mc[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fabs(particle1->Y()),particle1->Pt()};

	if(mom_pdg == fUtils->fPdgCodeEta){
	  if(fUtils->isDalitzProd()){
	    fSparseEtaDalitz->Fill(fill_mc);
	  }
	  if(fUtils->is2BodyProd()){
	    fSparseEta2Body->Fill(fill_mc);
	  }
	}
	else if(mom_pdg == fUtils->fPdgCodeRho){
	  double fill_mc_rho[]={dimuon->M(),fabs(dimuon->Y()),dimuon->Pt(),fabs(particle1->Y()),particle1->Pt(),particle1->GetCalcMass()};
	  if(fUtils->is2BodyProd()){
	    fSparseRho2Body->Fill(fill_mc_rho);
	  }
	}
	else if(mom_pdg == fUtils->fPdgCodeOmega){
	  if(fUtils->isDalitzProd()){
	    fSparseOmegaDalitz->Fill(fill_mc);
	  }
	  if(fUtils->is2BodyProd()){
	    fSparseOmega2Body->Fill(fill_mc);
	  }	  
	}
	else if(mom_pdg == fUtils->fPdgCodePhi){
	  if(fUtils->is2BodyProd()){
	    fSparsePhi2Body->Fill(fill_mc);
	  }	  
	}
	else if(mom_pdg == fUtils->fPdgCodeEtaPrime){
	  if(fUtils->isDalitzProd()){
	    fSparseEtaPrimeDalitz->Fill(fill_mc);
	  }
	}

      }


      delete dimuon;

    }//end of loop track2
  }//end of loop track1 
}
