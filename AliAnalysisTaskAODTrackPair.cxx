#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRefArray.h"
#include "TTree.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEventPoolManager.h"
#include "AliEventPoolMuon.h"
#include "AliGenEventHeader.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliAODRecoDecay.h"
#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODv0.h"
#include "AliExternalTrackParam.h"
#include "AliPIDCombined.h"
#include "AliTOFHeader.h"

#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskAODTrackPair.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)

    AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair()
    : AliAnalysisTaskSE(), fEvent(NULL), fPoolMuonTrackMgr(NULL), fUtils(NULL),
  fMCTrackArray(NULL), fIsMC(false), fIsMidTrackAna(false), fIsV0TrackPairAna(false),
      fIsPrimTrackPairAna(false), fIsMixingAnalysis(false), fRunNumber(-99999),
      fTrackDepth(1), fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackPairPtBalance(NULL), fHistTrackPairLocalBoardPair(NULL),

      fHistTrackThetaAbs(NULL), fHistTrackTriggerMatch(NULL),
      fHistTrackPDCA(NULL), fHistTrackChiSquare(NULL),
      fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeULSPair(NULL), fTreeLSppPair(NULL), fTreeLSmmPair(NULL),

      fSparseULSPairMassPt(NULL), fSparseLSppPairMassPt(NULL),
      fSparseLSmmPairMassPt(NULL),
  
  fSparseULSPairMassPt_LeadingTrack(NULL),
  fSparseULSPairMassPt_RejectKpmStar(NULL),
  
  fSparseULSPairMassPt_SideBandLeftRight(NULL),
  fSparseULSPairMassPt_SideBandLeft(NULL),
  fSparseULSPairMassPt_SideBandRight(NULL),
  fSparseULSPairMassPt_SideBand(NULL),

  fHistTrueK0sPtRapidity(NULL),
  fHistRecTrueK0sPtRapidity(NULL),

      fTreeULSPair_TightCut(NULL), fTreeLSppPair_TightCut(NULL),
      fTreeLSmmPair_TightCut(NULL),

      fTreeULSPair_ProngV0(NULL), fTreeLSppPair_ProngV0(NULL),
      fTreeLSmmPair_ProngV0(NULL),

      fTreeMixULSPair(NULL), fTreeMixLSppPair(NULL), fTreeMixLSmmPair(NULL),
    
      fHistULSPairMassPt(NULL), fHistLSppPairMassPt(NULL),
      fHistLSmmPairMassPt(NULL),

      fHistULSPairMassPt_TightCut(NULL), fHistLSppPairMassPt_TightCut(NULL),
      fHistLSmmPairMassPt_TightCut(NULL),

  fHistULSPairMassPt_ProngV0(NULL),
  fHistLSppPairMassPt_ProngV0(NULL),
      fHistLSmmPairMassPt_ProngV0(NULL),

  fHistULSPairMassPt_ProngV0_TrueK0s(NULL),

      fHistMixULSPairMassPt(NULL), fHistMixLSppPairMassPt(NULL),
      fHistMixLSmmPairMassPt(NULL),

      fSparseMixULSPairMassPt(NULL), fSparseMixLSppPairMassPt(NULL),
      fSparseMixLSmmPairMassPt(NULL),

  fSparseMixULSPairMassPt_LeadingTrack(NULL),
  fSparseMixULSPairMassPt_RejectKpmStar(NULL),

      fHistMassK0s1K0s2(NULL),
  fHistKpmStarCandMass(NULL),
  fHistPionPionCorrelationPlot(NULL),

      fHistTPCdEdxP(NULL), fHistBetaP(NULL), fHistTPCSigmaElectron(NULL),
      fHistTOFSigmaElectron(NULL), fHistTPCSigmaMuon(NULL),
      fHistTOFSigmaMuon(NULL), fHistTPCSigmaPion(NULL), fHistTOFSigmaPion(NULL),
      fHistTPCSigmaKaon(NULL), fHistTOFSigmaKaon(NULL),
      fHistTPCSigmaProton(NULL), fHistTOFSigmaProton(NULL),

      fHistSelTPCdEdxP(NULL), fHistSelBetaP(NULL),
      fHistSelTPCSigmaElectron(NULL), fHistSelTOFSigmaElectron(NULL),
      fHistSelTPCSigmaMuon(NULL), fHistSelTOFSigmaMuon(NULL),
      fHistSelTPCSigmaPion(NULL), fHistSelTOFSigmaPion(NULL),
      fHistSelTPCSigmaKaon(NULL), fHistSelTOFSigmaKaon(NULL),
      fHistSelTPCSigmaProton(NULL), fHistSelTOFSigmaProton(NULL),

      fHistTrackP(NULL), fHistTrackPt(NULL), fHistTrackEta(NULL),
      fHistTPCNClusts(NULL), fHistSPDNClusts(NULL),
      fHistTPCCrossRowsFindableRatio(NULL), fHistReducedChi2TPC(NULL),
  fHistReducedChi2ITS(NULL), fHistDCAz(NULL), fHistDCAxyPt(NULL), fHistOpeningAngle(NULL),
  
  fHistArmenteros(NULL), fHistV0MassDecayLength(NULL),
  fHistV0MassPointingAngle(NULL), fHistV0MassV0DCA(NULL),
  fHistV0MassV0TrackDCA(NULL), fHistV0MassV0DecayRadius(NULL),
  fHistV0MassV0PropLifeTime(NULL), 
  
  fHistSelArmenteros(NULL), fHistSelV0MassDecayLength(NULL),
  fHistSelV0MassPointingAngle(NULL), fHistSelV0MassV0DCA(NULL),
  fHistSelV0MassV0TrackDCA(NULL), fHistSelV0MassV0DecayRadius(NULL),
  fHistSelV0MassV0PropLifeTime(NULL),

  fHistArmenteros_TrueK0s(NULL), fHistV0MassDecayLength_TrueK0s(NULL),
  fHistV0MassPointingAngle_TrueK0s(NULL), fHistV0MassV0DCA_TrueK0s(NULL),
  fHistV0MassV0TrackDCA_TrueK0s(NULL), fHistV0MassV0DecayRadius_TrueK0s(NULL),
  fHistV0MassV0PropLifeTime_TrueK0s(NULL), 
  
  fHistSelArmenteros_TrueK0s(NULL), fHistSelV0MassDecayLength_TrueK0s(NULL),
  fHistSelV0MassPointingAngle_TrueK0s(NULL), fHistSelV0MassV0DCA_TrueK0s(NULL),
  fHistSelV0MassV0TrackDCA_TrueK0s(NULL), fHistSelV0MassV0DecayRadius_TrueK0s(NULL),
  fHistSelV0MassV0PropLifeTime_TrueK0s(NULL),
  
  fTreeK0s(NULL),
  TrueMass(0.),
  TrueRap(0.),
  TruePhi(0.),
  TruePt(0.),
  RecMass(0.),
  RecRap(0.),
  RecPhi(0.),
  RecPt(0.),
  isDetect(false),
  
      RecPairPt(0.), RecPairRap(0.), RecPairMass(0.),
      RecPairArmenterosArmPt(0.), RecPairArmenterosAlpha(0.), RecPairCent(0.),
      RecPairDS(0.) {}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char *name)
    : AliAnalysisTaskSE(name), fEvent(NULL), fPoolMuonTrackMgr(NULL),
      fUtils(NULL), fMCTrackArray(NULL), fIsMC(false), fIsMidTrackAna(false),
      fIsV0TrackPairAna(false), fIsPrimTrackPairAna(false),
      fIsMixingAnalysis(false), fRunNumber(-99999), fTrackDepth(1),
      fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackPairPtBalance(NULL), fHistTrackPairLocalBoardPair(NULL),

      fHistTrackThetaAbs(NULL), fHistTrackTriggerMatch(NULL),
      fHistTrackPDCA(NULL), fHistTrackChiSquare(NULL),
      fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeULSPair(NULL), fTreeLSppPair(NULL), fTreeLSmmPair(NULL),

      fSparseULSPairMassPt(NULL), fSparseLSppPairMassPt(NULL),
      fSparseLSmmPairMassPt(NULL),
  

  fSparseULSPairMassPt_LeadingTrack(NULL),
  fSparseULSPairMassPt_RejectKpmStar(NULL),


  fSparseULSPairMassPt_SideBandLeftRight(NULL),
  fSparseULSPairMassPt_SideBandLeft(NULL),
  fSparseULSPairMassPt_SideBandRight(NULL),
  fSparseULSPairMassPt_SideBand(NULL),

  fHistTrueK0sPtRapidity(NULL),
  fHistRecTrueK0sPtRapidity(NULL),

      fTreeULSPair_TightCut(NULL), fTreeLSppPair_TightCut(NULL),
      fTreeLSmmPair_TightCut(NULL),

      fTreeULSPair_ProngV0(NULL), fTreeLSppPair_ProngV0(NULL),
      fTreeLSmmPair_ProngV0(NULL),

      fTreeMixULSPair(NULL), fTreeMixLSppPair(NULL), fTreeMixLSmmPair(NULL),

      fHistULSPairMassPt(NULL), fHistLSppPairMassPt(NULL),
      fHistLSmmPairMassPt(NULL),

      fHistULSPairMassPt_TightCut(NULL), fHistLSppPairMassPt_TightCut(NULL),
      fHistLSmmPairMassPt_TightCut(NULL),

      fHistULSPairMassPt_ProngV0(NULL), fHistLSppPairMassPt_ProngV0(NULL),
      fHistLSmmPairMassPt_ProngV0(NULL),

  fHistULSPairMassPt_ProngV0_TrueK0s(NULL),

      fHistMixULSPairMassPt(NULL), fHistMixLSppPairMassPt(NULL),
      fHistMixLSmmPairMassPt(NULL),

      fSparseMixULSPairMassPt(NULL), fSparseMixLSppPairMassPt(NULL),
      fSparseMixLSmmPairMassPt(NULL),  

  fSparseMixULSPairMassPt_LeadingTrack(NULL),
  fSparseMixULSPairMassPt_RejectKpmStar(NULL),

      fHistMassK0s1K0s2(NULL),
  fHistKpmStarCandMass(NULL),
  fHistPionPionCorrelationPlot(NULL),

      fHistTPCdEdxP(NULL), fHistBetaP(NULL), fHistTPCSigmaElectron(NULL),
      fHistTOFSigmaElectron(NULL), fHistTPCSigmaMuon(NULL),
      fHistTOFSigmaMuon(NULL), fHistTPCSigmaPion(NULL), fHistTOFSigmaPion(NULL),
      fHistTPCSigmaKaon(NULL), fHistTOFSigmaKaon(NULL),
      fHistTPCSigmaProton(NULL), fHistTOFSigmaProton(NULL),

      fHistSelTPCdEdxP(NULL), fHistSelBetaP(NULL),
      fHistSelTPCSigmaElectron(NULL), fHistSelTOFSigmaElectron(NULL),
      fHistSelTPCSigmaMuon(NULL), fHistSelTOFSigmaMuon(NULL),
      fHistSelTPCSigmaPion(NULL), fHistSelTOFSigmaPion(NULL),
      fHistSelTPCSigmaKaon(NULL), fHistSelTOFSigmaKaon(NULL),
      fHistSelTPCSigmaProton(NULL), fHistSelTOFSigmaProton(NULL),

      fHistTrackP(NULL), fHistTrackPt(NULL), fHistTrackEta(NULL),
      fHistTPCNClusts(NULL), fHistSPDNClusts(NULL),
      fHistTPCCrossRowsFindableRatio(NULL), fHistReducedChi2TPC(NULL),
  fHistReducedChi2ITS(NULL), fHistDCAz(NULL), fHistDCAxyPt(NULL), fHistOpeningAngle(NULL),

  fHistArmenteros(NULL), fHistV0MassDecayLength(NULL),
  fHistV0MassPointingAngle(NULL), fHistV0MassV0DCA(NULL),
  fHistV0MassV0TrackDCA(NULL), fHistV0MassV0DecayRadius(NULL),
  fHistV0MassV0PropLifeTime(NULL), 
  
  fHistSelArmenteros(NULL), fHistSelV0MassDecayLength(NULL),
  fHistSelV0MassPointingAngle(NULL), fHistSelV0MassV0DCA(NULL),
  fHistSelV0MassV0TrackDCA(NULL), fHistSelV0MassV0DecayRadius(NULL),
  fHistSelV0MassV0PropLifeTime(NULL),

  fHistArmenteros_TrueK0s(NULL), fHistV0MassDecayLength_TrueK0s(NULL),
  fHistV0MassPointingAngle_TrueK0s(NULL), fHistV0MassV0DCA_TrueK0s(NULL),
  fHistV0MassV0TrackDCA_TrueK0s(NULL), fHistV0MassV0DecayRadius_TrueK0s(NULL),
  fHistV0MassV0PropLifeTime_TrueK0s(NULL), 
  
  fHistSelArmenteros_TrueK0s(NULL), fHistSelV0MassDecayLength_TrueK0s(NULL),
  fHistSelV0MassPointingAngle_TrueK0s(NULL), fHistSelV0MassV0DCA_TrueK0s(NULL),
  fHistSelV0MassV0TrackDCA_TrueK0s(NULL), fHistSelV0MassV0DecayRadius_TrueK0s(NULL),
  fHistSelV0MassV0PropLifeTime_TrueK0s(NULL),

  fTreeK0s(NULL),
  TrueMass(0.),
  TrueRap(0.),
  TruePhi(0.),
  TruePt(0.),
  RecMass(0.),
  RecRap(0.),
  RecPhi(0.),
  RecPt(0.),
  isDetect(false),
      RecPairPt(0.), RecPairRap(0.), RecPairMass(0.),
      RecPairArmenterosArmPt(0.), RecPairArmenterosAlpha(0.), RecPairCent(0.),
      RecPairDS(0.) {

  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPair::~AliAnalysisTaskAODTrackPair() {}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPair::UserCreateOutputObjects() {

  double fCentBins[] = {0,1,5,10,20,40,60,100};
  double fVtxBins[] = {-50, -10.5, -6, -2, 0, 2, 6, 10.5, 50};
  double fPsiBins[] = {-10, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 10};
  
  int fNCentBins = sizeof(fCentBins) / sizeof(double) - 1;
  int fNVtxZBins = sizeof(fVtxBins) / sizeof(double) - 1;
  int fNPsiBins = sizeof(fPsiBins) / sizeof(double) - 1;

  fPoolMuonTrackMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, fNCentBins, (double *)fCentBins, fNVtxZBins,
					      (double *)fVtxBins, fNPsiBins, (double *)fPsiBins);
  fPoolMuonTrackMgr->SetTargetValues(fTrackDepth, (double)fReadyFraction,fPoolSize);
  
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  double min_mass = 0.0;
  double max_mass = 6.0;
  double width_mass = 0.001;

  double min_pt = 0.0;
  double max_pt = 20.0;
  double width_pt = 0.5;

  double min_rap = -0.8;
  double max_rap = 0.8;
  double width_rap = 0.1;

  vector<double> bins_cent{0,1,5,10,20,40,60,80,100};
  int binnum_cent = bins_cent.size() - 1;

  vector<double> bins_event{0,60,120,180};
  int binnum_event = bins_event.size() - 1;

  if (!fIsMidTrackAna) {
    double bins_event_hist[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int binnum_event_hist = sizeof(bins_event_hist) / sizeof(double) - 1;

    std::string event_label[] = {"CMUL7",
				 "CMLL7",
				 "CMUL7orCMLL7",
				 "CMUL7andCMLL7",
				 "CMUL7withDS",
				 "CMLL7withDS",
				 "CMUL7orCMLL7withDS",
				 "CMUL7andCMLL7withDS"};
    fEventCounter = new TH2F("fEventCounter", "", 11, 0, 11, 200, 0, 200);
    for (unsigned int iname = 0; iname < sizeof(event_label) / sizeof(std::string); ++iname) {
      fEventCounter->GetXaxis()->SetBinLabel(iname + 1,event_label[iname].c_str());
    }

    fOutputList->Add(fEventCounter);

    if (!fIsMixingAnalysis) {
      fTreeULSPair = new TTree("fTreeULSPair", "");
      fTreeULSPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeULSPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeULSPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeULSPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeULSPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeULSPair);
      fTreeLSppPair = new TTree("fTreeLSppPair", "");
      fTreeLSppPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeLSppPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeLSppPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeLSppPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeLSppPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeLSppPair);
      fTreeLSmmPair = new TTree("fTreeLSmmPair", "");
      fTreeLSmmPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeLSmmPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeLSmmPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeLSmmPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeLSmmPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeLSmmPair);
    } else {
      fTreeMixULSPair = new TTree("fTreeMixULSPair", "");
      fTreeMixULSPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixULSPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixULSPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixULSPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixULSPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixULSPair);
      fTreeMixLSppPair = new TTree("fTreeMixLSppPair", "");
      fTreeMixLSppPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixLSppPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixLSppPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixLSppPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixLSppPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixLSppPair);
      fTreeMixLSmmPair = new TTree("fTreeMixLSmmPair", "");
      fTreeMixLSmmPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixLSmmPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixLSmmPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixLSmmPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixLSmmPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixLSmmPair);
    }
    fHistTrackPairPtBalance =
        new TH2F("fHistTrackPairPtBalance", "", 50, 0, 5, 50, 0, 5);
    fHistTrackPairLocalBoardPair =
        new TH2F("fHistTrackPairLocalBoardPair", "", 240, 0, 240, 240, 0, 240);
    fOutputList->Add(fHistTrackPairPtBalance);
    fOutputList->Add(fHistTrackPairLocalBoardPair);
    
    fHistTrackEta = new TH1F("fHistTrackEta", "", 30, 2.0,5.0);    
    fHistTrackThetaAbs =
        new TH2F("fHistTrackThetaAbs", "", 20, 0, 10, 60, 0, 15);
    fHistTrackTriggerMatch =
        new TH2F("fHistTrackTriggerMatch", "", 20, 0, 10, 5, 0, 5);
    fHistTrackPDCA = new TH2F("fHistTrackPDCA", "", 20, 0, 10, 200, 0, 20);
    fHistTrackChiSquare =
        new TH2F("fHistTrackChiSquare", "", 20, 0, 10, 100, 0, 10);
    fHistTriggerChiSquare =
        new TH2F("fHistTriggerChiSquare", "", 20, 0, 10, 100, 0, 10);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTrackThetaAbs);
    fOutputList->Add(fHistTrackTriggerMatch);
    fOutputList->Add(fHistTrackPDCA);
    fOutputList->Add(fHistTrackChiSquare);
    fOutputList->Add(fHistTriggerChiSquare);
  } else {
    if (!fIsMixingAnalysis) {

      int bins[3] = {int((max_mass - min_mass) / width_mass),
                     int((max_pt - min_pt) / width_pt), binnum_cent};
      double min_bins[3] = {min_mass, min_pt, 0};
      double max_bins[3] = {max_mass, max_pt, 100};      
      
      fSparseULSPairMassPt = new THnSparseF("fSparseULSPairMassPt", "", 3, bins,
                                            min_bins, max_bins);
      fSparseLSppPairMassPt = new THnSparseF("fSparseLSppPairMassPt", "", 3,
                                             bins, min_bins, max_bins);
      fSparseLSmmPairMassPt = new THnSparseF("fSparseLSmmPairMassPt", "", 3,
                                             bins, min_bins, max_bins);
      fSparseULSPairMassPt->SetBinEdges(2,bins_cent.data());
      fSparseLSppPairMassPt->SetBinEdges(2,bins_cent.data());
      fSparseLSmmPairMassPt->SetBinEdges(2,bins_cent.data());
      fOutputList->Add(fSparseULSPairMassPt);
      fOutputList->Add(fSparseLSppPairMassPt);
      fOutputList->Add(fSparseLSmmPairMassPt);

      fSparseULSPairMassPt_RejectKpmStar = new THnSparseF("fSparseULSPairMassPt_RejectKpmStar", "", 3, bins,
							  min_bins, max_bins);
      fSparseULSPairMassPt_RejectKpmStar->SetBinEdges(2,bins_cent.data());
      fOutputList->Add(fSparseULSPairMassPt_RejectKpmStar);

      fSparseULSPairMassPt_SideBandLeftRight = new THnSparseF("fSparseULSPairMassPt_SideBandLeftRight", "", 3, bins,
							 min_bins, max_bins);
      fSparseULSPairMassPt_SideBandLeft = new THnSparseF("fSparseULSPairMassPt_SideBandLeft", "", 3, bins,
							 min_bins, max_bins);
      fSparseULSPairMassPt_SideBandRight = new THnSparseF("fSparseULSPairMassPt_SideBandRight", "", 3, bins,
							 min_bins, max_bins);
      fSparseULSPairMassPt_SideBand = new THnSparseF("fSparseULSPairMassPt_SideBand", "", 3, bins,
						     min_bins, max_bins);
      fSparseULSPairMassPt_SideBandLeftRight->SetBinEdges(2,bins_cent.data());
      fSparseULSPairMassPt_SideBandLeft->SetBinEdges(2,bins_cent.data());
      fSparseULSPairMassPt_SideBandRight->SetBinEdges(2,bins_cent.data());
      fSparseULSPairMassPt_SideBand->SetBinEdges(2,bins_cent.data());
      fOutputList->Add(fSparseULSPairMassPt_SideBandLeftRight);
      fOutputList->Add(fSparseULSPairMassPt_SideBandLeft);
      fOutputList->Add(fSparseULSPairMassPt_SideBandRight);
      fOutputList->Add(fSparseULSPairMassPt_SideBand);

      int bins_leading[3] = {int((max_mass - min_mass) / width_mass),
			     int((max_pt - min_pt) / width_pt), binnum_event};
      double min_bins_leading[3] = {min_mass, min_pt, 0};
      double max_bins_leading[3] = {max_mass, max_pt, 360};      

      fSparseULSPairMassPt_LeadingTrack  = new THnSparseF("fSparseULSPairMassPt_LeadingTrack", "",
									    3, bins_leading,min_bins_leading, max_bins_leading);
      fSparseULSPairMassPt_LeadingTrack->SetBinEdges(2,bins_event.data());
      fOutputList->Add(fSparseULSPairMassPt_LeadingTrack);
    }
    if (fIsMixingAnalysis) {
      fHistMixULSPairMassPt =
          new TH2F("fHistMixULSPairMassPt", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistMixULSPairMassPt);
      
      int bins[3] = {int((max_mass - min_mass) / width_mass),
                     int((max_pt - min_pt) / width_pt), binnum_cent};
      double min_bins[3] = {min_mass, min_pt, 0};
      double max_bins[3] = {max_mass, max_pt, 100};
      fSparseMixULSPairMassPt = new THnSparseF("fSparseMixULSPairMassPt", "", 3,
                                               bins, min_bins, max_bins);
      fSparseMixLSppPairMassPt = new THnSparseF("fSparseMixLSppPairMassPt", "",
                                                3, bins, min_bins, max_bins);
      fSparseMixLSmmPairMassPt = new THnSparseF("fSparseMixLSmmPairMassPt", "",
                                                3, bins, min_bins, max_bins);      
      fSparseMixULSPairMassPt->SetBinEdges(2,bins_cent.data());
      fSparseMixLSppPairMassPt->SetBinEdges(2,bins_cent.data());
      fSparseMixLSmmPairMassPt->SetBinEdges(2,bins_cent.data());      
      fOutputList->Add(fSparseMixULSPairMassPt);
      fOutputList->Add(fSparseMixLSppPairMassPt);
      fOutputList->Add(fSparseMixLSmmPairMassPt);

      fSparseMixULSPairMassPt_RejectKpmStar = new THnSparseF("fSparseMixULSPairMassPt_RejectKpmStar", "", 3, bins,
							     min_bins, max_bins);
      fSparseMixULSPairMassPt_RejectKpmStar->SetBinEdges(2,bins_cent.data());
      fOutputList->Add(fSparseMixULSPairMassPt_RejectKpmStar);

      int bins_leading[3] = {int((max_mass - min_mass) / width_mass),
			     int((max_pt - min_pt) / width_pt), binnum_event};
      double min_bins_leading[3] = {min_mass, min_pt, 0};
      double max_bins_leading[3] = {max_mass, max_pt, 360};      

      fSparseMixULSPairMassPt_LeadingTrack  = new THnSparseF("fSparseMixULSPairMassPt_LeadingTrack", "",
									    3, bins_leading,min_bins_leading, max_bins_leading);
      fSparseMixULSPairMassPt_LeadingTrack->SetBinEdges(2,bins_event.data());
      fOutputList->Add(fSparseMixULSPairMassPt_LeadingTrack);
    }

    if (fIsV0TrackPairAna) {

      fHistPionPionCorrelationPlot = new TH2F("fHistPionPionCorrelationPlot", "",
				 int((max_mass - min_mass)/width_mass),min_mass,max_mass,
				 int((max_mass - min_mass)/width_mass),min_mass,max_mass);
      fOutputList->Add(fHistPionPionCorrelationPlot);

      min_mass = 0.45;
      max_mass = 0.55;
      width_mass = 0.001;

      fHistULSPairMassPt_ProngV0 =
          new TH2F("fHistULSPairMassPt_ProngV0", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistULSPairMassPt_ProngV0);

      fHistMassK0s1K0s2 =
	new TH2F("fHistMassK0s1K0s2", "",
		 int((max_mass - min_mass)/width_mass),min_mass,max_mass,
		 int((max_mass - min_mass)/width_mass),min_mass,max_mass);
      fOutputList->Add(fHistMassK0s1K0s2);

      fHistKpmStarCandMass = new TH1F("fHistKpsStarCandMass","",400,0.7,1.1);
      fOutputList->Add(fHistKpmStarCandMass);

      fHistArmenteros =
          new TH2F("fHistArmenteros", "", 200, -1, 1, 400, 0, 0.4);
      fHistV0MassDecayLength = new TH2F("fHistV0MassDecayLength", "",
                                        int((max_mass - min_mass) / width_mass),
                                        min_mass, max_mass, 500, 0., 100.);
      fHistV0MassPointingAngle =
          new TH2F("fHistV0MassPointingAngle", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0.95, 1.0);
      fHistV0MassV0DCA = new TH2F("fHistV0MassV0DCA", "",
                                  int((max_mass - min_mass) / width_mass),
                                  min_mass, max_mass, 400, 0., 80.);
      fHistV0MassV0TrackDCA = new TH2F("fHistV0MassV0TrackDCA", "",
                                       int((max_mass - min_mass) / width_mass),
                                       min_mass, max_mass, 100, 0, 10);
      fHistV0MassV0DecayRadius =
          new TH2F("fHistV0MassV0DecayRadius", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
      fHistV0MassV0PropLifeTime =
          new TH2F("fHistV0MassV0PropLifeTime", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   250, 0, 50.);
      fOutputList->Add(fHistArmenteros);     
      fOutputList->Add(fHistV0MassDecayLength);
      fOutputList->Add(fHistV0MassPointingAngle);
      fOutputList->Add(fHistV0MassV0DCA);
      fOutputList->Add(fHistV0MassV0TrackDCA);
      fOutputList->Add(fHistV0MassV0DecayRadius);
      fOutputList->Add(fHistV0MassV0PropLifeTime);

      fHistSelArmenteros =
          new TH2F("fHistSelArmenteros", "", 200, -1, 1, 400, 0, 0.4);
      fHistSelV0MassDecayLength =
          new TH2F("fHistSelV0MassDecayLength", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
      fHistSelV0MassPointingAngle =
          new TH2F("fHistSelV0MassPointingAngle", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0.95, 1.0);
      fHistSelV0MassV0DCA = new TH2F("fHistSelV0MassV0DCA", "",
                                     int((max_mass - min_mass) / width_mass),
                                     min_mass, max_mass, 400, 0., 80.);
      fHistSelV0MassV0TrackDCA =
          new TH2F("fHistSelV0MassV0TrackDCA", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   100, 0, 10);
      fHistSelV0MassV0DecayRadius =
          new TH2F("fHistSelV0MassV0DecayRadius", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
      fHistSelV0MassV0PropLifeTime =
          new TH2F("fHistSelV0MassV0PropLifeTime", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   250, 0, 50.);
      fOutputList->Add(fHistSelArmenteros);
      fOutputList->Add(fHistSelV0MassDecayLength);
      fOutputList->Add(fHistSelV0MassPointingAngle);
      fOutputList->Add(fHistSelV0MassV0DCA);
      fOutputList->Add(fHistSelV0MassV0TrackDCA);
      fOutputList->Add(fHistSelV0MassV0DecayRadius);
      fOutputList->Add(fHistSelV0MassV0PropLifeTime);

      if (fIsMC) {

	fHistTrueK0sPtRapidity = new TH2F("fHistTrueK0sPtRapidity", "",
					     int((max_rap - min_rap) / width_rap), min_rap, max_rap,
					     int((max_pt - min_pt) / width_pt), min_pt, max_pt);	
	fOutputList->Add(fHistTrueK0sPtRapidity);      

	fHistRecTrueK0sPtRapidity = new TH2F("fHistRecTrueK0sPtRapidity", "",
					     int((max_rap - min_rap) / width_rap), min_rap, max_rap,
					     int((max_pt - min_pt) / width_pt), min_pt, max_pt);	
	fOutputList->Add(fHistRecTrueK0sPtRapidity);      


	fHistULSPairMassPt_ProngV0_TrueK0s =
          new TH2F("fHistULSPairMassPt_ProngV0_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
	fOutputList->Add(fHistULSPairMassPt_ProngV0_TrueK0s);      
	
	fHistArmenteros_TrueK0s =
          new TH2F("fHistArmenteros_TrueK0s", "", 200, -1, 1, 400, 0, 0.4);
	fHistV0MassDecayLength_TrueK0s =
          new TH2F("fHistV0MassDecayLength_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
	fHistV0MassPointingAngle_TrueK0s =
          new TH2F("fHistV0MassPointingAngle_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0.95, 1.0);
	fHistV0MassV0DCA_TrueK0s = new TH2F("fHistV0MassV0DCA_TrueK0s", "",
					    int((max_mass - min_mass) / width_mass),
					    min_mass, max_mass, 400, 0., 80.);
	fHistV0MassV0TrackDCA_TrueK0s =
          new TH2F("fHistV0MassV0TrackDCA_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   100, 0, 10);
	fHistV0MassV0DecayRadius_TrueK0s =
          new TH2F("fHistV0MassV0DecayRadius_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
	fHistV0MassV0PropLifeTime_TrueK0s =
          new TH2F("fHistV0MassV0PropLifeTime_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   250, 0, 50.);
	fOutputList->Add(fHistArmenteros_TrueK0s);
	fOutputList->Add(fHistV0MassDecayLength_TrueK0s);
	fOutputList->Add(fHistV0MassPointingAngle_TrueK0s);
	fOutputList->Add(fHistV0MassV0DCA_TrueK0s);
	fOutputList->Add(fHistV0MassV0TrackDCA_TrueK0s);
	fOutputList->Add(fHistV0MassV0DecayRadius_TrueK0s);
	fOutputList->Add(fHistV0MassV0PropLifeTime_TrueK0s);

	fHistSelArmenteros_TrueK0s =
          new TH2F("fHistSelArmenteros_TrueK0s", "", 200, -1, 1, 400, 0, 0.4);
	fHistSelV0MassDecayLength_TrueK0s =
          new TH2F("fHistSelV0MassDecayLength_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
	fHistSelV0MassPointingAngle_TrueK0s =
          new TH2F("fHistSelV0MassPointingAngle_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0.95, 1.0);
	fHistSelV0MassV0DCA_TrueK0s = new TH2F("fHistSelV0MassV0DCA_TrueK0s", "",
					    int((max_mass - min_mass) / width_mass),
					    min_mass, max_mass, 400, 0., 80.);
	fHistSelV0MassV0TrackDCA_TrueK0s =
          new TH2F("fHistSelV0MassV0TrackDCA_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   100, 0, 10);
	fHistSelV0MassV0DecayRadius_TrueK0s =
          new TH2F("fHistSelV0MassV0DecayRadius_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   500, 0., 100.);
	fHistSelV0MassV0PropLifeTime_TrueK0s =
          new TH2F("fHistSelV0MassV0PropLifeTime_TrueK0s", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   250, 0, 50.);
	fOutputList->Add(fHistSelArmenteros_TrueK0s);
	fOutputList->Add(fHistSelV0MassDecayLength_TrueK0s);
	fOutputList->Add(fHistSelV0MassPointingAngle_TrueK0s);
	fOutputList->Add(fHistSelV0MassV0DCA_TrueK0s);
	fOutputList->Add(fHistSelV0MassV0TrackDCA_TrueK0s);
	fOutputList->Add(fHistSelV0MassV0DecayRadius_TrueK0s);
	fOutputList->Add(fHistSelV0MassV0PropLifeTime_TrueK0s);

	fTreeK0s = new TTree("fTreeK0s","");
	fTreeK0s->Branch("TrueMass",&TrueMass,"TrueMass/D");
	fTreeK0s->Branch("TrueRap",&TrueRap,"TrueRap/D");
	fTreeK0s->Branch("TruePhi",&TruePhi,"TruePhi/D");
	fTreeK0s->Branch("TruePt",&TruePt,"TruePt/D");
	fTreeK0s->Branch("RecMass",&RecMass,"RecMass/D");
	fTreeK0s->Branch("RecRap",&RecRap,"RecRap/D");
	fTreeK0s->Branch("RecPhi",&RecPhi,"RecPhi/D");
	fTreeK0s->Branch("RecPt",&RecPt,"RecPt/D");
	fTreeK0s->Branch("isDetect",&isDetect,"isDetect/O");
	fOutputList->Add(fTreeK0s);

      }
    }
    if (fIsPrimTrackPairAna) {
      fHistULSPairMassPt_TightCut =
          new TH2F("fHistULSPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistULSPairMassPt_TightCut);
      fHistLSppPairMassPt_TightCut =
          new TH2F("fHistLSppPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSppPairMassPt_TightCut);
      fHistLSmmPairMassPt_TightCut =
          new TH2F("fHistLSmmPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSmmPairMassPt_TightCut);
    }

    fHistTrackP = new TH1F("fHistTrackP", "", 100, 0, 10);
    fHistTrackPt = new TH1F("fHistTrackPt", "", 100, 0, 10);
    fHistTrackEta = new TH1F("fHistTrackEta", "", 20, -1, 1);
    fHistTPCNClusts = new TH1F("fHistTPCNClusts", "", 160, 0, 160);
    fHistSPDNClusts = new TH1F("fHistSPDNClusts", "", 6, 0, 6);
    fHistTPCCrossRowsFindableRatio =
        new TH1F("fHistTPCCrossRowsFindableRatio", "", 200, 0, 2);
    fHistReducedChi2TPC = new TH1F("fHistReducedChi2TPC", "", 100, 0, 10);
    fHistReducedChi2ITS = new TH1F("fHistReducedChi2ITS", "", 250, 0, 50);
    fHistDCAz = new TH1F("fHistDCAz", "", 100, 0, 10);
    fHistDCAxyPt = new TH2F("fHistDCAxyPt", "", 200, -1, 1, 100, 0, 10);
    fHistOpeningAngle = new TH1F("fHistOpeningAngle","",100,0,1);
    
    fOutputList->Add(fHistTrackP);
    fOutputList->Add(fHistTrackPt);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTPCNClusts);
    fOutputList->Add(fHistSPDNClusts);
    fOutputList->Add(fHistTPCCrossRowsFindableRatio);
    fOutputList->Add(fHistReducedChi2TPC);
    fOutputList->Add(fHistReducedChi2ITS);
    fOutputList->Add(fHistDCAz);
    fOutputList->Add(fHistDCAxyPt);
    fOutputList->Add(fHistOpeningAngle);

    double min_dEdx = 0.;
    double max_dEdx = 300.;
    double width_dEdx = 0.1;

    double min_beta = 0.;
    double max_beta = 1.2;
    double width_beta = 0.01;

    double min_p = 0.;
    double max_p = 5.;
    double width_p = 0.1;

    double min_sigma = -10;
    double max_sigma = +10;
    double width_sigma = 0.1;

    fHistTPCdEdxP =
        new TH2F("fHistTPCdEdxP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
    fHistBetaP =
        new TH2F("fHistBetaP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_beta - min_beta) / width_beta, min_beta, max_beta);
    fHistTPCSigmaElectron = new TH2F(
        "fHistTPCSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaElectron = new TH2F(
        "fHistTOFSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaMuon = new TH2F(
        "fHistTPCSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaMuon = new TH2F(
        "fHistTOFSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaPion = new TH2F(
        "fHistTPCSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaPion = new TH2F(
        "fHistTOFSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaKaon = new TH2F(
        "fHistTPCSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaKaon = new TH2F(
        "fHistTOFSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaProton = new TH2F(
        "fHistTPCSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaProton = new TH2F(
        "fHistTOFSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);

    fHistSelTPCdEdxP =
        new TH2F("fHistSelTPCdEdxP", "", (max_p - min_p) / width_p, min_p,
                 max_p, (max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
    fHistSelBetaP =
        new TH2F("fHistSelBetaP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_beta - min_beta) / width_beta, min_beta, max_beta);
    fHistSelTPCSigmaElectron = new TH2F(
        "fHistSelTPCSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaElectron = new TH2F(
        "fHistSelTOFSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaMuon = new TH2F(
        "fHistSelTPCSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaMuon = new TH2F(
        "fHistSelTOFSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaPion = new TH2F(
        "fHistSelTPCSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaPion = new TH2F(
        "fHistSelTOFSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaKaon = new TH2F(
        "fHistSelTPCSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaKaon = new TH2F(
        "fHistSelTOFSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaProton = new TH2F(
        "fHistSelTPCSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaProton = new TH2F(
        "fHistSelTOFSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);

    fOutputList->Add(fHistTPCdEdxP);
    fOutputList->Add(fHistBetaP);
    fOutputList->Add(fHistTPCSigmaElectron);
    fOutputList->Add(fHistTOFSigmaElectron);
    fOutputList->Add(fHistTPCSigmaMuon);
    fOutputList->Add(fHistTOFSigmaMuon);
    fOutputList->Add(fHistTPCSigmaPion);
    fOutputList->Add(fHistTOFSigmaPion);
    fOutputList->Add(fHistTPCSigmaKaon);
    fOutputList->Add(fHistTOFSigmaKaon);
    fOutputList->Add(fHistTPCSigmaProton);
    fOutputList->Add(fHistTOFSigmaProton);

    fOutputList->Add(fHistSelTPCdEdxP);
    fOutputList->Add(fHistSelBetaP);
    fOutputList->Add(fHistSelTPCSigmaElectron);
    fOutputList->Add(fHistSelTOFSigmaElectron);
    fOutputList->Add(fHistSelTPCSigmaMuon);
    fOutputList->Add(fHistSelTOFSigmaMuon);
    fOutputList->Add(fHistSelTPCSigmaPion);
    fOutputList->Add(fHistSelTOFSigmaPion);
    fOutputList->Add(fHistSelTPCSigmaKaon);
    fOutputList->Add(fHistSelTOFSigmaKaon);
    fOutputList->Add(fHistSelTPCSigmaProton);
    fOutputList->Add(fHistSelTOFSigmaProton);
  }

  fHistEventVtxZ = new TH1F("fHistEventVtxZ", "", 60, -30, 30);
  fHistEventCent = new TH1F("fHistEventCent", "", 100, 0, 100);
  fHistEventMulti = new TH1F("fHistEventMulti", "", 200, 0, 200);
  fHistEventVtxCont = new TH1F("fHistEventVtxCont", "", 100, 0, 100);
  fOutputList->Add(fHistEventVtxZ);
  fOutputList->Add(fHistEventCent);
  fOutputList->Add(fHistEventMulti);
  fOutputList->Add(fHistEventVtxCont);

  PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *) {
  
  if (!Initialize()) {
    return;
  }
  if (!fUtils->isAcceptEvent()) {
    return;
  }

  EventQA();
  
  if(fIsMC){
    ProcessMC();
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Fwd muon
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (!fIsMidTrackAna) {
    if (fIsMixingAnalysis) {
      FwdMuonPairAnalysisEveMixing();
    } else {
      FwdMuonPairAnalysis();
    }
  } else {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Mid rapidity analysis
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (!fIsMixingAnalysis) {
      if (fIsV0TrackPairAna) {
	MidV0Analysis(fUtils->getPairTargetPIDs(0),
                      fUtils->getPairTargetPIDs(1));
      } else if (fIsPrimTrackPairAna) {
        MidPairAnalysis(fUtils->getPairTargetPIDs(0),
                        fUtils->getPairTargetPIDs(1));
      }
    } else if (fIsMixingAnalysis) {
      if (fIsV0TrackPairAna) {
        MidV0AnalysisEventMixing(fUtils->getPairTargetPIDs(0),
                                 fUtils->getPairTargetPIDs(1));
      } else if (fIsPrimTrackPairAna) {
        MidPairAnalysisEventMixing(fUtils->getPairTargetPIDs(0),
				   fUtils->getPairTargetPIDs(1));
      }
    }
  }


}

bool AliAnalysisTaskAODTrackPair::Initialize() {
  fEvent = dynamic_cast<AliAODEvent *>(InputEvent());

  if (!fUtils->setEvent(fEvent, fInputHandler)) {
    return false;
  }
  
  if (fRunNumber != fEvent->GetRunNumber()) {
    fRunNumber = fUtils->getRunnumber();
    if (!fIsMidTrackAna) {
      AliMuonTrackCuts *trackCut = fUtils->getMuonTrackCuts();
      trackCut->SetRun(fInputHandler);
    }
  }

  if (!fIsMidTrackAna) {
    fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);
  }

  if ( fIsMC ) {
    
    fMCTrackArray = dynamic_cast<TClonesArray *>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));

    if (!fMCTrackArray) {
      return false;
    }

    fUtils->setMCArray(fMCTrackArray);
    
    AliAODMCHeader *mcHeader = (AliAODMCHeader *)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    
    if (!mcHeader){
      return false;
    }

    TList *headerList = mcHeader->GetCocktailHeaders();

    if (!headerList){
      return false;
    }

    for (Int_t i = 0; i < headerList->GetEntries(); i++) {
      AliGenEventHeader *eventHeader2 = (AliGenEventHeader *)headerList->At(i);
      TString name = eventHeader2->GetName();
      if (name.Contains("EtaDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("EtaDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("OmegaDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("OmegaDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("EtaPrimeDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("PhiDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("RhoDirect")) {
        fUtils->set2BodyProd(true);
      } else {
        fUtils->setDalitzProd(false);
        fUtils->set2BodyProd(false);
      }
    }
    
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::EventQA() {
  fHistEventVtxZ->Fill(fUtils->getVtxZ());
  fHistEventCent->Fill(fUtils->getCentClass());
  fHistEventMulti->Fill(fUtils->getNCorrSPDTrkInfo(1));
  fHistEventVtxCont->Fill(fUtils->getVtxCont());
  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonTrackQA(AliAODTrack *track) {
  fHistTrackEta->Fill(track->Pt(), fabs(track->Eta()));  
  fHistTrackThetaAbs->Fill(track->Pt(),AliAnalysisMuonUtility::GetThetaAbsDeg(track));  
  fHistTrackTriggerMatch->Fill(track->Pt(),AliAnalysisMuonUtility::GetMatchTrigger(track));  
  fHistTrackChiSquare->Fill(track->Pt(), AliAnalysisMuonUtility::GetChi2perNDFtracker(track));  
  fHistTriggerChiSquare->Fill(track->Pt(), AliAnalysisMuonUtility::GetChi2MatchTrigger(track));

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairQA(AliAODDimuon *dimuon) {

  AliAODTrack *track1 = dynamic_cast<AliAODTrack *>(dimuon->GetMu(0));
  AliAODTrack *track2 = dynamic_cast<AliAODTrack *>(dimuon->GetMu(1));

  int triggerLB1 = AliAnalysisMuonUtility::GetLoCircuit(track1);
  int triggerLB2 = AliAnalysisMuonUtility::GetLoCircuit(track2);

  double pt_max = track1->Pt();
  double pt_min = track2->Pt();

  if (track1->Pt() > track2->Pt()) {
    pt_max = track1->Pt();
    pt_min = track2->Pt();
  } else {
    pt_max = track2->Pt();
    pt_min = track1->Pt();
  }

  fHistTrackPairPtBalance->Fill(pt_max, pt_min);
  fHistTrackPairLocalBoardPair->Fill(triggerLB1, triggerLB2);

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysisEveMixing() {
  
  if (!(fInputHandler->IsEventSelected() & fTriggerMaskForMixing))
    return false;

  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  double poolCent = 0.;
  double poolVtxZ = 0.;
  double poolPsi = 0.;

  if (onEvtMixingPoolVtxZ) {
    poolVtxZ = fUtils->getVtxZ();
  }
  if (onEvtMixingPoolCent) {
    poolCent = fUtils->getCentClass();
  }
  if (onEvtMixingPoolPsi) {
    poolPsi = fUtils->getPsi();
  }

  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgr->GetEventPool(
								       poolCent, poolVtxZ, poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    AliAODTrack *track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptFwdMuonTrack(track1))
      continue;

    if (pool->IsReady()) {

      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {

        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

        for (Int_t iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast();
             ++iTrack2) {

          AliAODTrack *__track2__ = (AliAODTrack *)poolTracks->At(iTrack2);

          AliAODTrack *track2 = (AliAODTrack *)__track2__->Clone();

          if (!fUtils->isAcceptFwdMuonTrack(track2))
            continue;

          AliAODDimuon *dimuon = new AliAODDimuon();
          dimuon->SetMuons(track1, track2);

          if (!fUtils->isAcceptFwdDimuon(dimuon))
            continue;

          RecPairPt = dimuon->Pt();
          RecPairRap = fabs(dimuon->Y());
          RecPairMass = dimuon->M();
          RecPairCent = fUtils->getCentClass();
          RecPairDS = fUtils->getDS();

          string fFiredTrigName = string(fEvent->GetFiredTriggerClasses());

          if (dimuon->Charge() == 0) {
            fTreeMixULSPair->Fill();
          } else if (dimuon->Charge() > 0) {
            fTreeMixLSppPair->Fill();
          } else {
            fTreeMixLSmmPair->Fill();
          }

          delete track2;
          delete dimuon;

        } // end of loop track2
      }   // end of loop iMixEvt
    }     // poolPion->IsReady()

    fTrackArray->Add(track1);

  } // end of loop track1

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysis() {

  if (!fIsMC && !(fInputHandler->IsEventSelected() & fTriggerMaskForSame)) {
    return false;
  }
  
  if (fIsCMUL7) {
    fEventCounter->Fill(0., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(4., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMLL7) {
    fEventCounter->Fill(1., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(5., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 || fIsCMLL7) {
    fEventCounter->Fill(2., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(6., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 && fIsCMLL7) {
    fEventCounter->Fill(3., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(7., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }

  
  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;

  AliAODDimuon *dimuon;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {
    
    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptFwdMuonTrack(track1)){
      continue;
    }
    
    FwdMuonTrackQA(track1);

    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTrack; ++iTrack2) {
      
      track2 = (AliAODTrack *)fEvent->GetTrack(iTrack2);
      
      if (!fUtils->isAcceptFwdMuonTrack(track2)){
        continue;
      }
      
      dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1, track2);

      if (!fUtils->isAcceptFwdDimuon(dimuon)){
        continue;
      }
      
      FwdMuonPairQA(dimuon);

      RecPairPt = dimuon->Pt();
      RecPairRap = fabs(dimuon->Y());
      RecPairMass = dimuon->M();
      RecPairCent = fUtils->getCentClass();
      RecPairDS = fUtils->getDS();

      if (dimuon->Charge() == 0) {
        fTreeULSPair->Fill();
      } else if (dimuon->Charge() > 0) {
        fTreeLSppPair->Fill();
      } else {
        fTreeLSmmPair->Fill();
      }

      delete dimuon;

    } // end of loop track2
  }   // end of loop track1
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidTrackQualityChecker(AliAODTrack *track) {

  fHistTPCNClusts->Fill(track->GetTPCNcls());

  int nSPD = 0;

  if (track->HasPointOnITSLayer(0)) {
    ++nSPD;
  }
  if (track->HasPointOnITSLayer(1)) {
    ++nSPD;
  }
  fHistSPDNClusts->Fill(nSPD);
  fHistTPCCrossRowsFindableRatio->Fill(track->GetTPCNclsF() /
                                       track->GetTPCCrossedRows());
  fHistReducedChi2TPC->Fill(track->GetTPCchi2() / track->GetTPCNcls());

  if (track->GetITSNcls() > 0) {
    fHistReducedChi2ITS->Fill(track->GetITSchi2() / track->GetITSNcls());
  }

  float dca_xy = 9999;
  float dca_z = 9999;
  track->GetImpactParameters(dca_xy, dca_z);

  fHistDCAz->Fill(dca_z);
  fHistDCAxyPt->Fill(dca_xy, track->Pt());

  fHistTrackP->Fill(track->P());
  fHistTrackPt->Fill(track->Pt());
  fHistTrackEta->Fill(track->Eta());

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0Checker(AliAODv0 *v0, bool isSel) {

  if(!v0) {
    return false;
  }

  double vtx[] = {fUtils->getVtxX(), fUtils->getVtxY(), fUtils->getVtxZ()};

  if (isSel) {
    fHistULSPairMassPt_ProngV0->Fill(v0->MassK0Short(),v0->Pt());
    fHistSelArmenteros->Fill(v0->Alpha(), v0->PtArmV0());
    fHistSelV0MassDecayLength->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
    fHistSelV0MassPointingAngle->Fill(v0->MassK0Short(),v0->CosPointingAngle(vtx));
    fHistSelV0MassV0DCA->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
    fHistSelV0MassV0TrackDCA->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
    fHistSelV0MassV0DecayRadius->Fill(v0->MassK0Short(), v0->RadiusV0());
    fHistSelV0MassV0PropLifeTime->Fill(v0->MassK0Short(), fUtils->fPdgK0sMass * v0->DecayLengthV0(vtx) / v0->P());
    if(fIsMC && 
       fUtils->isSameMotherPair((AliAODTrack *)v0->GetDaughter(0),(AliAODTrack *)v0->GetDaughter(1)) &&
       fUtils->getMotherPdgCode((AliAODTrack *)v0->GetDaughter(0)) == 310){
      fHistULSPairMassPt_ProngV0_TrueK0s->Fill(v0->MassK0Short(),v0->Pt());
      fHistSelArmenteros_TrueK0s->Fill(v0->Alpha(), v0->PtArmV0());
      fHistSelV0MassDecayLength_TrueK0s->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
      fHistSelV0MassPointingAngle_TrueK0s->Fill(v0->MassK0Short(),v0->CosPointingAngle(vtx));
      fHistSelV0MassV0DCA_TrueK0s->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
      fHistSelV0MassV0TrackDCA_TrueK0s->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
      fHistSelV0MassV0DecayRadius_TrueK0s->Fill(v0->MassK0Short(), v0->RadiusV0());
      fHistSelV0MassV0PropLifeTime_TrueK0s->Fill(v0->MassK0Short(), fUtils->fPdgK0sMass * v0->DecayLengthV0(vtx) / v0->P());      
    }
  } else {
    fHistArmenteros->Fill(v0->Alpha(), v0->PtArmV0());
    fHistV0MassDecayLength->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
    fHistV0MassPointingAngle->Fill(v0->MassK0Short(),
                                   v0->CosPointingAngle(vtx));
    fHistV0MassV0DCA->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
    fHistV0MassV0TrackDCA->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
    fHistV0MassV0DecayRadius->Fill(v0->MassK0Short(), v0->RadiusV0());
    fHistV0MassV0PropLifeTime->Fill(v0->MassK0Short(),
                                    fUtils->fPdgK0sMass *
				    v0->DecayLengthV0(vtx) / v0->P());
    if(fIsMC && 
       fUtils->isSameMotherPair((AliAODTrack *)v0->GetDaughter(0),(AliAODTrack *)v0->GetDaughter(1)) &&
       fUtils->getMotherPdgCode((AliAODTrack *)v0->GetDaughter(0)) == 310){
      fHistArmenteros_TrueK0s->Fill(v0->Alpha(), v0->PtArmV0());
      fHistV0MassDecayLength_TrueK0s->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
      fHistV0MassPointingAngle_TrueK0s->Fill(v0->MassK0Short(),v0->CosPointingAngle(vtx));
      fHistV0MassV0DCA_TrueK0s->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
      fHistV0MassV0TrackDCA_TrueK0s->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
      fHistV0MassV0DecayRadius_TrueK0s->Fill(v0->MassK0Short(), v0->RadiusV0());
      fHistV0MassV0PropLifeTime_TrueK0s->Fill(v0->MassK0Short(), fUtils->fPdgK0sMass * v0->DecayLengthV0(vtx) / v0->P());      
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidTrackPIDChecker(AliAODTrack *track,
                                                     AliPID::EParticleType pid,
                                                     bool isSel) {

  double p = track->P();
  double sigTOF = track->GetTOFsignal();
  double length = track->GetIntegratedLength();
  double beta =
      (sigTOF > 0) ? (double)length / (2.99792457999999984e-02 * sigTOF) : -999;
  double dEdx = track->GetTPCsignal();
  
  if (isSel) {
    fHistSelTPCdEdxP->Fill(p, dEdx);
    fHistSelTPCSigmaElectron->Fill(p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistSelTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistSelTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistSelTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistSelTPCSigmaProton->Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistSelBetaP->Fill(p, beta);
      fHistSelTOFSigmaElectron->Fill(p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistSelTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistSelTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistSelTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistSelTOFSigmaProton->Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
    }
  } else {
    fHistTPCdEdxP->Fill(p, dEdx);
    fHistTPCSigmaElectron->Fill(p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistTPCSigmaProton->Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistBetaP->Fill(p, beta);
      fHistTOFSigmaElectron->Fill(p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistTOFSigmaProton->Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0Analysis(AliPID::EParticleType pid1,
                                                AliPID::EParticleType pid2) {  
  
  /*
  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliExternalTrackParam *exTrack;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {
    
    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);
    
    unique_ptr<AliExternalTrackParam> exTrack(new AliExternalTrackParam());
        
  }
  */

  Int_t nV0 = fEvent->GetNumberOfV0s();

  AliAODv0 *v0_1;
  AliAODv0 *v0_2;
  
  TLorentzVector lv1, lv2, lv12, lv_leading;
  
  int iLeading = -1;  
  AliAODTrack *leadingTrack = 0x0;
  
   if (fUtils->getLeadingTrack(iLeading)) {
  
     leadingTrack = (AliAODTrack*)fEvent->GetTrack(iLeading);
  
     lv_leading.SetPtEtaPhiM(leadingTrack->Pt(),leadingTrack->Eta(),leadingTrack->Phi(),leadingTrack->M());
   }
  
   for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {

     v0_1 = (AliAODv0 *)fEvent->GetV0(iV0_1);

     if (!fUtils->isAcceptV0Kinematics(v0_1)) {
       continue;
     }

     RecPairPt = v0_1->Pt();
     RecPairMass = v0_1->MassK0Short();
     RecPairRap = v0_1->RapK0Short();

     AliAODTrack *pTrack1 = (AliAODTrack *)v0_1->GetDaughter(0);
     AliAODTrack *nTrack1 = (AliAODTrack *)v0_1->GetDaughter(1);
     
     if (fUtils->isAcceptV0Basic(v0_1, 0)) {
       MidV0Checker(v0_1, false);    
       MidTrackPIDChecker(pTrack1, pid1, false);
       MidTrackPIDChecker(nTrack1, pid2, false);
     }
     if (0.4 > RecPairMass || RecPairMass > 0.6) {
       continue;
     }    
     if (!fUtils->isAcceptK0s(v0_1)) {
       continue;
     }
     
     bool isKmpCand1 = false;
     double MKpm = 0;
     if (fUtils->isAcceptedK0sFromKpmStar(v0_1,MKpm)) {
       fHistKpmStarCandMass->Fill(MKpm);
       isKmpCand1 = true;
     }
         
     MidTrackQualityChecker(pTrack1);
     MidTrackQualityChecker(nTrack1);
     
     MidV0Checker(v0_1, true);
     MidTrackPIDChecker(pTrack1, pid1, true);
     MidTrackPIDChecker(nTrack1, pid2, true);

     for (int iV0_2 = iV0_1 + 1; iV0_2 < nV0; ++iV0_2) {

       v0_2 = (AliAODv0 *)fEvent->GetV0(iV0_2);

       RecPairPt = v0_2->Pt();
       RecPairMass = v0_2->MassK0Short();
       RecPairRap = v0_2->RapK0Short();

       if (0.4 > RecPairMass || RecPairMass > 0.6) {
	 continue;
       }
       if (!fUtils->isAcceptV0Kinematics(v0_2)) {
	 continue;
       }
       if (!fUtils->isAcceptK0s(v0_2)) {
	 continue;
       }

       bool isKmpCand2 = false;
       if (fUtils->isAcceptedK0sFromKpmStar(v0_2,MKpm)) {
	 isKmpCand2 = true;
       }
       
       fHistMassK0s1K0s2->Fill(v0_1->MassK0Short(), v0_2->MassK0Short());
       
       if (!fUtils->isAcceptK0sK0sOpeningAngle(v0_1,v0_2)){
	 continue;
       }
       if (!fUtils->isAcceptNotSharingTracks(v0_1,v0_2)){
	 continue;
       }

       lv1.SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(),
			TDatabasePDG::Instance()->GetParticle(310)->Mass());
       lv2.SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(),
			TDatabasePDG::Instance()->GetParticle(310)->Mass());
       lv12 = lv1 + lv2;
       
       double angle = lv1.Angle(lv2.Vect());
       
       RecPairPt = lv12.Pt();
       RecPairMass = lv12.M();
       RecPairRap = lv12.Rapidity();

       double fill[] = {RecPairMass, RecPairPt, fUtils->getCentClass()};

       if (fUtils->isAcceptK0sCandidateMassRange(v0_1->MassK0Short()) && fUtils->isAcceptK0sCandidateMassRange(v0_2->MassK0Short())) {
	 
	 fSparseULSPairMassPt->Fill(fill);
	 if (!isKmpCand1 && !isKmpCand2) {
	   fSparseULSPairMassPt_RejectKpmStar->Fill(fill);
	 }
	 
	 fHistOpeningAngle->Fill(TMath::Cos(angle));
	 
	 if (iLeading>0) {	   
	   double leading_angle = lv_leading.Angle(lv12.Vect()) * 180./TMath::Pi();
	   double fill_leading[] = {RecPairMass, RecPairPt, leading_angle};	   
	   fSparseULSPairMassPt_LeadingTrack->Fill(fill_leading);
	 }
       } else { 
	 if (fUtils->isAcceptK0sCandidateSideBandRight(v0_1->MassK0Short()) && 
	     fUtils->isAcceptK0sCandidateSideBandRight(v0_2->MassK0Short())) {
	   fSparseULSPairMassPt_SideBandRight->Fill(fill);
	 } else if (fUtils->isAcceptK0sCandidateSideBandLeft(v0_1->MassK0Short()) && 
		    fUtils->isAcceptK0sCandidateSideBandLeft(v0_2->MassK0Short())) {
	   fSparseULSPairMassPt_SideBandLeft->Fill(fill);
	 } else if ((fUtils->isAcceptK0sCandidateSideBandLeft(v0_1->MassK0Short()) && 
		     fUtils->isAcceptK0sCandidateSideBandRight(v0_2->MassK0Short())) || 
		    (fUtils->isAcceptK0sCandidateSideBandRight(v0_1->MassK0Short()) &&
		     fUtils->isAcceptK0sCandidateSideBandLeft(v0_2->MassK0Short()))) {
	   fSparseULSPairMassPt_SideBandLeftRight->Fill(fill);
	 } 
	 if ((fUtils->isAcceptK0sCandidateSideBandRight(v0_1->MassK0Short()) || 
	      fUtils->isAcceptK0sCandidateSideBandLeft(v0_1->MassK0Short())) && 
	     (fUtils->isAcceptK0sCandidateSideBandRight(v0_2->MassK0Short()) ||
	      fUtils->isAcceptK0sCandidateSideBandLeft(v0_2->MassK0Short()))) {	    
	   fSparseULSPairMassPt_SideBand->Fill(fill);
	 }
       }
     }
   }
   return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0AnalysisEventMixing(
    AliPID::EParticleType pid1, AliPID::EParticleType pid2) {
  
  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  double poolCent = 0.;
  double poolVtxZ = 0.;
  double poolPsi = 0.;

  TLorentzVector lv1, lv2, lv12, lv_leading;
  
  int iLeading = -1;  
  AliAODTrack *leadingTrack = 0x0;

  if (fUtils->getLeadingTrack(iLeading)) {
    leadingTrack = (AliAODTrack*)fEvent->GetTrack(iLeading);
    lv_leading.SetPtEtaPhiM(leadingTrack->Pt(),leadingTrack->Eta(),leadingTrack->Phi(),leadingTrack->M());
  }
  
  if (onEvtMixingPoolVtxZ) {
    poolVtxZ = fUtils->getVtxZ();
  }
  if (onEvtMixingPoolCent) {
    poolCent = fUtils->getCentClass();
  }
  
  if (onEvtMixingPoolPsi) {
    if (iLeading>0) {
      poolPsi = leadingTrack->Phi();
    } else {
      poolPsi = fUtils->getPsi();
    }  
    poolPsi = fUtils->getPsi();
  }
  
  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgr->GetEventPool(poolCent, poolVtxZ, poolPsi);
  
  Int_t nV0 = fEvent->GetNumberOfV0s();

  AliAODv0 *v0_1;
  AliAODv0 *v0_2;
  
  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {
    
    v0_1 = static_cast<AliAODv0 *>(fEvent->GetV0(iV0_1));
    
    if (!fUtils->isAcceptV0Kinematics(v0_1)) {
      continue;
    }
    
    RecPairPt = v0_1->Pt();
    RecPairMass = v0_1->MassK0Short();
    RecPairRap = v0_1->RapK0Short();

    AliAODTrack *pTrack = (AliAODTrack *)v0_1->GetDaughter(0);
    AliAODTrack *nTrack = (AliAODTrack *)v0_1->GetDaughter(1);

    if (fUtils->isAcceptV0Basic(v0_1, 0)) {
      MidV0Checker(v0_1, false);    
      MidTrackPIDChecker(pTrack, pid1, false);
      MidTrackPIDChecker(nTrack, pid2, false);
    }

    if (0.4 > RecPairMass || RecPairMass > 0.6) {
      continue;
    }
    if (!fUtils->isAcceptK0s(v0_1)) {
      continue;
    }
    
    bool isKmpCand1 = false;
    double MKpm = 0;
    if (fUtils->isAcceptedK0sFromKpmStar(v0_1,MKpm)) {
      fHistKpmStarCandMass->Fill(MKpm);
      isKmpCand1 = true;
    }
    
    MidV0Checker(v0_1, true);

    MidTrackQualityChecker(pTrack);
    MidTrackQualityChecker(nTrack);
    
    if (!fUtils->isAcceptK0sCandidateMassRange(v0_1->MassK0Short())) {
      continue;
    }
    
    v0_1->SetOnFlyStatus(isKmpCand1);
    MidTrackPIDChecker(pTrack, pid1, true);
    MidTrackPIDChecker(nTrack, pid2, true);
    
    if (pool->IsReady()) {
      
      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {
	
        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);
	
        for (int iV0_2 = 0; iV0_2 < poolTracks->GetEntriesFast(); ++iV0_2) {
	  
	  v0_2 = static_cast<AliAODv0 *>(poolTracks->At(iV0_2));
	  
	  if (!fUtils->isAcceptK0sK0sOpeningAngle(v0_1,v0_2)){
	    continue;
	  }
	  
	  bool isKmpCand2 = false;	  
	  if (v0_2->GetOnFlyStatus()) {
	    isKmpCand2 = true;
	  }
	  
          lv1.SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(),
                           TDatabasePDG::Instance()->GetParticle(310)->Mass());
	  lv2.SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(),
                           TDatabasePDG::Instance()->GetParticle(310)->Mass());
          lv12 = lv1 + lv2;
	  
	  double angle = lv1.Angle(lv2.Vect());

          RecPairPt = lv12.Pt();
          RecPairMass = lv12.M();
          RecPairRap = lv12.Rapidity();

          double fill[] = {RecPairMass, RecPairPt, fUtils->getCentClass()};

	  fSparseMixULSPairMassPt->Fill(fill);
	  if (!isKmpCand1 && !isKmpCand2) {
	    fSparseMixULSPairMassPt_RejectKpmStar->Fill(fill);
	  }
	  fHistOpeningAngle->Fill(TMath::Cos(angle));
	  
	  if (iLeading>0) {
	    double leading_angle = lv_leading.Angle(lv12.Vect()) * 180./TMath::Pi();
	    double fill_leading[] = {RecPairMass, RecPairPt, leading_angle};	   
	    fSparseMixULSPairMassPt_LeadingTrack->Fill(fill_leading);
	  }
	    	  
        }
      }
    }
    
    fTrackArray->Add(v0_1);
  }

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidPairAnalysis(AliPID::EParticleType pid1,
                                                  AliPID::EParticleType pid2) {
  
  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;

  TLorentzVector lv1, lv2, lv12;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptTrackKinematics(track1)) {
      continue;
    }
    if (!fUtils->isAcceptMidPrimTrackQuality(track1)) {
      continue;
    }
    
    MidTrackPIDChecker(track1, pid1, false);

    if (!fUtils->isAcceptMidPid(track1, pid1)) {
      continue;
    }

    MidTrackPIDChecker(track1, pid1, true);
    MidTrackQualityChecker(track1);

    double mass1 = 0;
    
    if (pid1 == AliPID::kElectron) {
      mass1 = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    } else if (pid1 == AliPID::kPion) {
      mass1 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    } else if (pid1 == AliPID::kKaon) {
      mass1 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    } else if (pid1 == AliPID::kProton) {
      mass1 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    } else if (pid1 == AliPID::kMuon) {
      mass1 = TDatabasePDG::Instance()->GetParticle(13)->Mass();
    } else {
      continue;
    }
    
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTrack; ++iTrack2) {

      track2 = (AliAODTrack *)fEvent->GetTrack(iTrack2);

      if (!fUtils->isAcceptTrackKinematics(track2)) {
	continue;
      }
      if (!fUtils->isAcceptMidPrimTrackQuality(track2)) {
	continue;
      }
      if (!fUtils->isAcceptMidPid(track2, pid2)) {
	continue;
      }

      double mass2 = 0;
    
      if (pid2 == AliPID::kElectron) {
	mass2 = TDatabasePDG::Instance()->GetParticle(11)->Mass();
      } else if (pid2 == AliPID::kPion) {
	mass2 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
      } else if (pid2 == AliPID::kKaon) {
	mass2 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
      } else if (pid2 == AliPID::kProton) {
	mass2 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
      } else if (pid2 == AliPID::kMuon) {
	mass2 = TDatabasePDG::Instance()->GetParticle(13)->Mass();
      } else {
	continue;
      }
      	    
      lv1.SetPtEtaPhiM(track1->Pt(),track1->Eta(),track1->Phi(),mass1);
      lv2.SetPtEtaPhiM(track2->Pt(),track2->Eta(),track2->Phi(),mass2);
      
      lv12 = lv1 + lv2;

      double angle = lv1.Angle(lv2.Vect());

      double fill[] = {lv12.M(), lv12.Pt(), fUtils->getCentClass()};

      if (track1->Charge() + track2->Charge() == 0){
	fSparseULSPairMassPt->Fill(fill);
      } else if (track1->Charge() + track2->Charge() > 0){
	fSparseLSppPairMassPt->Fill(fill);
      } else {
	fSparseLSmmPairMassPt->Fill(fill);
      }
      
    } // end of loop track2
  }   // end of loop track1

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidPairAnalysisEventMixing(AliPID::EParticleType pid1, AliPID::EParticleType pid2) {
  //cout<<"MidPairAnalysisEventMixing"<<endl;
  double mass1 = 0;
  double mass2 = 0;

  if (pid1 == AliPID::kElectron) {
    mass1 = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  } else if (pid1 == AliPID::kPion) {
    mass1 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  } else if (pid1 == AliPID::kKaon) {
    mass1 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  } else if (pid1 == AliPID::kProton) {
    mass1 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  } else if (pid1 == AliPID::kMuon) {
    mass1 = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  } 

  if (pid2 == AliPID::kElectron) {
    mass2 = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  } else if (pid2 == AliPID::kPion) {
    mass2 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  } else if (pid2 == AliPID::kKaon) {
    mass2 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  } else if (pid2 == AliPID::kProton) {
    mass2 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  } else if (pid2 == AliPID::kMuon) {
    mass2 = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  } 

  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  double poolCent = 0.;
  double poolVtxZ = 0.;
  double poolPsi = 0.;

  if (onEvtMixingPoolVtxZ) {
    poolVtxZ = fUtils->getVtxZ();
  }
  if (onEvtMixingPoolCent) {
    poolCent = fUtils->getCentClass();
  }
  if (onEvtMixingPoolPsi) {
    poolPsi = fUtils->getPsi();
  }

  AliEventPool *pool = (AliEventPool*)fPoolMuonTrackMgr->GetEventPool(poolCent, poolVtxZ, poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;
  
  TLorentzVector lv1, lv2, lv12;

  for (int iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {
    
    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptMidPrimTrackQuality(track1)) {
      continue;
    }
    if (!fUtils->isAcceptTrackKinematics(track1)) {
      continue;
    }
    
    MidTrackPIDChecker(track1, pid1, false);
    
    if (!fUtils->isAcceptMidPid(track1, pid1)) {
      continue;
    }

    MidTrackPIDChecker(track1, pid1, true);
    MidTrackQualityChecker(track1);

    if (pool->IsReady()) {

      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {
	
        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

        for (int iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast(); iTrack2++) {

	  track2 = (AliAODTrack *)poolTracks->At(iTrack2);	  

          lv1.SetPtEtaPhiM(track1->Pt(), track1->Eta(), track1->Phi(),mass1);
	  lv2.SetPtEtaPhiM(track2->Pt(), track2->Eta(), track2->Phi(),mass2);
	  
          lv12 = lv1 + lv2;
	  
	  double angle = lv1.Angle(lv2.Vect());

          double fill[] = {lv12.M(), lv12.Pt(), fUtils->getCentClass()};
	  
	  if (track1->Charge() + track2->Charge() == 0){
	    fSparseMixULSPairMassPt->Fill(fill);
	  } else if (track1->Charge() + track2->Charge() > 0){
	    fSparseMixLSppPairMassPt->Fill(fill);
	  } else {
	    fSparseMixLSmmPairMassPt->Fill(fill);
	  }
	  
        }
      }
    }    

    fTrackArray->Add(track1);
  }

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::ProcessMC(){
  
  AliAODMCParticle *particle1;
  AliAODMCParticle *particle2;
  AliAODMCParticle *particle12;
  
  Int_t nV0 = fEvent->GetNumberOfV0s();
  
  vector<int> plabels1;
  vector<int> plabels2;
  vector<int> plabels12;
  vector<int> v0label;
  
  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {
    
    AliAODv0* v0_1 = (AliAODv0 *)fEvent->GetV0(iV0_1);
     
    if (!fUtils->isAcceptV0Kinematics(v0_1)) {
      continue;
    }
    if (0.4 > v0_1->MassK0Short() || v0_1->MassK0Short() > 0.6) {
      continue;
    }    
    if (!fUtils->isAcceptK0s(v0_1)) {
      continue;
    }

    AliAODTrack *Track1 = (AliAODTrack *)v0_1->GetDaughter(0);
    AliAODTrack *Track2 = (AliAODTrack *)v0_1->GetDaughter(1);
    
    plabels12.push_back(Track1->GetLabel());    
    plabels12.push_back(Track2->GetLabel());
    v0label.push_back(iV0_1);
    v0label.push_back(iV0_1);

    plabels1.push_back(Track1->GetLabel());    
    plabels2.push_back(Track2->GetLabel());    
  }
  
  for (Int_t iTrack1 = 0; iTrack1 < fMCTrackArray->GetEntries(); ++iTrack1) {
    
    TrueMass = 0;
    TrueRap = 0;
    TruePhi = 0;
    TruePt = 0;

    RecMass = 0;
    RecRap = 0;
    RecPhi = 0;
    RecPt = 0;

    isDetect = false;
    
    particle1 = (AliAODMCParticle *)fMCTrackArray->At(iTrack1);
    
    if(!particle1) {
      continue;
    }
    
    if (particle1->GetPdgCode() != 310 || particle1->GetNDaughters() != 2){
      continue;
    }
    
    int label1 = particle1->GetDaughterFirst();
    int label2 = particle1->GetDaughterLast();
    
    AliAODMCParticle* p1 = (AliAODMCParticle *)fMCTrackArray->At(label1);
    AliAODMCParticle* p2 = (AliAODMCParticle *)fMCTrackArray->At(label2);

    if(!p1 || !p2) {
      continue;
    }
    
    if (fabs(p1->GetPdgCode()) != 211 || fabs(p2->GetPdgCode()) != 211) {
      continue;
    }
    
    TrueMass = particle1->M();
    TrueRap = particle1->Y();
    TruePhi = particle1->Phi();
    TruePt = particle1->Pt();

    bool detect1=false;
    bool detect2=false;
    
    int irec1 = -1;
    int irec2 = -1;
    
    int v0label1=-1;
    int v0label2=-1;
    
    fHistTrueK0sPtRapidity->Fill(particle1->Y(),particle1->Pt());

    for (int iDetect = 0; iDetect<plabels12.size(); ++iDetect) {
      if (plabels12[iDetect] == label1) {
	detect1 = true;
	v0label1 = v0label[iDetect];
      }
      if (plabels12[iDetect] == label2) {
	detect2 = true;
	v0label2 = v0label[iDetect];
      }
    }
    
    if (detect1 && detect2 && v0label1==v0label2){

      fHistRecTrueK0sPtRapidity->Fill(particle1->Y(),particle1->Pt());

      AliAODv0* v0_1 = (AliAODv0 *)fEvent->GetV0(v0label1);
    
      RecMass = v0_1->MassK0Short();
      RecRap = v0_1->RapK0Short();
      RecPhi = v0_1->Phi();
      RecPt = v0_1->Pt();
      
      isDetect = true;
    }
    
    fTreeK0s->Fill();

  }
  

  return true;
}
