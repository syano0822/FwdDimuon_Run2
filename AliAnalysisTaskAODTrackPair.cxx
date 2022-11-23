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

#include "AliKFParticle.h"
#include "AliKFVertex.h"

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

#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"

#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskAODTrackPair.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

#include "iostream"
#include "memory"

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair() : AliAnalysisTaskSE(),
  
  fEvent(NULL),
  fPrimVtx(NULL),
  
  fPoolMuonTrackMgrK0s(NULL),  
  fPoolMuonTrackMgrPion(NULL),  
  
  fMultSelection(NULL),

  fUtils(NULL),
  
  fMCTrackArray(NULL),
  fArrayK0s(NULL),
  
  fHistDecayLengthXYCut(NULL),
  fHistPointingAngleXYCut(NULL),
  fHistChi2Cut(NULL),
  fHistDaughterPairDCAXYCut(NULL),
  fHistDaughterTrackDistanceXYCut(NULL),
  fHistDCAxyCut(NULL),
  fHistLifeTimeCut(NULL),
  
  fTrackDepth(100),
  fPoolSize(100),

  fReadyFraction(0.1),
  
  onVerbose(true),

  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  
  onDecayLengthXYCut(false),
  onPointingAngleXYCut(true),
  onChi2Cut(true),
  onDaughterPairDCAXYCut(false),
  onDaughterTrackDistanceXYCut(true),
  onDCAXYCut(false),
  onLifetimeCut(false),
  
  fIsMC(false),
  fIsMixingAnalysis(false),
  fIsManualV0Analysis(true),

  fMinK0sRap(-0.5),
  fMaxK0sRap(0.5),
  fMinK0sPt(0.4),
  fMaxK0sPt(99.),

  fMethodCent("SPDTracklets"),
  fPrimVtxPos(),
  fCent(0.),
  fPsi(0.),

  fMinNContPrimVtx(2),

  ////////////////////////////////////////////////
  // Output histos
  ////////////////////////////////////////////////

  fOutputList(NULL),
  
  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),
  
  fSparseULSPionPair(NULL),
  fSparseLSppPionPair(NULL),
  fSparseLSmmPionPair(NULL),

  fSparseMixULSPionPair(NULL),
  fSparseMixLSppPionPair(NULL),
  fSparseMixLSmmPionPair(NULL),
  
  fSparseULSPionPairBeforeCuts(NULL),
  fSparseLSppPionPairBeforeCuts(NULL),
  fSparseLSmmPionPairBeforeCuts(NULL),

  fSparseULSPionPair_PassChi2perNDFCut(NULL),
  fSparseULSPionPair_PassCPVXYCut(NULL),
  fSparseULSPionPair_PassDecayLengthXYCut(NULL),
  fSparseULSPionPair_PassDaughterDistanceXYCut(NULL),  
  fSparseULSPionPair_PassDCAXYCut(NULL),
  fSparseULSPionPair_PassLifetimeCut(NULL),
  
  fSparseNeutralK0sPair(NULL),
  fSparseNeutralNegativeK0sPair(NULL),
  fSparseNeutralPositiveK0sPair(NULL),
  fSparsePositiveK0sPair(NULL),
  fSparseNegativeK0sPair(NULL),
  fSparsePositiveNegativeK0sPair(NULL),

  fSparseTrueK0s(NULL),
  fSparseTrueK0sRecK0s(NULL),

  fHistTPCdEdxP(NULL),
  fHistTPCSigmaElectron(NULL),
  fHistTPCSigmaMuon(NULL),
  fHistTPCSigmaPion(NULL),
  fHistTPCSigmaKaon(NULL),
  fHistTPCSigmaProton(NULL),  
  fHistBetaP(NULL),  
  fHistTOFSigmaElectron(NULL),  
  fHistTOFSigmaMuon(NULL),  
  fHistTOFSigmaPion(NULL),
  fHistTOFSigmaKaon(NULL),
  fHistTOFSigmaProton(NULL),  
  
  fHistSelTPCdEdxP(NULL),
  fHistSelTPCSigmaElectron(NULL),
  fHistSelTPCSigmaMuon(NULL),
  fHistSelTPCSigmaPion(NULL),
  fHistSelTPCSigmaKaon(NULL),
  fHistSelTPCSigmaProton(NULL),  
  fHistSelBetaP(NULL),  
  fHistSelTOFSigmaElectron(NULL),  
  fHistSelTOFSigmaMuon(NULL),  
  fHistSelTOFSigmaPion(NULL),
  fHistSelTOFSigmaKaon(NULL),
  fHistSelTOFSigmaProton(NULL),
  
  fHistTPCNClusts(NULL),
  fHistSPDNClusts(NULL),
  fHistTPCCrossRowsFindableRatio(NULL),
  fHistReducedChi2TPC(NULL),
  fHistReducedChi2ITS(NULL),
  fHistDCAz(NULL),
  fHistDCAxyPt(NULL),
  fHistTrackP(NULL),
  fHistTrackPt(NULL),
  fHistTrackEta(NULL),    

  fHistV0PV0DecayLengthXY(NULL),
  fHistV0PV0DecayLength(NULL),
  fHistV0PV0PointingAngleXY(NULL),
  fHistV0PV0PointingAngle(NULL),
  fHistV0PV0DCAXY(NULL),
  fHistV0PV0DCA(NULL),
  fHistV0PV0TrackDistanceXY(NULL),
  fHistV0PV0TrackDistance(NULL),
  fHistV0PV0Chi2perNDF(NULL),
  fHistV0PV0PropLifeTime(NULL),
  fHistV0TrackDCAXY(NULL),
  fHistV0TrackDCA(NULL),

  fHistSelV0PV0DecayLengthXY(NULL),
  fHistSelV0PV0DecayLength(NULL),
  fHistSelV0PV0PointingAngleXY(NULL),
  fHistSelV0PV0PointingAngle(NULL),
  fHistSelV0PV0DCAXY(NULL),
  fHistSelV0PV0DCA(NULL),
  fHistSelV0PV0TrackDistanceXY(NULL),
  fHistSelV0PV0TrackDistance(NULL),
  fHistSelV0PV0Chi2perNDF(NULL),
  fHistSelV0PV0PropLifeTime(NULL),
  fHistSelV0TrackDCAXY(NULL),
  fHistSelV0TrackDCA(NULL),

  fHistV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistV0PV0DecayLength_TrueK0s(NULL),
  fHistV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistV0PV0PointingAngle_TrueK0s(NULL),
  fHistV0PV0DCAXY_TrueK0s(NULL),
  fHistV0PV0DCA_TrueK0s(NULL),
  fHistV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistV0PV0TrackDistance_TrueK0s(NULL),
  fHistV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistV0PV0PropLifeTime_TrueK0s(NULL),
  fHistV0TrackDCAXY_TrueK0s(NULL),
  fHistV0TrackDCA_TrueK0s(NULL),

  fHistSelV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistSelV0PV0DecayLength_TrueK0s(NULL),
  fHistSelV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistSelV0PV0PointingAngle_TrueK0s(NULL),
  fHistSelV0PV0DCAXY_TrueK0s(NULL),
  fHistSelV0PV0DCA_TrueK0s(NULL),
  fHistSelV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistSelV0PV0TrackDistance_TrueK0s(NULL),
  fHistSelV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistSelV0PV0PropLifeTime_TrueK0s(NULL),
  fHistSelV0TrackDCAXY_TrueK0s(NULL),
  fHistSelV0TrackDCA_TrueK0s(NULL),
  
  fHistLSV0PV0DecayLengthXY(NULL),
  fHistLSV0PV0DecayLength(NULL),
  fHistLSV0PV0PointingAngleXY(NULL),
  fHistLSV0PV0PointingAngle(NULL),
  fHistLSV0PV0DCAXY(NULL),
  fHistLSV0PV0DCA(NULL),
  fHistLSV0PV0TrackDistanceXY(NULL),
  fHistLSV0PV0TrackDistance(NULL),
  fHistLSV0PV0Chi2perNDF(NULL),
  fHistLSV0PV0PropLifeTime(NULL),
  fHistLSV0TrackDCAXY(NULL),
  fHistLSV0TrackDCA(NULL),

  fHistLSSelV0PV0DecayLengthXY(NULL),
  fHistLSSelV0PV0DecayLength(NULL),
  fHistLSSelV0PV0PointingAngleXY(NULL),
  fHistLSSelV0PV0PointingAngle(NULL),
  fHistLSSelV0PV0DCAXY(NULL),
  fHistLSSelV0PV0DCA(NULL),
  fHistLSSelV0PV0TrackDistanceXY(NULL),
  fHistLSSelV0PV0TrackDistance(NULL),
  fHistLSSelV0PV0Chi2perNDF(NULL),
  fHistLSSelV0PV0PropLifeTime(NULL),
  fHistLSSelV0TrackDCAXY(NULL),
  fHistLSSelV0TrackDCA(NULL),

  fHistLSV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistLSV0PV0DecayLength_TrueK0s(NULL),
  fHistLSV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistLSV0PV0PointingAngle_TrueK0s(NULL),
  fHistLSV0PV0DCAXY_TrueK0s(NULL),
  fHistLSV0PV0DCA_TrueK0s(NULL),
  fHistLSV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistLSV0PV0TrackDistance_TrueK0s(NULL),
  fHistLSV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistLSV0PV0PropLifeTime_TrueK0s(NULL),
  fHistLSV0TrackDCAXY_TrueK0s(NULL),
  fHistLSV0TrackDCA_TrueK0s(NULL),

  fHistLSSelV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistLSSelV0PV0DecayLength_TrueK0s(NULL),
  fHistLSSelV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistLSSelV0PV0PointingAngle_TrueK0s(NULL),
  fHistLSSelV0PV0DCAXY_TrueK0s(NULL),
  fHistLSSelV0PV0DCA_TrueK0s(NULL),
  fHistLSSelV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistLSSelV0PV0TrackDistance_TrueK0s(NULL),
  fHistLSSelV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistLSSelV0PV0PropLifeTime_TrueK0s(NULL),
  fHistLSSelV0TrackDCAXY_TrueK0s(NULL),
  fHistLSSelV0TrackDCA_TrueK0s(NULL),

  all(0),
  abord(0)
  
{

}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char *name) : AliAnalysisTaskSE(name), 

  fEvent(NULL),
  fPrimVtx(NULL),
  
  fPoolMuonTrackMgrK0s(NULL),  
  fPoolMuonTrackMgrPion(NULL),  
  
  fMultSelection(NULL),

  fUtils(NULL),
  
  fMCTrackArray(NULL),
  fArrayK0s(NULL),
  
  fHistDecayLengthXYCut(NULL),
  fHistPointingAngleXYCut(NULL),
  fHistChi2Cut(NULL),
  fHistDaughterPairDCAXYCut(NULL),
  fHistDaughterTrackDistanceXYCut(NULL),
  fHistDCAxyCut(NULL),
  fHistLifeTimeCut(NULL),

  fTrackDepth(100),
  fPoolSize(100),

  fReadyFraction(0.1),
  
  onVerbose(true),

  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),

  onDecayLengthXYCut(false),
  onPointingAngleXYCut(true),
  onChi2Cut(true),
  onDaughterPairDCAXYCut(false),
  onDaughterTrackDistanceXYCut(true),
  onDCAXYCut(false),
  onLifetimeCut(false),

  fIsMC(false),
  fIsMixingAnalysis(false),
  fIsManualV0Analysis(true),

  fMinK0sRap(-0.5),
  fMaxK0sRap(0.5),
  fMinK0sPt(0.4),
  fMaxK0sPt(99.),

  fMethodCent("SPDTracklets"),
  fPrimVtxPos(),
  fCent(0.),
  fPsi(0.),

  fMinNContPrimVtx(2),

  ////////////////////////////////////////////////
  // Output histos
  ////////////////////////////////////////////////

  fOutputList(NULL),
  
  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),
  
  fSparseULSPionPair(NULL),
  fSparseLSppPionPair(NULL),
  fSparseLSmmPionPair(NULL),

  fSparseMixULSPionPair(NULL),
  fSparseMixLSppPionPair(NULL),
  fSparseMixLSmmPionPair(NULL),
  
  fSparseULSPionPairBeforeCuts(NULL),
  fSparseLSppPionPairBeforeCuts(NULL),
  fSparseLSmmPionPairBeforeCuts(NULL),

  fSparseULSPionPair_PassChi2perNDFCut(NULL),
  fSparseULSPionPair_PassCPVXYCut(NULL),
  fSparseULSPionPair_PassDecayLengthXYCut(NULL),
  fSparseULSPionPair_PassDaughterDistanceXYCut(NULL),  
  fSparseULSPionPair_PassDCAXYCut(NULL),
  fSparseULSPionPair_PassLifetimeCut(NULL),

  fSparseNeutralK0sPair(NULL),
  fSparseNeutralNegativeK0sPair(NULL),
  fSparseNeutralPositiveK0sPair(NULL),
  fSparsePositiveK0sPair(NULL),
  fSparseNegativeK0sPair(NULL),
  fSparsePositiveNegativeK0sPair(NULL),

  fSparseTrueK0s(NULL),
  fSparseTrueK0sRecK0s(NULL),

  fHistTPCdEdxP(NULL),
  fHistTPCSigmaElectron(NULL),
  fHistTPCSigmaMuon(NULL),
  fHistTPCSigmaPion(NULL),
  fHistTPCSigmaKaon(NULL),
  fHistTPCSigmaProton(NULL),  
  fHistBetaP(NULL),  
  fHistTOFSigmaElectron(NULL),  
  fHistTOFSigmaMuon(NULL),  
  fHistTOFSigmaPion(NULL),
  fHistTOFSigmaKaon(NULL),
  fHistTOFSigmaProton(NULL),  
  
  fHistSelTPCdEdxP(NULL),
  fHistSelTPCSigmaElectron(NULL),
  fHistSelTPCSigmaMuon(NULL),
  fHistSelTPCSigmaPion(NULL),
  fHistSelTPCSigmaKaon(NULL),
  fHistSelTPCSigmaProton(NULL),  
  fHistSelBetaP(NULL),  
  fHistSelTOFSigmaElectron(NULL),  
  fHistSelTOFSigmaMuon(NULL),  
  fHistSelTOFSigmaPion(NULL),
  fHistSelTOFSigmaKaon(NULL),
  fHistSelTOFSigmaProton(NULL),
  
  fHistTPCNClusts(NULL),
  fHistSPDNClusts(NULL),
  fHistTPCCrossRowsFindableRatio(NULL),
  fHistReducedChi2TPC(NULL),
  fHistReducedChi2ITS(NULL),
  fHistDCAz(NULL),
  fHistDCAxyPt(NULL),
  fHistTrackP(NULL),
  fHistTrackPt(NULL),
  fHistTrackEta(NULL),    

  fHistV0PV0DecayLengthXY(NULL),
  fHistV0PV0DecayLength(NULL),
  fHistV0PV0PointingAngleXY(NULL),
  fHistV0PV0PointingAngle(NULL),
  fHistV0PV0DCAXY(NULL),
  fHistV0PV0DCA(NULL),
  fHistV0PV0TrackDistanceXY(NULL),
  fHistV0PV0TrackDistance(NULL),
  fHistV0PV0Chi2perNDF(NULL),
  fHistV0PV0PropLifeTime(NULL),
  fHistV0TrackDCAXY(NULL),
  fHistV0TrackDCA(NULL),

  fHistSelV0PV0DecayLengthXY(NULL),
  fHistSelV0PV0DecayLength(NULL),
  fHistSelV0PV0PointingAngleXY(NULL),
  fHistSelV0PV0PointingAngle(NULL),
  fHistSelV0PV0DCAXY(NULL),
  fHistSelV0PV0DCA(NULL),
  fHistSelV0PV0TrackDistanceXY(NULL),
  fHistSelV0PV0TrackDistance(NULL),
  fHistSelV0PV0Chi2perNDF(NULL),
  fHistSelV0PV0PropLifeTime(NULL),
  fHistSelV0TrackDCAXY(NULL),
  fHistSelV0TrackDCA(NULL),

  fHistV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistV0PV0DecayLength_TrueK0s(NULL),
  fHistV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistV0PV0PointingAngle_TrueK0s(NULL),
  fHistV0PV0DCAXY_TrueK0s(NULL),
  fHistV0PV0DCA_TrueK0s(NULL),
  fHistV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistV0PV0TrackDistance_TrueK0s(NULL),
  fHistV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistV0PV0PropLifeTime_TrueK0s(NULL),
  fHistV0TrackDCAXY_TrueK0s(NULL),
  fHistV0TrackDCA_TrueK0s(NULL),

  fHistSelV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistSelV0PV0DecayLength_TrueK0s(NULL),
  fHistSelV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistSelV0PV0PointingAngle_TrueK0s(NULL),
  fHistSelV0PV0DCAXY_TrueK0s(NULL),
  fHistSelV0PV0DCA_TrueK0s(NULL),
  fHistSelV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistSelV0PV0TrackDistance_TrueK0s(NULL),
  fHistSelV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistSelV0PV0PropLifeTime_TrueK0s(NULL),
  fHistSelV0TrackDCAXY_TrueK0s(NULL),
  fHistSelV0TrackDCA_TrueK0s(NULL),
  
  fHistLSV0PV0DecayLengthXY(NULL),
  fHistLSV0PV0DecayLength(NULL),
  fHistLSV0PV0PointingAngleXY(NULL),
  fHistLSV0PV0PointingAngle(NULL),
  fHistLSV0PV0DCAXY(NULL),
  fHistLSV0PV0DCA(NULL),
  fHistLSV0PV0TrackDistanceXY(NULL),
  fHistLSV0PV0TrackDistance(NULL),
  fHistLSV0PV0Chi2perNDF(NULL),
  fHistLSV0PV0PropLifeTime(NULL),
  fHistLSV0TrackDCAXY(NULL),
  fHistLSV0TrackDCA(NULL),

  fHistLSSelV0PV0DecayLengthXY(NULL),
  fHistLSSelV0PV0DecayLength(NULL),
  fHistLSSelV0PV0PointingAngleXY(NULL),
  fHistLSSelV0PV0PointingAngle(NULL),
  fHistLSSelV0PV0DCAXY(NULL),
  fHistLSSelV0PV0DCA(NULL),
  fHistLSSelV0PV0TrackDistanceXY(NULL),
  fHistLSSelV0PV0TrackDistance(NULL),
  fHistLSSelV0PV0Chi2perNDF(NULL),
  fHistLSSelV0PV0PropLifeTime(NULL),
  fHistLSSelV0TrackDCAXY(NULL),
  fHistLSSelV0TrackDCA(NULL),

  fHistLSV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistLSV0PV0DecayLength_TrueK0s(NULL),
  fHistLSV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistLSV0PV0PointingAngle_TrueK0s(NULL),
  fHistLSV0PV0DCAXY_TrueK0s(NULL),
  fHistLSV0PV0DCA_TrueK0s(NULL),
  fHistLSV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistLSV0PV0TrackDistance_TrueK0s(NULL),
  fHistLSV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistLSV0PV0PropLifeTime_TrueK0s(NULL),
  fHistLSV0TrackDCAXY_TrueK0s(NULL),
  fHistLSV0TrackDCA_TrueK0s(NULL),

  fHistLSSelV0PV0DecayLengthXY_TrueK0s(NULL),
  fHistLSSelV0PV0DecayLength_TrueK0s(NULL),
  fHistLSSelV0PV0PointingAngleXY_TrueK0s(NULL),
  fHistLSSelV0PV0PointingAngle_TrueK0s(NULL),
  fHistLSSelV0PV0DCAXY_TrueK0s(NULL),
  fHistLSSelV0PV0DCA_TrueK0s(NULL),
  fHistLSSelV0PV0TrackDistanceXY_TrueK0s(NULL),
  fHistLSSelV0PV0TrackDistance_TrueK0s(NULL),
  fHistLSSelV0PV0Chi2perNDF_TrueK0s(NULL),
  fHistLSSelV0PV0PropLifeTime_TrueK0s(NULL),
  fHistLSSelV0TrackDCAXY_TrueK0s(NULL),
  fHistLSSelV0TrackDCA_TrueK0s(NULL),

  all(0),
  abord(0)
  
{

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
   double fVtxBins[]  = {-50, -10.5, -6, -2, 0, 2, 6, 10.5, 50};
   double fPsiBins[]  = {-10, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 10};

   int fNCentBins = sizeof(fCentBins) / sizeof(double) - 1;
   int fNVtxZBins = sizeof(fVtxBins) / sizeof(double) - 1;
   int fNPsiBins  = sizeof(fPsiBins) / sizeof(double) - 1;

   fPoolMuonTrackMgrK0s = new AliEventPoolManager(fPoolSize, fTrackDepth,
					       fNCentBins, (double *)fCentBins,
					       fNVtxZBins, (double *)fVtxBins,
					       fNPsiBins, (double *)fPsiBins);
   fPoolMuonTrackMgrK0s -> SetTargetValues(fTrackDepth, (double)fReadyFraction,fPoolSize);

   fPoolMuonTrackMgrPion = new AliEventPoolManager(fPoolSize, fTrackDepth,
						   fNCentBins, (double *)fCentBins,
						   fNVtxZBins, (double *)fVtxBins,
						   fNPsiBins, (double *)fPsiBins);
   fPoolMuonTrackMgrPion -> SetTargetValues(fTrackDepth, (double)fReadyFraction,fPoolSize);
   
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

   double min_angle = 0;
   double max_angle = 180.;
   double width_angle = 60.;

   vector<double> bins_cent{0,1,5,10,20,40,60,80,100};
   int binnum_cent = bins_cent.size() - 1;
   
   int bins[4] = {int((max_mass - min_mass) / width_mass), int((max_pt - min_pt) / width_pt), binnum_cent,int((max_angle-min_angle)/width_angle)};
   double min_bins[4] = {min_mass, min_pt, 0, min_angle};
   double max_bins[4] = {max_mass, max_pt, 100, max_angle};   

   fSparseULSPionPair  = new THnSparseF("fSparseULSPionPair", "", 4, bins,min_bins, max_bins);
   fSparseLSppPionPair = new THnSparseF("fSparseLSppPionPair", "", 4,bins, min_bins, max_bins);
   fSparseLSmmPionPair = new THnSparseF("fSparseLSmmPionPair", "", 4,bins, min_bins, max_bins);
   fSparseULSPionPair  -> SetBinEdges(2,bins_cent.data());
   fSparseLSppPionPair -> SetBinEdges(2,bins_cent.data());
   fSparseLSmmPionPair -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseULSPionPair);
   fOutputList -> Add(fSparseLSppPionPair);
   fOutputList -> Add(fSparseLSmmPionPair);

   fSparseMixULSPionPair  = new THnSparseF("fSparseMixULSPionPair", "", 4, bins,min_bins, max_bins);
   fSparseMixLSppPionPair = new THnSparseF("fSparseMixLSppPionPair", "", 4,bins, min_bins, max_bins);
   fSparseMixLSmmPionPair = new THnSparseF("fSparseMixLSmmPionPair", "", 4,bins, min_bins, max_bins);
   fSparseMixULSPionPair  -> SetBinEdges(2,bins_cent.data());
   fSparseMixLSppPionPair -> SetBinEdges(2,bins_cent.data());
   fSparseMixLSmmPionPair -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseMixULSPionPair);
   fOutputList -> Add(fSparseMixLSppPionPair);
   fOutputList -> Add(fSparseMixLSmmPionPair);

   fSparseULSPionPairBeforeCuts  = new THnSparseF("fSparseULSPionPairBeforeCuts", "", 4, bins,min_bins, max_bins);
   fSparseLSppPionPairBeforeCuts = new THnSparseF("fSparseLSppPionPairBeforeCuts", "", 4,bins, min_bins, max_bins);
   fSparseLSmmPionPairBeforeCuts = new THnSparseF("fSparseLSmmPionPairBeforeCuts", "", 4,bins, min_bins, max_bins);
   fSparseULSPionPairBeforeCuts  -> SetBinEdges(2,bins_cent.data());
   fSparseLSppPionPairBeforeCuts -> SetBinEdges(2,bins_cent.data());
   fSparseLSmmPionPairBeforeCuts -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseULSPionPairBeforeCuts);
   fOutputList -> Add(fSparseLSppPionPairBeforeCuts);
   fOutputList -> Add(fSparseLSmmPionPairBeforeCuts);
   
   fSparseULSPionPair_PassChi2perNDFCut         = new THnSparseF("fSparseULSPionPair_PassChi2perNDFCut", "", 4, bins,min_bins, max_bins);
   fSparseULSPionPair_PassCPVXYCut              = new THnSparseF("fSparseULSPionPair_PassCPVXYCut", "", 4, bins,min_bins, max_bins);
   fSparseULSPionPair_PassDecayLengthXYCut      = new THnSparseF("fSparseULSPionPair_PassDecayLengthXYCut", "", 4, bins,min_bins, max_bins);
   fSparseULSPionPair_PassDaughterDistanceXYCut = new THnSparseF("fSparseULSPionPair_PassDaughterDistanceXYCut", "", 4, bins,min_bins, max_bins);  
   fSparseULSPionPair_PassDCAXYCut              = new THnSparseF("fSparseULSPionPair_PassDCAXYCut", "", 4, bins,min_bins, max_bins);
   fSparseULSPionPair_PassLifetimeCut           = new THnSparseF("fSparseULSPionPair_PassLifetimeCut", "", 4, bins,min_bins, max_bins);
   fOutputList -> Add(fSparseULSPionPair_PassChi2perNDFCut);
   fOutputList -> Add(fSparseULSPionPair_PassCPVXYCut);
   fOutputList -> Add(fSparseULSPionPair_PassDecayLengthXYCut);
   fOutputList -> Add(fSparseULSPionPair_PassDaughterDistanceXYCut);
   fOutputList -> Add(fSparseULSPionPair_PassDCAXYCut);
   fOutputList -> Add(fSparseULSPionPair_PassLifetimeCut);


   fSparseNeutralK0sPair          = new THnSparseF("fSparseNeutralK0sPair", "", 4, bins,min_bins, max_bins);
   fSparseNeutralNegativeK0sPair  = new THnSparseF("fSparseNeutralNegativeK0sPair", "", 4, bins,min_bins, max_bins);
   fSparseNeutralPositiveK0sPair  = new THnSparseF("fSparseNeutralPositiveK0sPair", "", 4, bins,min_bins, max_bins);
   fSparsePositiveK0sPair         = new THnSparseF("fSparsePositiveK0sPair", "", 4, bins,min_bins, max_bins);
   fSparseNegativeK0sPair         = new THnSparseF("fSparseNegativeK0sPair", "", 4, bins,min_bins, max_bins);
   fSparsePositiveNegativeK0sPair = new THnSparseF("fSparsePositiveNegativeK0sPair", "", 4, bins,min_bins, max_bins);
   fSparseNeutralK0sPair         -> SetBinEdges(2,bins_cent.data());
   fSparseNeutralNegativeK0sPair -> SetBinEdges(2,bins_cent.data());
   fSparseNeutralPositiveK0sPair -> SetBinEdges(2,bins_cent.data());
   fSparsePositiveK0sPair        -> SetBinEdges(2,bins_cent.data());
   fSparseNegativeK0sPair        -> SetBinEdges(2,bins_cent.data());
   fSparsePositiveNegativeK0sPair->SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseNeutralK0sPair);
   fOutputList -> Add(fSparseNeutralNegativeK0sPair);
   fOutputList -> Add(fSparseNeutralPositiveK0sPair);
   fOutputList -> Add(fSparsePositiveK0sPair);
   fOutputList -> Add(fSparseNegativeK0sPair);
   fOutputList -> Add(fSparsePositiveNegativeK0sPair);
      
   fSparseTrueK0s = new THnSparseF("fSparseTrueK0s", "", 4, bins,min_bins, max_bins);
   fSparseTrueK0s -> SetBinEdges(2,bins_cent.data());
   fOutputList    -> Add(fSparseTrueK0s);
   
   int bins_K0sParams[4]        = {10, 400, 400, 600};
   double min_bins_K0sParams[4] = {0,    0,   0,  -1};
   double max_bins_K0sParams[4] = {10,   2,   2,   2};
   
   fSparseTrueK0sRecK0s = new THnSparseF("fSparseTrueK0sRecK0s", "", 4, bins_K0sParams, min_bins_K0sParams, max_bins_K0sParams);
   fOutputList -> Add(fSparseTrueK0sRecK0s);
   
   width_pt = 0.1;
   
   fHistV0PV0DecayLengthXY   = new TH2F("fHistV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistV0PV0DecayLength     = new TH2F("fHistV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistV0PV0PointingAngleXY = new TH2F("fHistV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistV0PV0PointingAngle   = new TH2F("fHistV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistV0PV0DCAXY           = new TH2F("fHistV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistV0PV0DCA             = new TH2F("fHistV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistV0PV0TrackDistanceXY = new TH2F("fHistV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistV0PV0TrackDistance   = new TH2F("fHistV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistV0PV0Chi2perNDF      = new TH2F("fHistV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistV0PV0PropLifeTime    = new TH2F("fHistV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
   fHistV0TrackDCAXY         = new TH2F("fHistV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fHistV0TrackDCA           = new TH2F("fHistV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fOutputList -> Add(fHistV0PV0DecayLengthXY);
   fOutputList -> Add(fHistV0PV0DecayLength);
   fOutputList -> Add(fHistV0PV0PointingAngleXY);
   fOutputList -> Add(fHistV0PV0PointingAngle);
   fOutputList -> Add(fHistV0PV0DCAXY);
   fOutputList -> Add(fHistV0PV0DCA);
   fOutputList -> Add(fHistV0PV0TrackDistanceXY);
   fOutputList -> Add(fHistV0PV0TrackDistance);
   fOutputList -> Add(fHistV0PV0Chi2perNDF);
   fOutputList -> Add(fHistV0PV0PropLifeTime);
   fOutputList -> Add(fHistV0TrackDCAXY);
   fOutputList -> Add(fHistV0TrackDCA);

   fHistSelV0PV0DecayLengthXY   = new TH2F("fHistSelV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistSelV0PV0DecayLength     = new TH2F("fHistSelV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistSelV0PV0PointingAngleXY = new TH2F("fHistSelV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistSelV0PV0PointingAngle   = new TH2F("fHistSelV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistSelV0PV0DCAXY           = new TH2F("fHistSelV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistSelV0PV0DCA             = new TH2F("fHistSelV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistSelV0PV0TrackDistanceXY = new TH2F("fHistSelV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistSelV0PV0TrackDistance   = new TH2F("fHistSelV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistSelV0PV0Chi2perNDF      = new TH2F("fHistSelV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistSelV0PV0PropLifeTime    = new TH2F("fHistSelV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
   fHistSelV0TrackDCAXY         = new TH2F("fHistSelV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fHistSelV0TrackDCA           = new TH2F("fHistSelV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fOutputList -> Add(fHistSelV0PV0DecayLengthXY);
   fOutputList -> Add(fHistSelV0PV0DecayLength);
   fOutputList -> Add(fHistSelV0PV0PointingAngleXY);
   fOutputList -> Add(fHistSelV0PV0PointingAngle);
   fOutputList -> Add(fHistSelV0PV0DCAXY);
   fOutputList -> Add(fHistSelV0PV0DCA);
   fOutputList -> Add(fHistSelV0PV0TrackDistanceXY);
   fOutputList -> Add(fHistSelV0PV0TrackDistance);
   fOutputList -> Add(fHistSelV0PV0Chi2perNDF);
   fOutputList -> Add(fHistSelV0PV0PropLifeTime);
   fOutputList -> Add(fHistSelV0TrackDCAXY);
   fOutputList -> Add(fHistSelV0TrackDCA);

   fHistLSV0PV0DecayLengthXY   = new TH2F("fHistLSV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistLSV0PV0DecayLength     = new TH2F("fHistLSV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistLSV0PV0PointingAngleXY = new TH2F("fHistLSV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistLSV0PV0PointingAngle   = new TH2F("fHistLSV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistLSV0PV0DCAXY           = new TH2F("fHistLSV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistLSV0PV0DCA             = new TH2F("fHistLSV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistLSV0PV0TrackDistanceXY = new TH2F("fHistLSV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistLSV0PV0TrackDistance   = new TH2F("fHistLSV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistLSV0PV0Chi2perNDF      = new TH2F("fHistLSV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistLSV0PV0PropLifeTime    = new TH2F("fHistLSV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
   fHistLSV0TrackDCAXY         = new TH2F("fHistLSV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fHistLSV0TrackDCA           = new TH2F("fHistLSV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fOutputList -> Add(fHistLSV0PV0DecayLengthXY);
   fOutputList -> Add(fHistLSV0PV0DecayLength);
   fOutputList -> Add(fHistLSV0PV0PointingAngleXY);
   fOutputList -> Add(fHistLSV0PV0PointingAngle);
   fOutputList -> Add(fHistLSV0PV0DCAXY);
   fOutputList -> Add(fHistLSV0PV0DCA);
   fOutputList -> Add(fHistLSV0PV0TrackDistanceXY);
   fOutputList -> Add(fHistLSV0PV0TrackDistance);
   fOutputList -> Add(fHistLSV0PV0Chi2perNDF);
   fOutputList -> Add(fHistLSV0PV0PropLifeTime);
   fOutputList -> Add(fHistLSV0TrackDCAXY);
   fOutputList -> Add(fHistLSV0TrackDCA);

   fHistLSSelV0PV0DecayLengthXY   = new TH2F("fHistLSSelV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistLSSelV0PV0DecayLength     = new TH2F("fHistLSSelV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
   fHistLSSelV0PV0PointingAngleXY = new TH2F("fHistLSSelV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistLSSelV0PV0PointingAngle   = new TH2F("fHistLSSelV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
   fHistLSSelV0PV0DCAXY           = new TH2F("fHistLSSelV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistLSSelV0PV0DCA             = new TH2F("fHistLSSelV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
   fHistLSSelV0PV0TrackDistanceXY = new TH2F("fHistLSSelV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistLSSelV0PV0TrackDistance   = new TH2F("fHistLSSelV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistLSSelV0PV0Chi2perNDF      = new TH2F("fHistLSSelV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistLSSelV0PV0PropLifeTime    = new TH2F("fHistLSSelV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
   fHistLSSelV0TrackDCAXY         = new TH2F("fHistLSSelV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fHistLSSelV0TrackDCA           = new TH2F("fHistLSSelV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
   fOutputList -> Add(fHistLSSelV0PV0DecayLengthXY);
   fOutputList -> Add(fHistLSSelV0PV0DecayLength);
   fOutputList -> Add(fHistLSSelV0PV0PointingAngleXY);
   fOutputList -> Add(fHistLSSelV0PV0PointingAngle);
   fOutputList -> Add(fHistLSSelV0PV0DCAXY);
   fOutputList -> Add(fHistLSSelV0PV0DCA);
   fOutputList -> Add(fHistLSSelV0PV0TrackDistanceXY);
   fOutputList -> Add(fHistLSSelV0PV0TrackDistance);
   fOutputList -> Add(fHistLSSelV0PV0Chi2perNDF);
   fOutputList -> Add(fHistLSSelV0PV0PropLifeTime);
   fOutputList -> Add(fHistLSSelV0TrackDCAXY);
   fOutputList -> Add(fHistLSSelV0TrackDCA);

   if (fIsMC) {
     fHistV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistV0PV0DecayLength_TrueK0s     = new TH2F("fHistV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistV0PV0PointingAngle_TrueK0s   = new TH2F("fHistV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistV0PV0DCAXY_TrueK0s           = new TH2F("fHistV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistV0PV0DCA_TrueK0s             = new TH2F("fHistV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistV0PV0TrackDistance_TrueK0s   = new TH2F("fHistV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
     fHistV0TrackDCAXY_TrueK0s         = new TH2F("fHistV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fHistV0TrackDCA_TrueK0s           = new TH2F("fHistV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fOutputList -> Add(fHistV0PV0DecayLengthXY_TrueK0s);
     fOutputList -> Add(fHistV0PV0DecayLength_TrueK0s);
     fOutputList -> Add(fHistV0PV0PointingAngleXY_TrueK0s);
     fOutputList -> Add(fHistV0PV0PointingAngle_TrueK0s);
     fOutputList -> Add(fHistV0PV0DCAXY_TrueK0s);
     fOutputList -> Add(fHistV0PV0DCA_TrueK0s);
     fOutputList -> Add(fHistV0PV0TrackDistanceXY_TrueK0s);
     fOutputList -> Add(fHistV0PV0TrackDistance_TrueK0s);
     fOutputList -> Add(fHistV0PV0Chi2perNDF_TrueK0s);
     fOutputList -> Add(fHistV0PV0PropLifeTime_TrueK0s);
     fOutputList -> Add(fHistV0TrackDCAXY_TrueK0s);
     fOutputList -> Add(fHistV0TrackDCA_TrueK0s);

     fHistSelV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistSelV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistSelV0PV0DecayLength_TrueK0s     = new TH2F("fHistSelV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistSelV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistSelV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistSelV0PV0PointingAngle_TrueK0s   = new TH2F("fHistSelV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistSelV0PV0DCAXY_TrueK0s           = new TH2F("fHistSelV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistSelV0PV0DCA_TrueK0s             = new TH2F("fHistSelV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistSelV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistSelV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistSelV0PV0TrackDistance_TrueK0s   = new TH2F("fHistSelV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistSelV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistSelV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistSelV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistSelV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
     fHistSelV0TrackDCAXY_TrueK0s         = new TH2F("fHistSelV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fHistSelV0TrackDCA_TrueK0s           = new TH2F("fHistSelV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fOutputList -> Add(fHistSelV0PV0DecayLengthXY_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0DecayLength_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0PointingAngleXY_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0PointingAngle_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0DCAXY_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0DCA_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0TrackDistanceXY_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0TrackDistance_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0Chi2perNDF_TrueK0s);
     fOutputList -> Add(fHistSelV0PV0PropLifeTime_TrueK0s);
     fOutputList -> Add(fHistSelV0TrackDCAXY_TrueK0s);
     fOutputList -> Add(fHistSelV0TrackDCA_TrueK0s);


     fHistLSV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistLSV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistLSV0PV0DecayLength_TrueK0s     = new TH2F("fHistLSV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistLSV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistLSV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistLSV0PV0PointingAngle_TrueK0s   = new TH2F("fHistLSV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistLSV0PV0DCAXY_TrueK0s           = new TH2F("fHistLSV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistLSV0PV0DCA_TrueK0s             = new TH2F("fHistLSV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistLSV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistLSV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistLSV0PV0TrackDistance_TrueK0s   = new TH2F("fHistLSV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistLSV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistLSV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistLSV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistLSV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
     fHistLSV0TrackDCAXY_TrueK0s         = new TH2F("fHistLSV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fHistLSV0TrackDCA_TrueK0s           = new TH2F("fHistLSV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fOutputList -> Add(fHistLSV0PV0DecayLengthXY_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0DecayLength_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0PointingAngleXY_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0PointingAngle_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0DCAXY_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0DCA_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0TrackDistanceXY_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0TrackDistance_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0Chi2perNDF_TrueK0s);
     fOutputList -> Add(fHistLSV0PV0PropLifeTime_TrueK0s);
     fOutputList -> Add(fHistLSV0TrackDCAXY_TrueK0s);
     fOutputList -> Add(fHistLSV0TrackDCA_TrueK0s);

     fHistLSSelV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistLSSelV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistLSSelV0PV0DecayLength_TrueK0s     = new TH2F("fHistLSSelV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 500, 0., 100.);
     fHistLSSelV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistLSSelV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistLSSelV0PV0PointingAngle_TrueK0s   = new TH2F("fHistLSSelV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0.95, 1.0);
     fHistLSSelV0PV0DCAXY_TrueK0s           = new TH2F("fHistLSSelV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistLSSelV0PV0DCA_TrueK0s             = new TH2F("fHistLSSelV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 400, 0., 80.);
     fHistLSSelV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistLSSelV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistLSSelV0PV0TrackDistance_TrueK0s   = new TH2F("fHistLSSelV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistLSSelV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistLSSelV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistLSSelV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistLSSelV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,250, 0, 50.);
     fHistLSSelV0TrackDCAXY_TrueK0s         = new TH2F("fHistLSSelV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fHistLSSelV0TrackDCA_TrueK0s           = new TH2F("fHistLSSelV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 1000, 0., 100.);
     fOutputList -> Add(fHistLSSelV0PV0DecayLengthXY_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0DecayLength_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0PointingAngleXY_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0PointingAngle_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0DCAXY_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0DCA_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0TrackDistanceXY_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0TrackDistance_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0Chi2perNDF_TrueK0s);
     fOutputList -> Add(fHistLSSelV0PV0PropLifeTime_TrueK0s);
     fOutputList -> Add(fHistLSSelV0TrackDCAXY_TrueK0s);
     fOutputList -> Add(fHistLSSelV0TrackDCA_TrueK0s);
   }
   
   fHistTrackP                    = new TH1F("fHistTrackP", "", 100, 0, 10);
   fHistTrackPt                   = new TH1F("fHistTrackPt", "", 100, 0, 10);
   fHistTrackEta                  = new TH1F("fHistTrackEta", "", 20, -1, 1);   
   fHistTPCNClusts                = new TH1F("fHistTPCNClusts", "", 160, 0, 160);
   fHistSPDNClusts                = new TH1F("fHistSPDNClusts", "", 6, 0, 6);
   fHistTPCCrossRowsFindableRatio = new TH1F("fHistTPCCrossRowsFindableRatio", "", 200, 0, 2);
   fHistReducedChi2TPC            = new TH1F("fHistReducedChi2TPC", "", 100, 0, 10);
   fHistReducedChi2ITS            = new TH1F("fHistReducedChi2ITS", "", 250, 0, 50);
   fHistDCAz                      = new TH1F("fHistDCAz", "", 100, 0, 10);
   fHistDCAxyPt                   = new TH2F("fHistDCAxyPt", "", 200, -1, 1, 100, 0, 10);
   
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

   fHistTPCdEdxP         = new TH2F("fHistTPCdEdxP", "", (max_p-min_p)/width_p,min_p,max_p,(max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
   fHistBetaP            = new TH2F("fHistBetaP", "", (max_p-min_p)/width_p,min_p,max_p,(max_beta - min_beta) / width_beta, min_beta, max_beta);
   fHistTPCSigmaElectron = new TH2F("fHistTPCSigmaElectron", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTOFSigmaElectron = new TH2F("fHistTOFSigmaElectron", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTPCSigmaMuon     = new TH2F("fHistTPCSigmaMuon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTOFSigmaMuon     = new TH2F("fHistTOFSigmaMuon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTPCSigmaPion     = new TH2F("fHistTPCSigmaPion", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTOFSigmaPion     = new TH2F("fHistTOFSigmaPion", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTPCSigmaKaon     = new TH2F("fHistTPCSigmaKaon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTOFSigmaKaon     = new TH2F("fHistTOFSigmaKaon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTPCSigmaProton   = new TH2F("fHistTPCSigmaProton", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistTOFSigmaProton   = new TH2F("fHistTOFSigmaProton", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fOutputList -> Add(fHistTPCdEdxP);
   fOutputList -> Add(fHistBetaP);
   fOutputList -> Add(fHistTPCSigmaElectron);
   fOutputList -> Add(fHistTOFSigmaElectron);
   fOutputList -> Add(fHistTPCSigmaMuon);
   fOutputList -> Add(fHistTOFSigmaMuon);
   fOutputList -> Add(fHistTPCSigmaPion);
   fOutputList -> Add(fHistTOFSigmaPion);
   fOutputList -> Add(fHistTPCSigmaKaon);
   fOutputList -> Add(fHistTOFSigmaKaon);
   fOutputList -> Add(fHistTPCSigmaProton);
   fOutputList -> Add(fHistTOFSigmaProton);

   fHistSelTPCdEdxP         = new TH2F("fHistSelTPCdEdxP", "", (max_p-min_p)/width_p,min_p,max_p,(max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
   fHistSelBetaP            = new TH2F("fHistSelBetaP", "", (max_p-min_p)/width_p,min_p,max_p,(max_beta - min_beta) / width_beta, min_beta, max_beta);
   fHistSelTPCSigmaElectron = new TH2F("fHistSelTPCSigmaElectron", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTOFSigmaElectron = new TH2F("fHistSelTOFSigmaElectron", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTPCSigmaMuon     = new TH2F("fHistSelTPCSigmaMuon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTOFSigmaMuon     = new TH2F("fHistSelTOFSigmaMuon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTPCSigmaPion     = new TH2F("fHistSelTPCSigmaPion", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTOFSigmaPion     = new TH2F("fHistSelTOFSigmaPion", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTPCSigmaKaon     = new TH2F("fHistSelTPCSigmaKaon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTOFSigmaKaon     = new TH2F("fHistSelTOFSigmaKaon", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTPCSigmaProton   = new TH2F("fHistSelTPCSigmaProton", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fHistSelTOFSigmaProton   = new TH2F("fHistSelTOFSigmaProton", "", (max_p-min_p)/width_p,min_p,max_p,(max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
   fOutputList -> Add(fHistSelTPCdEdxP);
   fOutputList -> Add(fHistSelBetaP);
   fOutputList -> Add(fHistSelTPCSigmaElectron);
   fOutputList -> Add(fHistSelTOFSigmaElectron);
   fOutputList -> Add(fHistSelTPCSigmaMuon);
   fOutputList -> Add(fHistSelTOFSigmaMuon);
   fOutputList -> Add(fHistSelTPCSigmaPion);
   fOutputList -> Add(fHistSelTOFSigmaPion);
   fOutputList -> Add(fHistSelTPCSigmaKaon);
   fOutputList -> Add(fHistSelTOFSigmaKaon);
   fOutputList -> Add(fHistSelTPCSigmaProton);
   fOutputList -> Add(fHistSelTOFSigmaProton);

   fHistEventVtxZ    = new TH1F("fHistEventVtxZ", "", 60, -30, 30);
   fHistEventCent    = new TH1F("fHistEventCent", "", 100, 0, 100);
   fHistEventMulti   = new TH1F("fHistEventMulti", "", 200, 0, 200);
   fHistEventVtxCont = new TH1F("fHistEventVtxCont", "", 100, 0, 100);
   fOutputList -> Add(fHistEventVtxZ);
   fOutputList -> Add(fHistEventCent);
   fOutputList -> Add(fHistEventMulti);
   fOutputList -> Add(fHistEventVtxCont);

   PostData(1, fOutputList);
 }

 //________________________________________________________________________

 void AliAnalysisTaskAODTrackPair::UserExec(Option_t *) {
   
   fArrayK0s = new TObjArray();
   fArrayK0s -> SetOwner(true);
   
   if (!Initialize())            return ;
   if (!fUtils->isAcceptEvent()) return ;
   if (!fIsMixingAnalysis) {
     AddK0sArray(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));
     K0sPairAnalysis(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));
   } else {
     //K0sPairAnalysisEventMixing(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));
   }
   if(fIsMC){
     ProcessMC();
   }
   if (fArrayK0s) delete fArrayK0s;
   
   cout<<abord<<"  outof "<<all<<endl;
   
 }

 bool AliAnalysisTaskAODTrackPair::Initialize() {
   
   fEvent = dynamic_cast<AliAODEvent *>(InputEvent());
   
   if (!fUtils->setEvent(fEvent, fInputHandler)) return false;
   
   fPrimVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();
   
   if (!fPrimVtx)                                       return false;   
   if (fPrimVtx->GetNContributors() < fMinNContPrimVtx) return false;

   fPrimVtx->GetXYZ(fPrimVtxPos);
   
   fMultSelection = (AliMultSelection *)fEvent->FindListObject("MultSelection");
   
   if (!fMultSelection) return false;

   fCent = fMultSelection->GetMultiplicityPercentile(fMethodCent.c_str(), false);
   
   fHistEventVtxZ    -> Fill(fPrimVtxPos[2]);
   fHistEventCent    -> Fill(fCent);
   fHistEventMulti   -> Fill(fUtils->getNCorrSPDTrkInfo(1));
   fHistEventVtxCont -> Fill(fPrimVtx->GetNContributors());
   
   if ( fIsMC ) {
     
     fMCTrackArray = dynamic_cast<TClonesArray *>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
     
     if (!fMCTrackArray) return false;     
     
     fUtils->setMCArray(fMCTrackArray);
     
   }   
   return true;
 }

bool AliAnalysisTaskAODTrackPair::TrackPIDChecker(AliAODTrack *track, AliPID::EParticleType pid, bool isSel) {
  
  double p      = track->P();
  double sigTOF = track->GetTOFsignal();
  double length = track->GetIntegratedLength();
  double beta   = sigTOF > 0 ? length / (2.99792457999999984e-02 * sigTOF) : -999;
  double dEdx   = track->GetTPCsignal();
  
  if (isSel) {
    fHistSelTPCdEdxP         -> Fill(p, dEdx);
    fHistSelTPCSigmaElectron -> Fill(p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistSelTPCSigmaMuon     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistSelTPCSigmaPion     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistSelTPCSigmaKaon     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistSelTPCSigmaProton   -> Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistSelBetaP            -> Fill(p, beta);
      fHistSelTOFSigmaElectron -> Fill(p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistSelTOFSigmaMuon     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistSelTOFSigmaPion     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistSelTOFSigmaKaon     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistSelTOFSigmaProton   -> Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
    }
  } else {
    fHistTPCdEdxP         -> Fill(p, dEdx);
    fHistTPCSigmaElectron -> Fill(p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistTPCSigmaMuon     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistTPCSigmaPion     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistTPCSigmaKaon     -> Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistTPCSigmaProton   -> Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistBetaP            -> Fill(p, beta);
      fHistTOFSigmaElectron -> Fill(p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistTOFSigmaMuon     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistTOFSigmaPion     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistTOFSigmaKaon     -> Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistTOFSigmaProton   -> Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::TrackQualityChecker(AliAODTrack *track) {

  fHistTPCNClusts->Fill(track->GetTPCNcls());

  int nSPD = 0;

  if (track->HasPointOnITSLayer(0)) ++nSPD;
  if (track->HasPointOnITSLayer(1)) ++nSPD;

  fHistSPDNClusts->Fill(nSPD);
  
  double fr       = track->GetTPCCrossedRows() > 0 ? track->GetTPCNclsF() / track->GetTPCCrossedRows() : 0;
  double rchi2tpc = track->GetTPCNcls() > 0 ? track->GetTPCchi2() / track->GetTPCNcls() : 0;
  double rchi2its = track->GetITSNcls() > 0 ? track->GetITSchi2() / track->GetITSNcls() : 0;
  
  fHistTPCCrossRowsFindableRatio -> Fill(fr);
  fHistReducedChi2TPC            -> Fill(rchi2tpc);
  fHistReducedChi2ITS            -> Fill(rchi2its);
  
  double dca_xy, dca_z;

  fHistDCAz    -> Fill(dca_z);
  fHistDCAxyPt -> Fill(dca_xy, track->Pt());  
  fHistTrackP   -> Fill(track->P());
  fHistTrackPt  -> Fill(track->Pt());
  fHistTrackEta -> Fill(track->Eta());

  return true;
}

bool AliAnalysisTaskAODTrackPair::V0QualityChecker(K0sContainer *v0, bool isSel) {
    
  double p  = v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1] + v0->P[2]*v0->P[2] > 0 ? 
    sqrt(v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1] + v0->P[2]*v0->P[2]) : 0;
  
  double DCApairXY = v0->dDcaXY[0]*v0->dDcaXY[0] + v0->dDcaXY[1]*v0->dDcaXY[1] > 0 ?
    sqrt(v0->dDcaXY[0]*v0->dDcaXY[0] + v0->dDcaXY[1]*v0->dDcaXY[1]) : 999;
  double DCApair   = v0->dDcaXYZ[0]*v0->dDcaXYZ[0] + v0->dDcaXYZ[1]*v0->dDcaXYZ[1] > 0 ?
    sqrt(v0->dDcaXYZ[0]*v0->dDcaXYZ[0] + v0->dDcaXYZ[1]*v0->dDcaXYZ[1]) : 999;

  if (isSel) {    
    
    fHistSelV0PV0DecayLengthXY   -> Fill(p, v0->DecayLengthXY);
    fHistSelV0PV0DecayLength     -> Fill(p, v0->DecayLengthXYZ);
    fHistSelV0PV0PointingAngleXY -> Fill(p, v0->CpvXY);
    fHistSelV0PV0PointingAngle   -> Fill(p, v0->CpvXYZ);
    fHistSelV0PV0DCAXY           -> Fill(p, v0->DcaXY);
    fHistSelV0PV0DCA             -> Fill(p, v0->DcaXYZ);
    fHistSelV0PV0TrackDistanceXY -> Fill(p, v0->DcaDaughtersXY);
    fHistSelV0PV0TrackDistance    -> Fill(p, v0->DcaDaughtersXYZ);
    fHistSelV0PV0Chi2perNDF      -> Fill(p, v0->rChi2V0);
    fHistSelV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
    fHistSelV0TrackDCAXY         -> Fill(p, DCApairXY);
    fHistSelV0TrackDCA           -> Fill(p, DCApair);

    if (fabs(v0->Charge) > 0) {
      fHistLSSelV0PV0DecayLengthXY   -> Fill(p, v0->DecayLengthXY);
      fHistLSSelV0PV0DecayLength     -> Fill(p, v0->DecayLengthXYZ);
      fHistLSSelV0PV0PointingAngleXY -> Fill(p, v0->CpvXY);
      fHistLSSelV0PV0PointingAngle   -> Fill(p, v0->CpvXYZ);
      fHistLSSelV0PV0DCAXY           -> Fill(p, v0->DcaXY);
      fHistLSSelV0PV0DCA             -> Fill(p, v0->DcaXYZ);
      fHistLSSelV0PV0TrackDistanceXY -> Fill(p, v0->DcaDaughtersXY);
      fHistLSSelV0PV0TrackDistance    -> Fill(p, v0->DcaDaughtersXYZ);
      fHistLSSelV0PV0Chi2perNDF      -> Fill(p, v0->rChi2V0);
      fHistLSSelV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
      fHistLSSelV0TrackDCAXY         -> Fill(p, DCApairXY);
      fHistLSSelV0TrackDCA           -> Fill(p, DCApair);
    }

    if(fIsMC && 
       fUtils->isSameMotherPair(v0->dLabel[0],v0->dLabel[1]) &&
       fUtils->getMotherPdgCode(v0->dLabel[0]) == 310 &&
       fUtils->getGrandMotherLabel(v0->dLabel[0]) == -1){      
    
      fHistSelV0PV0DecayLengthXY_TrueK0s   -> Fill(p, v0->DecayLengthXY);
      fHistSelV0PV0DecayLength_TrueK0s     -> Fill(p, v0->DecayLengthXYZ);
      fHistSelV0PV0PointingAngleXY_TrueK0s -> Fill(p, v0->CpvXY);
      fHistSelV0PV0PointingAngle_TrueK0s   -> Fill(p, v0->CpvXYZ);
      fHistSelV0PV0DCAXY_TrueK0s           -> Fill(p, v0->DcaXY);
      fHistSelV0PV0DCA_TrueK0s             -> Fill(p, v0->DcaXYZ);
      fHistSelV0PV0TrackDistanceXY_TrueK0s -> Fill(p, v0->DcaDaughtersXY);
      fHistSelV0PV0TrackDistance_TrueK0s   -> Fill(p, v0->DcaDaughtersXYZ);
      fHistSelV0PV0Chi2perNDF_TrueK0s      -> Fill(p, v0->rChi2V0);
      fHistSelV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
      fHistSelV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
      fHistSelV0TrackDCA_TrueK0s           -> Fill(p, DCApair);      

      if (fabs(v0->Charge) > 0) {
	fHistLSSelV0PV0DecayLengthXY_TrueK0s   -> Fill(p, v0->DecayLengthXY);
	fHistLSSelV0PV0DecayLength_TrueK0s     -> Fill(p, v0->DecayLengthXYZ);
	fHistLSSelV0PV0PointingAngleXY_TrueK0s -> Fill(p, v0->CpvXY);
	fHistLSSelV0PV0PointingAngle_TrueK0s   -> Fill(p, v0->CpvXYZ);
	fHistLSSelV0PV0DCAXY_TrueK0s           -> Fill(p, v0->DcaXY);
	fHistLSSelV0PV0DCA_TrueK0s             -> Fill(p, v0->DcaXYZ);
	fHistLSSelV0PV0TrackDistanceXY_TrueK0s -> Fill(p, v0->DcaDaughtersXY);
	fHistLSSelV0PV0TrackDistance_TrueK0s    -> Fill(p, v0->DcaDaughtersXYZ);
	fHistLSSelV0PV0Chi2perNDF_TrueK0s      -> Fill(p, v0->rChi2V0);
	fHistLSSelV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
	fHistLSSelV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
	fHistLSSelV0TrackDCA_TrueK0s           -> Fill(p, DCApair);
      }

    }    
  } else {
    
    fHistV0PV0DecayLengthXY   -> Fill(p, v0->DecayLengthXY);
    fHistV0PV0DecayLength     -> Fill(p, v0->DecayLengthXYZ);
    fHistV0PV0PointingAngleXY -> Fill(p, v0->CpvXY);
    fHistV0PV0PointingAngle   -> Fill(p, v0->CpvXYZ);
    fHistV0PV0DCAXY           -> Fill(p, v0->DcaXY);
    fHistV0PV0DCA             -> Fill(p, v0->DcaXYZ);
    fHistV0PV0TrackDistanceXY -> Fill(p, v0->DcaDaughtersXY);
    fHistV0PV0TrackDistance   -> Fill(p, v0->DcaDaughtersXYZ);
    fHistV0PV0Chi2perNDF      -> Fill(p, v0->rChi2V0);
    fHistV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
    fHistV0TrackDCAXY         -> Fill(p, DCApairXY);
    fHistV0TrackDCA           -> Fill(p, DCApair);

    if (fabs(v0->Charge) > 0) {
      fHistLSV0PV0DecayLengthXY   -> Fill(p, v0->DecayLengthXY);
      fHistLSV0PV0DecayLength     -> Fill(p, v0->DecayLengthXYZ);
      fHistLSV0PV0PointingAngleXY -> Fill(p, v0->CpvXY);
      fHistLSV0PV0PointingAngle   -> Fill(p, v0->CpvXYZ);
      fHistLSV0PV0DCAXY           -> Fill(p, v0->DcaXY);
      fHistLSV0PV0DCA             -> Fill(p, v0->DcaXYZ);
      fHistLSV0PV0TrackDistanceXY -> Fill(p, v0->DcaDaughtersXY);
      fHistLSV0PV0TrackDistance    -> Fill(p, v0->DcaDaughtersXYZ);
      fHistLSV0PV0Chi2perNDF      -> Fill(p, v0->rChi2V0);
      fHistLSV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
      fHistLSV0TrackDCAXY         -> Fill(p, DCApairXY);
      fHistLSV0TrackDCA           -> Fill(p, DCApair);
    }

    if(fIsMC && 
       fUtils->isSameMotherPair(v0->dLabel[0],v0->dLabel[1]) &&
       fUtils->getMotherPdgCode(v0->dLabel[0]) == 310 &&
       fUtils->getGrandMotherLabel(v0->dLabel[0]) == -1){      
    
      fHistV0PV0DecayLengthXY_TrueK0s   -> Fill(p, v0->DecayLengthXY);
      fHistV0PV0DecayLength_TrueK0s     -> Fill(p, v0->DecayLengthXYZ);
      fHistV0PV0PointingAngleXY_TrueK0s -> Fill(p, v0->CpvXY);
      fHistV0PV0PointingAngle_TrueK0s   -> Fill(p, v0->CpvXYZ);
      fHistV0PV0DCAXY_TrueK0s           -> Fill(p, v0->DcaXY);
      fHistV0PV0DCA_TrueK0s             -> Fill(p, v0->DcaXYZ);
      fHistV0PV0TrackDistanceXY_TrueK0s -> Fill(p, v0->DcaDaughtersXY);
      fHistV0PV0TrackDistance_TrueK0s   -> Fill(p, v0->DcaDaughtersXYZ);
      fHistV0PV0Chi2perNDF_TrueK0s      -> Fill(p, v0->rChi2V0);
      fHistV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
      fHistV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
      fHistV0TrackDCA_TrueK0s           -> Fill(p, DCApair);      

      if (fabs(v0->Charge) > 0) {
	fHistLSV0PV0DecayLengthXY_TrueK0s   -> Fill(p, v0->DecayLengthXY);
	fHistLSV0PV0DecayLength_TrueK0s     -> Fill(p, v0->DecayLengthXYZ);
	fHistLSV0PV0PointingAngleXY_TrueK0s -> Fill(p, v0->CpvXY);
	fHistLSV0PV0PointingAngle_TrueK0s   -> Fill(p, v0->CpvXYZ);
	fHistLSV0PV0DCAXY_TrueK0s           -> Fill(p, v0->DcaXY);
	fHistLSV0PV0DCA_TrueK0s             -> Fill(p, v0->DcaXYZ);
	fHistLSV0PV0TrackDistanceXY_TrueK0s -> Fill(p, v0->DcaDaughtersXY);
	fHistLSV0PV0TrackDistance_TrueK0s    -> Fill(p, v0->DcaDaughtersXYZ);
	fHistLSV0PV0Chi2perNDF_TrueK0s      -> Fill(p, v0->rChi2V0);
	fHistLSV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p);    
	fHistLSV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
	fHistLSV0TrackDCA_TrueK0s           -> Fill(p, DCApair);
      }

    }    
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::K0sPairAnalysis(AliPID::EParticleType pid1,AliPID::EParticleType pid2) {  
  

  Int_t nV0 = fArrayK0s->GetEntriesFast();
  
  int iLeading = -1;  

  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {

    K0sContainer* v0_1 = static_cast<K0sContainer*>(fArrayK0s->At(iV0_1));

    unique_ptr<TVector3> vec1(new TVector3());
    vec1->SetXYZ(v0_1->P[0],v0_1->P[1],v0_1->P[2]);

    unique_ptr<TLorentzVector> lv1(new TLorentzVector());    
    lv1->SetPtEtaPhiM(vec1->Pt(), vec1->PseudoRapidity(), vec1->Phi(), TDatabasePDG::Instance()->GetParticle(310)->Mass());
    
    for (int iV0_2 = iV0_1 + 1; iV0_2 < nV0; ++iV0_2) {
      
      K0sContainer *v0_2 = static_cast<K0sContainer*>(fArrayK0s->At(iV0_2));
      
      if (!isAcceptK0sPair(v0_1,v0_2)) continue; 
      
      unique_ptr<TVector3> vec2(new TVector3());
      vec2->SetXYZ(v0_2->P[0],v0_2->P[1],v0_2->P[2]);
      
      unique_ptr<TLorentzVector> lv2(new TLorentzVector());
      lv2->SetPtEtaPhiM(vec2->Pt(), vec2->PseudoRapidity(), vec2->Phi(), TDatabasePDG::Instance()->GetParticle(310)->Mass());
      
      unique_ptr<TLorentzVector> lv12(new TLorentzVector(lv1->Vect()+lv2->Vect(),lv1->T()+lv2->T()));
      
      double fill[] = {lv12->M(), lv12->Pt(), fCent};      
      
      if (fUtils->isAcceptK0sCandidateMassRange(v0_1->Mass) && fUtils->isAcceptK0sCandidateMassRange(v0_2->Mass)) {
	
	if ( v0_1->Charge == 0 && v0_2->Charge == 0 ) {
	  fSparseNeutralK0sPair->Fill(fill);
	} else if ( (v0_1->Charge == 0 && v0_2->Charge < 0) || (v0_1->Charge < 0 && v0_2->Charge == 0) ) {
	  fSparseNeutralNegativeK0sPair->Fill(fill);
	} else if ( (v0_1->Charge == 0 && v0_2->Charge > 0) || (v0_1->Charge > 0 && v0_2->Charge == 0) ) {
	  fSparseNeutralPositiveK0sPair->Fill(fill);
	} else if ( v0_1->Charge > 0 && v0_2->Charge > 0 ) {
	  fSparsePositiveK0sPair->Fill(fill);
	} else if ( v0_1->Charge < 0 && v0_2->Charge < 0 ) {
	  fSparseNegativeK0sPair->Fill(fill);
	} else if ( (v0_1->Charge < 0 && v0_2->Charge > 0) || (v0_1->Charge > 0 && v0_2->Charge < 0) ) {
	  fSparsePositiveNegativeK0sPair->Fill(fill);
	}
      
      }
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::K0sPairAnalysisEventMixing(AliPID::EParticleType pid1, AliPID::EParticleType pid2) {
  
  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner(true);
  
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
  
  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgrK0s->GetEventPool(poolCent, poolVtxZ, poolPsi);
  
  Int_t nV0 = fEvent->GetNumberOfV0s();

  AliAODv0 *v0_1;
  AliAODv0 *v0_2;
  
  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {
    
    v0_1 = static_cast<AliAODv0 *>(fEvent->GetV0(iV0_1));
    
    if (pool->IsReady()) {
      
      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {
	
        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);
	
        for (int iV0_2 = 0; iV0_2 < poolTracks->GetEntriesFast(); ++iV0_2) {
	  
        }
      }
    }
    
    fTrackArray->Add(v0_1);
  }

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner(true);
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  //delete fTrackArray;
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::AddK0sArray(AliPID::EParticleType pid1,AliPID::EParticleType pid2) {  
  
  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner(true);

  double poolCent = 0.;
  double poolVtxZ = 0.;
  double poolPsi  = 0.;
  
  if (onEvtMixingPoolVtxZ) poolVtxZ = fPrimVtxPos[2];
  if (onEvtMixingPoolCent) poolCent = fCent;
  if (onEvtMixingPoolPsi)  poolPsi  = fPsi;
 
  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgrPion->GetEventPool(poolCent, poolVtxZ, poolPsi);
   
  int nTrack = fEvent->GetNumberOfTracks();
  
  if ( fIsManualV0Analysis ) {
    
    for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {
      
      AliAODTrack* aodTrack1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);
      if(!aodTrack1) continue;

      TrackPIDChecker(aodTrack1, pid1, false);
      
      if ( !fUtils->isAcceptTrackKinematics(aodTrack1) ) continue;
      if ( !fUtils->isAcceptV0TrackQuality(aodTrack1) )  continue;
      if ( !fUtils->isAcceptMidPid(aodTrack1,pid1) )     continue;
      
      TrackQualityChecker(aodTrack1);      
      TrackPIDChecker(aodTrack1, pid1, true);
      
      for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTrack; ++iTrack2) {
	
	AliAODTrack* aodTrack2 = (AliAODTrack *)fEvent->GetTrack(iTrack2);
	if(!aodTrack2) continue;

	if ( !fUtils->isAcceptTrackKinematics(aodTrack2) ) continue;
	if ( !fUtils->isAcceptV0TrackQuality(aodTrack2) )  continue;
	if ( !fUtils->isAcceptMidPid(aodTrack2,pid2) )     continue;
	
	
	bool isSignalFlag = false;;	
	if (fUtils->isSameMotherPair(aodTrack1,aodTrack2)){
	  if(fUtils->getMotherPdgCode(aodTrack1) == 310){	 
	    isSignalFlag = true;
	  }
	}
	
	aodTrack1->SetID(iTrack1);
	aodTrack2->SetID(iTrack2);
	
	K0sContainer* v0 = calcK0sFromTracks(aodTrack1,aodTrack2);	
	if (!v0) continue;
	
	V0QualityChecker(v0, false);	
	
	unique_ptr<TVector3> vec(new TVector3());
	vec->SetXYZ(v0->P[0],v0->P[1],v0->P[2]);
	
	unique_ptr<TLorentzVector> lv(new TLorentzVector());
	lv->SetPtEtaPhiM(vec->Pt(), vec->PseudoRapidity(), vec->Phi(), v0->Mass);
	
	double fill[] = {lv->M(), lv->Pt(), fCent, 0};
	
	if ( !isAcceptK0sKinematics(v0) )              continue;	
	
	if ( v0->Charge == 0 )     fSparseULSPionPairBeforeCuts  -> Fill(fill);
	else if ( v0->Charge > 0 ) fSparseLSppPionPairBeforeCuts -> Fill(fill);
	else if ( v0->Charge < 0 ) fSparseLSmmPionPairBeforeCuts -> Fill(fill);
	
	bool *isAcceptCuts = new bool[6];
	
	isAcceptV0QualityCuts(v0,isAcceptCuts);
	
	/*
	if ( v0->Charge == 0 ) {
	  if (isAcceptCuts[0]) fSparseULSPionPair_PassChi2perNDFCut         -> Fill(fill);
	  if (isAcceptCuts[1]) fSparseULSPionPair_PassCPVXYCut              -> Fill(fill);
	  if (isAcceptCuts[2]) fSparseULSPionPair_PassDecayLengthXYCut      -> Fill(fill);
	  if (isAcceptCuts[3]) fSparseULSPionPair_PassDaughterDistanceXYCut -> Fill(fill);
	  if (isAcceptCuts[4]) fSparseULSPionPair_PassDCAXYCut              -> Fill(fill);
	  if (isAcceptCuts[5]) fSparseULSPionPair_PassLifetimeCut           -> Fill(fill);
	}
	*/
	
	//if ( !isAcceptV0QualityCuts(v0,isAcceptCuts) ) continue;
	
	V0QualityChecker(v0, true);
	
	delete [] isAcceptCuts;
	
	if ( v0->Charge == 0 )     fSparseULSPionPair  -> Fill(fill);
	else if ( v0->Charge > 0 ) fSparseLSppPionPair -> Fill(fill);
	else if ( v0->Charge < 0 ) fSparseLSmmPionPair -> Fill(fill);

	fArrayK0s->Add(v0);
      }
      
      if (pool->IsReady()) {
	continue;
	for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {
	
	  TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);
	  
	  for (int iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast(); ++iTrack2) {
	    
	    AliAODTrack* aodTrack2 = dynamic_cast<AliAODTrack*>(poolTracks->At(iTrack2));
	    
	    K0sContainer* v0 = calcK0sFromTracks(aodTrack1,aodTrack2);
	    
	    unique_ptr<TVector3> vec(new TVector3());
	    vec->SetXYZ(v0->P[0],v0->P[1],v0->P[2]);
	    
	    double pt = v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1] > 0 ? sqrt(v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1]) : 0;
	    
	    unique_ptr<TLorentzVector> lv(new TLorentzVector());
	    lv->SetPtEtaPhiM(pt, vec->PseudoRapidity(), vec->Phi(), v0->Mass);
	    
	    double fill[] = {lv->M(), lv->Pt(), fUtils->getCentClass(), 0};
	
	    if ( v0->Charge == 0 )     fSparseMixULSPionPair->Fill(fill);
	    else if ( v0->Charge > 0 ) fSparseMixLSppPionPair->Fill(fill);
	    else if ( v0->Charge < 0 ) fSparseMixLSmmPionPair->Fill(fill);
	    
	  }//end of loop iTrack2
	}//end of loop iMixEvt	  	
      }
      fTrackArray->Add(aodTrack1);
    }//end of loop iTrack1
  
  } else {
    
    Int_t nV0 = fEvent->GetNumberOfV0s();
    
    for (int iV0 = 0; iV0 < nV0; ++iV0) {
      /*
      AliAODv0* _v0_ = (AliAODv0 *)fEvent->GetV0(iV0);
      
      AliAODTrack *aodTrack1 = (AliAODTrack *)_v0_ -> GetDaughter(0);
      AliAODTrack *aodTrack2 = (AliAODTrack *)_v0_ -> GetDaughter(1);      
      
      TrackPIDChecker(aodTrack1, pid1,false);
      
      if ( !fUtils->isAcceptTrackKinematics(aodTrack1) ) continue;
      if ( !fUtils->isAcceptV0TrackQuality(aodTrack1) )  continue;
      if ( !fUtils->isAcceptMidPid(aodTrack1,pid1) )     continue;
      
      TrackQualityChecker(aodTrack1);
      TrackPIDChecker(aodTrack1, pid1, true);

      if ( !fUtils->isAcceptTrackKinematics(aodTrack2) ) continue;
      if ( !fUtils->isAcceptV0TrackQuality(aodTrack2) )  continue;
      if ( !fUtils->isAcceptMidPid(aodTrack2,pid2) )     continue;

      AliAODVertex *primeVtx = new AliAODVertex();
      primeVtx->SetPosition(fUtils->getVtxX(),fUtils->getVtxY(),fUtils->getVtxZ());
      
      double Mom1[] = {aodTrack1->Px(),aodTrack1->Py(),aodTrack1->Pz()};
      double Mom2[] = {aodTrack2->Px(),aodTrack2->Py(),aodTrack2->Pz()};
      
      double rDcaDaughterToPrimVertex[] = {0,0};
      double rDcaV0ToPrimVertex = _v0_->DcaV0ToPrimVertex();
      
      AliAODVertex *secondVtx = new AliAODVertex();
      secondVtx->SetPosition(_v0_->Xv(),_v0_->Yv(),_v0_->Zv());
      secondVtx->AddDaughter(aodTrack1);
      secondVtx->AddDaughter(aodTrack2);
      
      AliAODv0* v0=(AliAODv0*)_v0_->Clone();      
      v0->SetCharge(aodTrack1->Charge() + aodTrack2->Charge());
      v0->SetSecondaryVtx(secondVtx);
      */
      /*
      if(!updateAODv0(v0)) {
	continue;
      }
      */
      
    }//end of loop iV0
   
  }
  ///*
  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner(true);
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }
  //*/

  return true;

}

K0sContainer* AliAnalysisTaskAODTrackPair::calcK0sFromTracks(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2){
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  bool isSignalFlag = false;
  if(fUtils->getMotherPdgCode(aodTrack1) == 310){
    if (fUtils->isSameMotherPair(aodTrack1,aodTrack2)){
      all++;
      cout<<"============================================================================================================"<<endl;
      isSignalFlag = true;
    }
  }
  
  double MagF = fEvent->GetMagneticField();
  AliKFParticle::SetField(MagF);

  KFPVertex kfVtxPrimeVertex;
  double fPrimVtxCov[6]={};
  fPrimVtx->GetCovarianceMatrix(fPrimVtxCov);
  
  float fPrimVtxPosF[3],fPrimVtxCovF[6];    
  for(int iEl = 0; iEl < 3; iEl++)
    fPrimVtxPosF[iEl] = (float)fPrimVtxPos[iEl];
  for(int iEl = 0; iEl < 6; iEl++)
    fPrimVtxCovF[iEl] = (float)fPrimVtxCov[iEl];  
  kfVtxPrimeVertex.SetXYZ((float)fPrimVtxPos[0], (float)fPrimVtxPos[1], (float)fPrimVtxPos[2]);
  kfVtxPrimeVertex.SetCovarianceMatrix(fPrimVtxCovF);    
  kfVtxPrimeVertex.SetChi2(fPrimVtx->GetChi2());
  kfVtxPrimeVertex.SetNDF(fPrimVtx->GetNDF());
  kfVtxPrimeVertex.SetNContributors(fPrimVtx->GetNContributors());
  
  KFParticle kfParticlePrimeVertex(kfVtxPrimeVertex);

  int pdg1=211;
  int pdg2=211;
  
  KFParticle kfTrack[2];
  kfTrack[0] = CreateKFParticle(aodTrack1, pdg1);
  kfTrack[1] = CreateKFParticle(aodTrack2, pdg1);

  bool checkval=true;
  bool checkval1=true;
  bool checkval2=true;
  int fRecoMethod = 0;

  KFParticle kfTrack1 = CreateKFParticle(aodTrack1, pdg1);
  KFParticle kfTrack2 = CreateKFParticle(aodTrack2, pdg1);
    
  KFParticleK0s kfCheckIfSafePair;
  kfCheckIfSafePair.AddDaughter(kfTrack[1]);
  if (!kfCheckIfSafePair.CheckIfSafePair(kfTrack[0])) {
    if (isSignalFlag) ++abord;
    return NULL;
  }
    
  KFParticle kfK0s;
  kfK0s.Initialize();
  kfK0s.SetConstructMethod(2); 
  //kfK0s.SetConstructMethod(0); 
  
  const KFParticle *kfKs0Daughters[2] = {&kfTrack[0],&kfTrack[1]};  
  cout<<"[1.0.0.0]"<<endl;
  kfK0s.Construct(kfKs0Daughters, 2);
  cout<<"[2.0.0.0]"<<endl;
  float chi2K0s   = kfK0s.GetChi2() / kfK0s.GetNDF();
  cout<<"[3.0.0.0]"<<endl;
  float momK0s[3] = {kfK0s.GetPx(),kfK0s.GetPy(),kfK0s.GetPz()};
  cout<<"[4.0.0.0]"<<endl;
  float posK0s[3] = {};
  float covK0s[36]= {};
  posK0s[0] = kfK0s.GetX();
  posK0s[1] = kfK0s.GetY();
  posK0s[2] = kfK0s.GetZ();
  cout<<"[5.0.0.0]"<<endl;
  float QtAlfa[2]={};
  kfK0s.GetArmenterosPodolanski(kfTrack[0],kfTrack[1],QtAlfa);
  cout<<"[6.0.0.0]"<<endl;
  for (int index=0; index<36; ++index) covK0s[index] = kfK0s.GetCovariance(index);
  cout<<"[7.0.0.0]"<<endl;
  double decay_length    = DecayLengthFromKF(kfK0s,kfParticlePrimeVertex);
  double decay_length_xy = DecayLengthXYFromKF(kfK0s,kfParticlePrimeVertex);
  cout<<"[8.0.0.0]"<<endl;
  float dcaPointK0s[8]={}, dcaPointK0sCov[36]={};  
  //kfK0s.GetParametersAtPoint(fPrimVtxPosF, fPrimVtxCovF, dcaPointK0s, dcaPointK0sCov);
  
  double dcaK0s    = sqrt(pow(fPrimVtxPos[0]-dcaPointK0s[0],2) + pow(fPrimVtxPos[1]-dcaPointK0s[1],2) + pow(fPrimVtxPos[2]-dcaPointK0s[2],2));
  double dcaK0s_xy = sqrt(pow(fPrimVtxPos[0]-dcaPointK0s[0],2) + pow(fPrimVtxPos[1]-dcaPointK0s[1],2));
  cout<<"[9.0.0.0]"<<endl;
  float dcaPointDaughter1[8]={}, dcaPointDaughterCov1[36]={};
  float dcaPointDaughter2[8]={}, dcaPointDaughterCov2[36]={};
  //kfTrack[0].GetParametersAtPoint(fPrimVtxPosF, fPrimVtxCovF, dcaPointDaughter1, dcaPointDaughterCov1);
  //kfTrack[1].GetParametersAtPoint(fPrimVtxPosF, fPrimVtxCovF, dcaPointDaughter2, dcaPointDaughterCov2);
  cout<<"[10.0.0.0]"<<endl;
  double dcaDaughter1;
  double dcaDaughter2;
  double dcaDaughter_xy1;
  double dcaDaughter_xy2;
  
  dcaDaughter1 = sqrt(pow(fPrimVtxPos[0]-dcaPointDaughter1[0],2)+
		      pow(fPrimVtxPos[1]-dcaPointDaughter1[1],2)+
		      pow(fPrimVtxPos[2]-dcaPointDaughter1[2],2));
  dcaDaughter2 = sqrt(pow(fPrimVtxPos[0]-dcaPointDaughter2[0],2)+
		      pow(fPrimVtxPos[1]-dcaPointDaughter2[1],2)+
		      pow(fPrimVtxPos[2]-dcaPointDaughter2[2],2));
  cout<<"[11.0.0.0]"<<endl;
  dcaDaughter_xy1 = sqrt(pow(fPrimVtxPos[0]-dcaPointDaughter1[0],2)+
			 pow(fPrimVtxPos[1]-dcaPointDaughter1[1],2));
  dcaDaughter_xy2 = sqrt(pow(fPrimVtxPos[0]-dcaPointDaughter2[0],2)+
			 pow(fPrimVtxPos[1]-dcaPointDaughter2[1],2));
  cout<<"[12.0.0.0]"<<endl;
  double daughterDistance    = sqrt(pow(dcaPointDaughter1[0]-dcaPointDaughter2[0],2)+
				    pow(dcaPointDaughter1[1]-dcaPointDaughter2[1],2)+
				    pow(dcaPointDaughter1[2]-dcaPointDaughter2[2],2));
  double daughterDistance_xy = sqrt(pow(dcaPointDaughter1[0]-dcaPointDaughter2[0],2)+
				    pow(dcaPointDaughter1[1]-dcaPointDaughter2[1],2));
  cout<<"[13.0.0.0]"<<endl;
  //kfTrack[0].TransportToPoint(posK0s);
  //kfTrack[1].TransportToPoint(posK0s);

  /*
  TLorentzVector lv1;
  TLorentzVector lv2;  
  lv1.SetPtEtaPhiM(sqrt(kfTrack[0].GetPx()*kfTrack[0].GetPx() + kfTrack[0].GetPy()*kfTrack[0].GetPy()),kfTrack[0].GetEta(),kfTrack[0].GetPhi(),TDatabasePDG::Instance()->GetParticle(211)->Mass());
  lv2.SetPtEtaPhiM(sqrt(kfTrack[1].GetPx()*kfTrack[1].GetPx() + kfTrack[1].GetPy()*kfTrack[1].GetPy()),kfTrack[1].GetEta(),kfTrack[1].GetPhi(),TDatabasePDG::Instance()->GetParticle(211)->Mass());
  TLorentzVector lv12 = lv1 + lv2;
  */

  unique_ptr<TVector3> momS(new TVector3(momK0s[0],momK0s[1],momK0s[2]));
  unique_ptr<TVector3> vecS(new TVector3(posK0s[0]-fPrimVtxPos[0],posK0s[1]-fPrimVtxPos[1],posK0s[2]-fPrimVtxPos[2]));
  double kfSpv  = momS->Angle(*vecS.get());
  double kfScpv = TMath::Cos(kfSpv);  
  momS->SetZ(0.);
  vecS->SetZ(0.); 
  double kfSpvXY                 = momS->Angle(*vecS.get());
  double kfScpvXY                = TMath::Cos(kfSpvXY);
  cout<<"[14.0.0.0]"<<endl;
  double kfMassK0s               = kfK0s.GetMass();
  cout<<"[15.0.0.0]"<<endl;
  double kfEtaK0s                = 0;//kfK0s.GetEta();  
  cout<<"[16.0.0.0]"<<endl;
  double kfMomK0s[]              = {kfK0s.GetPx(),kfK0s.GetPy(),kfK0s.GetPz()};
  double kfVtxK0s[]              = {kfK0s.GetX(),kfK0s.GetY(),kfK0s.GetZ()};
  cout<<"[17.0.0.0]"<<endl;
  int kfQK0s                     = kfK0s.Q();
  double kfDecayLengthXYK0s      = decay_length_xy;
  double kfDecayLengthK0s        = decay_length;
  double kfImpParXYK0s           = dcaK0s_xy;
  double kfImpParK0s             = dcaK0s;
  double kfCosPointVectorXYK0s   = kfScpvXY;
  double kfCosPointVectorK0s     = kfScpv;
  double kfDaughterDistanceXYK0s = kfTrack[0].GetDistanceFromParticleXY(kfTrack[1]);
  double kfDaughterDistanceK0s   = 0;//daughterDistance;
  double kfChi2perNdfK0s         = chi2K0s;
  int kfTrackIdDaughters[]       = {aodTrack1->GetID(),aodTrack2->GetID()};
  int kfTrackLabelDaughters[]    = {aodTrack1->GetLabel(),aodTrack2->GetLabel()};
  double kfMomDaughter1[]        = {kfTrack[0].Px(),kfTrack[0].Py(),kfTrack[0].Pz()};
  double kfMomDaughter2[]        = {kfTrack[1].Px(),kfTrack[1].Py(),kfTrack[1].Pz()};
  double kfImpParXYDaughters[]   = {dcaDaughter_xy1,dcaDaughter_xy2};
  double kfImpParDaughters[]     = {dcaDaughter1,dcaDaughter2};
  cout<<"[18.0.0.0]"<<endl;
  K0sContainer* container = new K0sContainer(kfMassK0s,kfEtaK0s,kfMomK0s,kfVtxK0s,kfQK0s,kfDecayLengthXYK0s,kfDecayLengthK0s,
					     kfImpParXYK0s,kfImpParK0s,kfCosPointVectorXYK0s,kfCosPointVectorK0s,
					     kfDaughterDistanceXYK0s,kfDaughterDistanceK0s,
					     kfChi2perNdfK0s,kfTrackIdDaughters,kfTrackLabelDaughters,
					     kfMomDaughter1,kfMomDaughter2,kfImpParXYDaughters,kfImpParDaughters);
  cout<<"[19.0.0.0]"<<endl;
  if (0) {

    AliAODMCParticle* particle1 = (AliAODMCParticle *)fMCTrackArray->At(aodTrack1->GetLabel());
    
    if(particle1) {
      AliAODMCParticle* mother = (AliAODMCParticle *)fMCTrackArray->At(particle1->GetMother());
      if (mother) {
	if (mother->GetMother()<0) {
	  if (isSignalFlag) cout<<"=================================================================================================================="<<endl;
	  cout<<"GM Label  :"<<fUtils->getGrandMotherLabel(kfTrackLabelDaughters[0])<<endl;
	  cout<<"MC VTX    : "<<particle1->Xv()<<"   "<<particle1->Yv()<<"   "<<particle1->Zv()<<endl;
	  cout<<"KF VTX    : "<<kfK0s.GetX()<<"   "<<kfK0s.GetY()<<"   "<<kfK0s.GetZ()<<endl;
	  cout<<"MC MOM    : "<<mother->Px()<<"   "<<mother->Py()<<"   "<<mother->Pz()<<endl;
	  cout<<"KF MOM    : "<<kfK0s.Px()<<"   "<<kfK0s.Py()<<"   "<<kfK0s.Pz()<<endl;
	  cout<<"KF MASS   : "<<kfMassK0s<<endl;
	  cout<<"KF dDis   : "<<kfDaughterDistanceXYK0s<<endl;
	  cout<<endl;	  
	}
      }
    }
  }
  
  return container;
}

bool AliAnalysisTaskAODTrackPair::updateAODv0(AliAODv0* v0){
  /*
  AliAODTrack *aodTrack1 = (AliAODTrack *)v0->GetDaughter(0);
  AliAODTrack *aodTrack2 = (AliAODTrack *)v0->GetDaughter(1);
  
  AliAODv0 *new_v0 = calcK0sFromTracks(aodTrack1,aodTrack2);
  
  double dca_daughters = new_v0->DcaV0Daughters();
  
  double Mom1[] = {new_v0->MomPosX(),new_v0->MomPosY(),new_v0->MomPosZ()};
  double Mom2[] = {new_v0->MomNegX(),new_v0->MomNegY(),new_v0->MomNegZ()};

  double rDcaDaughterToPrimVertex[] = {new_v0->DcaPosToPrimVertex(),new_v0->DcaNegToPrimVertex()};
  double rDcaV0ToPrimVertex = new_v0->DcaV0ToPrimVertex();
  
  v0->Fill(new_v0->GetSecondaryVtx(),dca_daughters,rDcaV0ToPrimVertex,Mom1,Mom2,rDcaDaughterToPrimVertex);

  delete new_v0;
  */
  return true;
}

int AliAnalysisTaskAODTrackPair::findLeadingTrack(){
  
  return 1;
  int nTrack = fEvent->GetNumberOfTracks();
  
  for (int iTrack = 0; iTrack < nTrack; ++iTrack){
    
  }

}

KFParticle AliAnalysisTaskAODTrackPair::CreateKFParticle(AliAODTrack* track, int pdg){
  

  AliKFParticle alikfp(*track,pdg);

  float fkP[8]  = {};
  float fkC[21] = {};

  for (int index=0; index<6; index++) {
    fkP[index] = alikfp.GetParameter(index);
  }
  for (int index=0; index<21; index++) {
    fkC[index] = alikfp.GetCovariance(index);
  }
  
  double fP[6]={};
  double fC[21]={};

  track->GetXYZ(fP);
  track->GetPxPyPz(&fP[3]);
  track->GetCovarianceXYZPxPyPz(fC);
  
  KFPTrack kfpt;
  kfpt.SetParameters(fP);
  kfpt.SetCovarianceMatrix(fC);
  kfpt.SetCharge(track->Charge());
  kfpt.SetNDF(1);
  kfpt.SetChi2(track->Chi2perNDF());
  
  KFParticle kfp;
  kfp.Initialize();  
  //kfp = KFParticle(kfpt);
  kfp.Create(fkP,fkC,track->Charge(),  TDatabasePDG::Instance()->GetParticle(pdg)->Mass());

  return kfp;
}

bool AliAnalysisTaskAODTrackPair::CheckTrackCovariance(AliAODTrack* track){
  
  Double_t covMatrix[21];
  track->GetCovarianceXYZPxPyPz(covMatrix);
  Double_t cov[6][6]={0.};
  cov[0][0] = covMatrix[0];
  cov[1][0] = covMatrix[1];
  cov[1][1] = covMatrix[2];
  cov[2][0] = covMatrix[3];
  cov[2][1] = covMatrix[4];
  cov[2][2] = covMatrix[5];
  cov[3][0] = covMatrix[6];
  cov[3][1] = covMatrix[7];
  cov[3][2] = covMatrix[8];
  cov[3][3] = covMatrix[9];
  cov[4][0] = covMatrix[10];
  cov[4][1] = covMatrix[11];
  cov[4][2] = covMatrix[12];
  cov[4][3] = covMatrix[13];
  cov[4][4] = covMatrix[14];
  cov[5][0] = covMatrix[15];
  cov[5][1] = covMatrix[16];
  cov[5][2] = covMatrix[17];
  cov[5][3] = covMatrix[18];
  cov[5][4] = covMatrix[19];
  cov[5][5] = covMatrix[20];

  if ( cov[0][0]<0 || cov[1][1]<0 || cov[2][2]<0 || cov[3][3]<0 || cov[4][4]<0 || cov[5][5]<0 ) return kFALSE;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(cov[i][j]) > TMath::Sqrt(cov[i][i]*cov[j][j]) ) return false;
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::CheckTrackCovariance(KFParticle kfp){
  if ( kfp.GetCovariance(0,0)<0 || kfp.GetCovariance(1,1)<0 || kfp.GetCovariance(2,2)<0 || 
       kfp.GetCovariance(3,3)<0 || kfp.GetCovariance(4,4)<0 || kfp.GetCovariance(5,5)<0 ) return false;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(kfp.GetCovariance(i,j)) > TMath::Sqrt(kfp.GetCovariance(i,i)*kfp.GetCovariance(j,j)) ) return false;
    }
  }
  return true;
}

double AliAnalysisTaskAODTrackPair::DecayLengthFromKF(KFParticle kfpParticle, KFParticle PV){
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t dz_particle = PV.GetZ()-kfpParticle.GetZ();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle + dz_particle*dz_particle);
  return l_particle;
}
double AliAnalysisTaskAODTrackPair::DecayLengthXYFromKF(KFParticle kfpParticle, KFParticle PV){
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle);
  return l_particle;
}

bool AliAnalysisTaskAODTrackPair::ProcessMC(){
  
  double primRecVtx[3]={fUtils->getVtxX(),fUtils->getVtxY(),fUtils->getVtxZ()};

  AliAODMCParticle *particle1;
  AliAODMCParticle *particle2;
  AliAODMCParticle *particle12;
  
  Int_t nV0 = fArrayK0s->GetEntriesFast();
  
  vector<int> plabels1;
  vector<int> plabels2;
  vector<int> plabels12;
  vector<int> v0label1;
  vector<int> v0label2;
  vector<int> v0label12;
  
  for (int iV0 = 0; iV0 < nV0; ++iV0) {
    
    K0sContainer* v0 = (K0sContainer *)fArrayK0s->At(iV0);
    
    plabels12.push_back(v0->dLabel[0]);
    plabels12.push_back(v0->dLabel[1]);
    v0label12.push_back(iV0);
    v0label12.push_back(iV0);
    
    plabels1.push_back(v0->dLabel[0]);
    plabels2.push_back(v0->dLabel[1]);
    v0label1.push_back(iV0);
    v0label2.push_back(iV0);
  }
  
  for (Int_t iTrack1 = 0; iTrack1 < fMCTrackArray->GetEntries(); ++iTrack1) {
    
    particle1 = (AliAODMCParticle *)fMCTrackArray->At(iTrack1);
    
    if (!particle1)                       continue;        
    if (particle1->GetMother()>-1)        continue;
    if (particle1->GetPdgCode() != 310)   continue; 
    if (particle1->GetNDaughters() != 2 ) continue;
    if(!isAcceptK0sKinematics(particle1)) continue;
    
    int label1 = particle1->GetDaughterFirst();
    int label2 = particle1->GetDaughterLast();
    
    AliAODMCParticle* p1 = (AliAODMCParticle *)fMCTrackArray->At(label1);
    AliAODMCParticle* p2 = (AliAODMCParticle *)fMCTrackArray->At(label2);

    if(!p1 || !p2)                                                      continue;
    if (fabs(p1->GetPdgCode()) != 211 || fabs(p2->GetPdgCode()) != 211) continue;
    
    double fill_k0s[] = {particle1->M(), particle1->Pt(), fUtils->getCentClass(), 0};
    fSparseTrueK0s->Fill(fill_k0s);

    bool det = false;
    int idv0 = 0;
    
    for (int iV0 = 0; iV0 < nV0; ++iV0) {
      if ((plabels1[iV0]==label1 && plabels2[iV0]==label2) || (plabels1[iV0]==label2 && plabels2[iV0]==label1)) {
	det  = true;
	idv0 = iV0;
      }
    }
    
    if (det){      
      K0sContainer* v0 = (K0sContainer *)fArrayK0s->At(idv0);

      double rx = fabs(p1->Xv())>0 ? v0->Vtx[0]/p1->Xv() : 0;
      double ry = fabs(p1->Yv())>0 ? v0->Vtx[1]/p1->Yv() : 0;
      double rz = fabs(p1->Zv())>0 ? v0->Vtx[2]/p1->Zv() : 0;

      double fill_secondary[] = {particle1->P(), rx, ry, rz};
      
      fSparseTrueK0sRecK0s->Fill(fill_secondary);      
    }
    
  }
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sPair(K0sContainer* v0_1, K0sContainer* v0_2){

  int iTrack1[2];  
  iTrack1[0] = v0_1->dTrackid[0];
  iTrack1[1] = v0_1->dTrackid[1];

  int iTrack2[2];  
  iTrack2[0] = v0_2->dTrackid[0];
  iTrack2[1] = v0_2->dTrackid[1];
  
  if (iTrack1[0] == iTrack2[0] || iTrack1[0] == iTrack2[1] || iTrack1[1] == iTrack2[0] || iTrack1[1] == iTrack2[1]) {
    return false;
  } else {
    return true;
  }

}

bool AliAnalysisTaskAODTrackPair::isAcceptV0QualityCuts(K0sContainer* v0, bool* isAcceptCuts){  
  
  isAcceptCuts[0] = false;
  isAcceptCuts[1] = false;
  isAcceptCuts[2] = false;
  isAcceptCuts[3] = false;
  isAcceptCuts[4] = false;
  isAcceptCuts[5] = false;

  double p  = v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1] + v0->P[2]*v0->P[2] > 0 ? 
    sqrt(v0->P[0]*v0->P[0] + v0->P[1]*v0->P[1] + v0->P[2]*v0->P[2]) : 0;
  
  double max_chi2            = fHistChi2Cut->GetBinContent(fHistChi2Cut->GetXaxis()->FindBin(p));
  double min_cpv             = fHistPointingAngleXYCut->GetBinContent(fHistPointingAngleXYCut->GetXaxis()->FindBin(p));
  double max_lengthXY        = fHistDecayLengthXYCut->GetBinContent(fHistDecayLengthXYCut->GetXaxis()->FindBin(p));
  double max_trackDistanceXY = fHistDaughterTrackDistanceXYCut->GetBinContent(fHistDaughterTrackDistanceXYCut->GetXaxis()->FindBin(p));
  double max_DCAxy           = fHistDCAxyCut->GetBinContent(fHistDCAxyCut->GetXaxis()->FindBin(p));
  double max_lifetime        = fHistLifeTimeCut->GetBinContent(fHistLifeTimeCut->GetXaxis()->FindBin(p));

  double chi2            = v0->rChi2V0;
  double cpv             = v0->CpvXY;  
  double lengthXY        = v0->DecayLengthXY;
  double trackDistanceXY = v0->DcaDaughtersXY;
  double DCAxy           = v0->DcaXY;
  double lifetime        = fUtils->fPdgK0sMass * v0->DecayLengthXYZ / p;

  if ( onChi2Cut ) {
    if ( chi2 < max_chi2 ) isAcceptCuts[0] = true;
    else                   isAcceptCuts[0] = false;
  } else {
    isAcceptCuts[0] = true;
  }
  
  if ( onPointingAngleXYCut ) {
    if ( cpv > min_cpv ) isAcceptCuts[1] = true;
    else                 isAcceptCuts[1] = false;
  } else {
    isAcceptCuts[1] = true;
  }    
  
  if ( onDecayLengthXYCut ) {
    if ( lengthXY < max_lengthXY ) isAcceptCuts[2] = true;
    else                           isAcceptCuts[2] = false;
  } else {
    isAcceptCuts[2] = true;
  }
  
  if ( onDaughterTrackDistanceXYCut ) {
    if ( trackDistanceXY < max_trackDistanceXY ) isAcceptCuts[3] = true;
    else                                         isAcceptCuts[3] = false;
  } else {
    isAcceptCuts[3] = true;
  }
  
  if ( onDCAXYCut ) {
    if ( DCAxy < max_DCAxy )  isAcceptCuts[4] = true;
    else                      isAcceptCuts[4] = false;
  } else {
    isAcceptCuts[4] = true;
  }

  if ( onLifetimeCut ) {
    if ( lifetime < max_lifetime ) isAcceptCuts[5] = true;
    else                           isAcceptCuts[5] = false;
  } else {
    isAcceptCuts[5] = true;
  }
  
  if (isAcceptCuts[0]==true && isAcceptCuts[1]==true && 
      isAcceptCuts[2]==true && isAcceptCuts[3]==true && 
      isAcceptCuts[4]==true && isAcceptCuts[5]==true) {
    return true;
  } else {
    return false;
  }

}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sCuts(K0sContainer* v0){
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sKinematics(K0sContainer* v0) {

  unique_ptr<TVector3> vec(new TVector3());
  vec->SetXYZ(v0->P[0],v0->P[1],v0->P[2]);
  unique_ptr<TLorentzVector> lv(new TLorentzVector());
  lv->SetPtEtaPhiM(vec->Pt(), vec->Eta(), vec->Phi(), v0->Mass);
  if (fMinK0sRap > lv->Rapidity() || fMaxK0sRap < lv->Rapidity()) return false;
  if (fMinK0sPt > vec->Pt() || fMaxK0sPt < vec->Pt())             return false;
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sKinematics(AliAODMCParticle* v0) {

  unique_ptr<TLorentzVector> lv(new TLorentzVector());
  lv->SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), v0->M());  
  
  if (fMinK0sRap > lv->Rapidity() || fMaxK0sRap < lv->Rapidity()) return false;
  if (fMinK0sPt > lv->Pt() || fMaxK0sPt < lv->Pt())             return false;
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptPrimaryTrack(AliAODTrack* track){
  return true;
}

/*
void AliAnalysisTaskAODTrackPair::Construct( const KFParticle *vDaughters[], int nDaughters){
  
  fNDF      = -3;
  
  double fC[36] = {};  
  double fP[8]  = {};
  double fChi2  =  0;
  
  fC[35] = 1;
  
  for( Int_t index=0; index<7;  index++) fP[index] = vDaughters[0]->GetParameter(index);
  for( Int_t index=0; index<28; index++) fC[index] = vDaughters[0]->GetCovariance(index);
  
  //for (int index=0; index<nDaughters; ++index) {
  //AddDaughter( *vDaughters[index] );
  //}
  
}

void AliAnalysisTaskAODTrackPair::AddDaughter(){
  
  

}
*/
