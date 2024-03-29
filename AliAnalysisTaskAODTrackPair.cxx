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
#include <time.h>

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair() : AliAnalysisTaskSE(),
  
  fEvent(NULL),
  fPrimVtx(NULL),
  fLeadingTrack(NULL),

  fPoolMuonTrackMgrK0s(NULL),  
  fPoolMuonTrackMgrPion(NULL),  
  
  fMultSelection(NULL),
  fPIDResponse(NULL),

  fUtils(NULL),
  
  fMCTrackArray(NULL),
  fArrayK0s(NULL),
  fArrayK0sDaughterTrack(NULL),
  fArrayPrimaryPionTrack(NULL),
  
  fHistDecayLengthCut(NULL),
  fHistPointingAngleCut(NULL),
  fHistChi2perNDFCut(NULL),
  fHistDaughterPairDCAXYCut(NULL),
  fHistDaughterTrackDistanceXYCut(NULL),
  fHistDCACut(NULL),
  fHistProperLifeTimeCut(NULL),
  
  fTrackDepth(100),
  fPoolSize(100),

  fReadyFraction(0.0),
  
  onVerbose(true),

  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),
  
  onDecayLengthCut(false),
  onPointingAngleCut(true),
  onChi2perNDFCut(true),
  onDaughterPairDCAXYCut(false),
  onDaughterTrackDistanceXYCut(true),
  onDCACut(false),
  onLifetimeCut(false),
  onArmenterosCut(true),
  
  fIsMC(false),
  fIsMixingAnalysis(false),
  fIsManualV0Analysis(true),

  fArmenterosAlpha(0.16),

  fMinVertexCutZ(-10),
  fMaxVertexCutZ(10),

  fMinK0sRap(-1.5),
  fMaxK0sRap(1.5),
  fMinK0sPt(0.0),
  fMaxK0sPt(99.),
  fMinK0sMassRange(0.48833812),
  fMaxK0sMassRange(0.50741788),

  fMinTrackTPCNClusts(70),
  fMinTrackSPDNClusts(0),
  fMaxReducedChi2TPC(4.),
  fMaxReducedChi2ITS(36.),
  fMinCrossRowsFindableRatio(0.8),

  fMinTrackP(0.05),
  fMaxTrackP(999.),
  fMinTrackPt(0.05),
  fMaxTrackPt(999.),
  fMinTrackEta(-0.8),
  fMaxTrackEta(+0.8),  
  fMinLeadingTrackPt(5.),

  fMethodCent("V0A"),
  fPrimVtxPos(),
  fPrimVtxCov(),
  fCent(0.),
  fPsi(0.),
  
  fNK0s(0),

  fMinNContPrimVtx(2),

  ////////////////////////////////////////////////
  // Output histos
  ////////////////////////////////////////////////

  fOutputList(NULL),
  
  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),

  fSparseULSPionPairBeforeCuts(NULL),

  fSparseULSPionPair_PCAcut997_TrkDist1(NULL),
  fSparseULSPionPair_PCAcut999_TrkDist1(NULL),
  fSparseULSPionPair_PCAcut97_TrkDist1(NULL),  
  fSparseULSPionPair_PCAcut95_TrkDist1(NULL),  
  fSparseULSPionPair_PCAcut997_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut999_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut97_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut95_TrkDist04(NULL),
  fSparseULSPionPair_KstarSel(NULL),

  fSparseNeutralK0sPair_PCAcut997_TrkDist1(NULL),
  fSparseNeutralK0sPair_PCAcut999_TrkDist1(NULL),
  fSparseNeutralK0sPair_PCAcut97_TrkDist1(NULL),  
  fSparseNeutralK0sPair_PCAcut95_TrkDist1(NULL),  
  fSparseNeutralK0sPair_PCAcut997_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut999_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut97_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut95_TrkDist04(NULL),
  fSparseNeutralK0sPair_KstarSel(NULL),

  fSparseMixNeutralK0sPair_PCAcut997_TrkDist1(NULL),
  fSparseMixNeutralK0sPair_PCAcut999_TrkDist1(NULL),
  fSparseMixNeutralK0sPair_PCAcut97_TrkDist1(NULL),  
  fSparseMixNeutralK0sPair_PCAcut95_TrkDist1(NULL),  
  fSparseMixNeutralK0sPair_PCAcut997_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut999_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut97_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut95_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_KstarSel(NULL),

  fSparseNeutralK0sPairRejectKstar(NULL),
  fSparseMixNeutralK0sPairRejectKstar(NULL),

  fSparseNeutralK0sPairTransverse(NULL),
  fSparseNeutralK0sPairToward(NULL),
  fSparseNeutralK0sPairAway(NULL),
  fSparseNeutralK0sPairInJet(NULL),

  fSparseK0sPion(NULL),
  fSparseULSPionPairRejectKstar(NULL),

  fSparseRecK0s_PCAcut95_TrkDist1(NULL),
  fSparseRecK0s_PCAcut97_TrkDist1(NULL),
  fSparseRecK0s_PCAcut997_TrkDist1(NULL),
  fSparseRecK0s_PCAcut999_TrkDist1(NULL),
  fSparseRecK0s_PCAcut95_TrkDist04(NULL),
  fSparseRecK0s_PCAcut97_TrkDist04(NULL),
  fSparseRecK0s_PCAcut997_TrkDist04(NULL),
  fSparseRecK0s_PCAcut999_TrkDist04(NULL),  
  fSparseRecK0s_KstarSel(NULL),
  
  fSparseTrueK0s(NULL),
  

  fHistResoK0sMassPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut95_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut97_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut997_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut999_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut95_TrkDist04(NULL),
  
  fHistResoK0sMassPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut97_TrkDist04(NULL),

  fHistResoK0sMassPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut997_TrkDist04(NULL),

  fHistResoK0sMassPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut999_TrkDist04(NULL),

  fHistResoK0sMassPt_KstarSel(NULL),
  fHistResoK0sPtPt_KstarSel(NULL),
  fHistResoK0sEtaPt_KstarSel(NULL),
  fHistResoK0sPhiPt_KstarSel(NULL),

  fHistNumK0sCandidates(NULL),

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
  fHistDCAxy(NULL),
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
  fHistV0Armenteros(NULL),

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
  fHistSelV0Armenteros(NULL),

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
  fHistV0Armenteros_TrueK0s(NULL),
  
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
  fHistSelV0Armenteros_TrueK0s(NULL),
  nEvtAnalyzed(0),
  start(0),
  end(0)
{

}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char *name) :
  AliAnalysisTaskSE(name), 

  fEvent(NULL),
  fPrimVtx(NULL),
  fLeadingTrack(NULL),

  fPoolMuonTrackMgrK0s(NULL),  
  fPoolMuonTrackMgrPion(NULL),  
  
  fMultSelection(NULL),
  fPIDResponse(NULL),

  fUtils(NULL),
  
  fMCTrackArray(NULL),
  fArrayK0s(NULL),
  fArrayK0sDaughterTrack(NULL),
  fArrayPrimaryPionTrack(NULL),

  fHistDecayLengthCut(NULL),
  fHistPointingAngleCut(NULL),
  fHistChi2perNDFCut(NULL),
  fHistDaughterPairDCAXYCut(NULL),
  fHistDaughterTrackDistanceXYCut(NULL),
  fHistDCACut(NULL),
  fHistProperLifeTimeCut(NULL),

  fTrackDepth(100),
  fPoolSize(100),

  fReadyFraction(0.0),
  
  onVerbose(true),

  onEvtMixingPoolVtxZ(true),
  onEvtMixingPoolCent(true),
  onEvtMixingPoolPsi(true),

  onDecayLengthCut(false),
  onPointingAngleCut(true),
  onChi2perNDFCut(true),
  onDaughterPairDCAXYCut(false),
  onDaughterTrackDistanceXYCut(true),
  onDCACut(false),
  onLifetimeCut(false),
  onArmenterosCut(true),

  fIsMC(false),
  fIsMixingAnalysis(false),
  fIsManualV0Analysis(true),

  fArmenterosAlpha(0.16),

  fMinVertexCutZ(-10),
  fMaxVertexCutZ(10),

  fMinK0sRap(-1.5),
  fMaxK0sRap(1.5),
  fMinK0sPt(0.0),
  fMaxK0sPt(99.),
  fMinK0sMassRange(0.48833812),
  fMaxK0sMassRange(0.50741788),

  fMinTrackTPCNClusts(70),
  fMinTrackSPDNClusts(0),
  fMaxReducedChi2TPC(4.),
  fMaxReducedChi2ITS(36.),
  fMinCrossRowsFindableRatio(0.8),
  
  fMinTrackP(0.05),
  fMaxTrackP(999.),
  fMinTrackPt(0.05),
  fMaxTrackPt(999.),
  fMinTrackEta(-0.8),
  fMaxTrackEta(+0.8),
  fMinLeadingTrackPt(5.),

  fMethodCent("V0A"),
  fPrimVtxPos(),
  fPrimVtxCov(),
  fCent(0.),
  fPsi(0.),

  fNK0s(0),

  fMinNContPrimVtx(2),

  ////////////////////////////////////////////////
  // Output histos
  ////////////////////////////////////////////////

  fOutputList(NULL),
  
  fHistEventVtxZ(NULL),
  fHistEventCent(NULL),
  fHistEventMulti(NULL),
  fHistEventVtxCont(NULL),

  fSparseULSPionPairBeforeCuts(NULL),

  fSparseULSPionPair_PCAcut997_TrkDist1(NULL),
  fSparseULSPionPair_PCAcut999_TrkDist1(NULL),
  fSparseULSPionPair_PCAcut97_TrkDist1(NULL),  
  fSparseULSPionPair_PCAcut95_TrkDist1(NULL),  
  fSparseULSPionPair_PCAcut997_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut999_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut97_TrkDist04(NULL),
  fSparseULSPionPair_PCAcut95_TrkDist04(NULL),
  fSparseULSPionPair_KstarSel(NULL),

  fSparseNeutralK0sPair_PCAcut997_TrkDist1(NULL),
  fSparseNeutralK0sPair_PCAcut999_TrkDist1(NULL),
  fSparseNeutralK0sPair_PCAcut97_TrkDist1(NULL),  
  fSparseNeutralK0sPair_PCAcut95_TrkDist1(NULL),  
  fSparseNeutralK0sPair_PCAcut997_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut999_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut97_TrkDist04(NULL),
  fSparseNeutralK0sPair_PCAcut95_TrkDist04(NULL),
  fSparseNeutralK0sPair_KstarSel(NULL),

  fSparseMixNeutralK0sPair_PCAcut997_TrkDist1(NULL),
  fSparseMixNeutralK0sPair_PCAcut999_TrkDist1(NULL),
  fSparseMixNeutralK0sPair_PCAcut97_TrkDist1(NULL),  
  fSparseMixNeutralK0sPair_PCAcut95_TrkDist1(NULL),  
  fSparseMixNeutralK0sPair_PCAcut997_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut999_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut97_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_PCAcut95_TrkDist04(NULL),
  fSparseMixNeutralK0sPair_KstarSel(NULL),

  fSparseNeutralK0sPairRejectKstar(NULL),
  fSparseMixNeutralK0sPairRejectKstar(NULL),

  fSparseNeutralK0sPairTransverse(NULL),
  fSparseNeutralK0sPairToward(NULL),
  fSparseNeutralK0sPairAway(NULL),
  fSparseNeutralK0sPairInJet(NULL),

  fSparseK0sPion(NULL),
  fSparseULSPionPairRejectKstar(NULL),

  fSparseRecK0s_PCAcut95_TrkDist1(NULL),
  fSparseRecK0s_PCAcut97_TrkDist1(NULL),
  fSparseRecK0s_PCAcut997_TrkDist1(NULL),
  fSparseRecK0s_PCAcut999_TrkDist1(NULL),

  fSparseRecK0s_PCAcut95_TrkDist04(NULL),
  fSparseRecK0s_PCAcut97_TrkDist04(NULL),
  fSparseRecK0s_PCAcut997_TrkDist04(NULL),
  fSparseRecK0s_PCAcut999_TrkDist04(NULL),
  fSparseRecK0s_KstarSel(NULL),

  fSparseTrueK0s(NULL),

  fHistResoK0sMassPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut95_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut95_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut97_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut97_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut997_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut997_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sPtPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sEtaPt_PCAcut999_TrkDist1(NULL),
  fHistResoK0sPhiPt_PCAcut999_TrkDist1(NULL),

  fHistResoK0sMassPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut95_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut95_TrkDist04(NULL),

  fHistResoK0sMassPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut97_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut97_TrkDist04(NULL),

  fHistResoK0sMassPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut997_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut997_TrkDist04(NULL),

  fHistResoK0sMassPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sPtPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sEtaPt_PCAcut999_TrkDist04(NULL),
  fHistResoK0sPhiPt_PCAcut999_TrkDist04(NULL),

  fHistResoK0sMassPt_KstarSel(NULL),
  fHistResoK0sPtPt_KstarSel(NULL),
  fHistResoK0sEtaPt_KstarSel(NULL),
  fHistResoK0sPhiPt_KstarSel(NULL),

  fHistNumK0sCandidates(NULL),

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
  fHistDCAxy(NULL),
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
  fHistV0Armenteros(NULL),

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
  fHistSelV0Armenteros(NULL),

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
  fHistV0Armenteros_TrueK0s(NULL),

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
  fHistSelV0Armenteros_TrueK0s(NULL),
  nEvtAnalyzed(0),
  start(0),
  end(0)
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
   
   //double fCentBins[] = {0,1,5,10,20,40,60,100};
   double fCentBins[] = {-1,5,10,20,40,60,101};
   //double fCentBins[] = {1,2,3,4,6,10,9999};
   //double fVtxBins[]  = {-50, -10.5, -6, -2, 0, 2, 6, 10.5, 50};
   double fVtxBins[]  = {-10.5, -5, 0, 5, 10.5};
   double fPsiBins[]  = {-7,7};

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
   double max_mass = 1.0;
   double width_mass = 0.001;

   double min_pt = 0.0;
   double max_pt = 10.0;
   double width_pt = 0.5;
   
   double min_angle = 0;
   double max_angle = 180.;
   double width_angle = 60.;

   vector<double> bins_cent{0,1,5,10,20,40,60,80,100};
   int binnum_cent = bins_cent.size() - 1;
   
   int bins_K0s[3] = {int((max_mass - min_mass) / width_mass), int((max_pt - min_pt) / width_pt), binnum_cent};
   double min_bins_K0s[3] = {min_mass, min_pt,   0};
   double max_bins_K0s[3] = {max_mass, max_pt, 100};

   fSparseULSPionPairBeforeCuts  = new THnSparseF("fSparseULSPionPairBeforeCuts", "", 3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPairBeforeCuts-> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseULSPionPairBeforeCuts);

   fSparseULSPionPair_PCAcut997_TrkDist1  = 
     new THnSparseF("fSparseULSPionPair_PCAcut997_TrkDist1", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut997_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut997_TrkDist1);
   
   fSparseULSPionPair_PCAcut999_TrkDist1  = 
     new THnSparseF("fSparseULSPionPair_PCAcut999_TrkDist1", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut999_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut999_TrkDist1);

   fSparseULSPionPair_PCAcut95_TrkDist1  = 
     new THnSparseF("fSparseULSPionPair_PCAcut95_TrkDist1", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut95_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut95_TrkDist1);

   fSparseULSPionPair_PCAcut97_TrkDist1  = 
     new THnSparseF("fSparseULSPionPair_PCAcut97_TrkDist1", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut97_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut97_TrkDist1);


   fSparseULSPionPair_PCAcut999_TrkDist04  = 
     new THnSparseF("fSparseULSPionPair_PCAcut999_TrkDist04", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut999_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut999_TrkDist04);

   fSparseULSPionPair_PCAcut997_TrkDist04  = 
     new THnSparseF("fSparseULSPionPair_PCAcut997_TrkDist04", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut997_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut997_TrkDist04);


   fSparseULSPionPair_PCAcut95_TrkDist04  = 
     new THnSparseF("fSparseULSPionPair_PCAcut95_TrkDist04", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut95_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut95_TrkDist04);

   fSparseULSPionPair_PCAcut97_TrkDist04  = 
     new THnSparseF("fSparseULSPionPair_PCAcut97_TrkDist04", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_PCAcut97_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_PCAcut97_TrkDist04);

   fSparseULSPionPair_KstarSel  = 
     new THnSparseF("fSparseULSPionPair_KstarSel", "", 
		    3, bins_K0s,min_bins_K0s, max_bins_K0s);
   fSparseULSPionPair_KstarSel -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseULSPionPair_KstarSel);

      
   min_mass = 0.0;
   max_mass = 6.0;
   width_mass = 0.005;

   min_pt = 0.0;
   max_pt = 20.0;
   width_pt = 0.5;

   min_angle = 0;
   max_angle = 180.;
   width_angle = 60.;
   
   int bins_K0sK0s[4] = {int((max_mass - min_mass) / width_mass), int((max_pt - min_pt) / width_pt), 
			 binnum_cent, int((max_angle - min_angle)/width_angle)};
   double min_bins_K0sK0s[4] = {min_mass, min_pt,   0, min_angle};
   double max_bins_K0sK0s[4] = {max_mass, max_pt, 100, max_angle};

   fSparseNeutralK0sPair_PCAcut997_TrkDist1 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut997_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut997_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut997_TrkDist1);

   fSparseNeutralK0sPair_PCAcut999_TrkDist1 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut999_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut999_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut999_TrkDist1);
   
   fSparseNeutralK0sPair_PCAcut95_TrkDist1 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut95_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut95_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut95_TrkDist1);

   fSparseNeutralK0sPair_PCAcut97_TrkDist1 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut97_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut97_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut97_TrkDist1);

   

   
   fSparseNeutralK0sPair_PCAcut997_TrkDist04 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut997_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut997_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut997_TrkDist04);

   fSparseNeutralK0sPair_PCAcut999_TrkDist04 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut999_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut999_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut999_TrkDist04);

   fSparseNeutralK0sPair_PCAcut95_TrkDist04 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut95_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut95_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut95_TrkDist04);

   fSparseNeutralK0sPair_PCAcut97_TrkDist04 = 
     new THnSparseF("fSparseNeutralK0sPair_PCAcut97_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_PCAcut97_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_PCAcut97_TrkDist04);

   fSparseNeutralK0sPair_KstarSel = 
     new THnSparseF("fSparseNeutralK0sPair_KstarSel", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseNeutralK0sPair_KstarSel -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseNeutralK0sPair_KstarSel);


   fSparseMixNeutralK0sPair_PCAcut997_TrkDist1 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut997_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut997_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut997_TrkDist1);

   fSparseMixNeutralK0sPair_PCAcut999_TrkDist1 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut999_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut999_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut999_TrkDist1);

   fSparseMixNeutralK0sPair_PCAcut95_TrkDist1 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut95_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut95_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut95_TrkDist1);

   fSparseMixNeutralK0sPair_PCAcut97_TrkDist1 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut97_TrkDist1", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut97_TrkDist1 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut97_TrkDist1);
 
   fSparseMixNeutralK0sPair_PCAcut997_TrkDist04 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut997_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut997_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut997_TrkDist04);

   fSparseMixNeutralK0sPair_PCAcut999_TrkDist04 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut999_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut999_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut999_TrkDist04);

   fSparseMixNeutralK0sPair_PCAcut95_TrkDist04 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut95_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut95_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut95_TrkDist04);

   fSparseMixNeutralK0sPair_PCAcut97_TrkDist04 = 
     new THnSparseF("fSparseMixNeutralK0sPair_PCAcut97_TrkDist04", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_PCAcut97_TrkDist04 -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_PCAcut97_TrkDist04);


   fSparseMixNeutralK0sPair_KstarSel = 
     new THnSparseF("fSparseMixNeutralK0sPair_KstarSel", "", 
		    3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);   
   fSparseMixNeutralK0sPair_KstarSel -> SetBinEdges(2,bins_cent.data());   
   fOutputList -> Add(fSparseMixNeutralK0sPair_KstarSel);
   
   fSparseNeutralK0sPairRejectKstar   = new THnSparseF("fSparseNeutralK0sPairRejectKstar", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);  
   fSparseNeutralK0sPairRejectKstar -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseNeutralK0sPairRejectKstar);

   fSparseMixNeutralK0sPairRejectKstar = new THnSparseF("fSparseMixNeutralK0sPairRejectKstar", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fSparseMixNeutralK0sPairRejectKstar -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseMixNeutralK0sPairRejectKstar);
   
   fSparseNeutralK0sPairTransverse = new THnSparseF("fSparseNeutralK0sPairTransverse", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fSparseNeutralK0sPairToward     = new THnSparseF("fSparseNeutralK0sPairToward", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fSparseNeutralK0sPairAway       = new THnSparseF("fSparseNeutralK0sPairAway", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fSparseNeutralK0sPairInJet      = new THnSparseF("fSparseNeutralK0sPairInJet", "", 3, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fSparseNeutralK0sPairTransverse -> SetBinEdges(2,bins_cent.data());
   fSparseNeutralK0sPairToward     -> SetBinEdges(2,bins_cent.data());
   fSparseNeutralK0sPairAway       -> SetBinEdges(2,bins_cent.data());
   fSparseNeutralK0sPairInJet      -> SetBinEdges(2,bins_cent.data());
   fOutputList -> Add(fSparseNeutralK0sPairTransverse);
   fOutputList -> Add(fSparseNeutralK0sPairToward);
   fOutputList -> Add(fSparseNeutralK0sPairAway);
   fOutputList -> Add(fSparseNeutralK0sPairInJet);

   fSparseK0sPion = new THnSparseF("fSparseK0sPion", "", 2, bins_K0sK0s,min_bins_K0sK0s, max_bins_K0sK0s);
   fOutputList -> Add(fSparseK0sPion);

   fSparseULSPionPairRejectKstar = new THnSparseF("fSparseULSPionPairRejectKstar", "", 2, bins_K0s,min_bins_K0s, max_bins_K0s);
   fOutputList -> Add(fSparseULSPionPairRejectKstar);

   if (fIsMC) {     
     int bins_K0s_MC[3]        = {200, 160, 80};
     double min_bins_K0s_MC[3] = {0,  -0.8, 0};
     double max_bins_K0s_MC[3] = {20,  0.8, 8};
     fSparseRecK0s_PCAcut95_TrkDist1   = new THnSparseF("fSparseRecK0s_PCAcut95_TrkDist1", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut97_TrkDist1   = new THnSparseF("fSparseRecK0s_PCAcut97_TrkDist1", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut997_TrkDist1  = new THnSparseF("fSparseRecK0s_PCAcut997_TrkDist1", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut999_TrkDist1  = new THnSparseF("fSparseRecK0s_PCAcut999_TrkDist1", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fOutputList -> Add(fSparseRecK0s_PCAcut95_TrkDist1);
     fOutputList -> Add(fSparseRecK0s_PCAcut97_TrkDist1);
     fOutputList -> Add(fSparseRecK0s_PCAcut997_TrkDist1);
     fOutputList -> Add(fSparseRecK0s_PCAcut999_TrkDist1);

     fSparseRecK0s_PCAcut95_TrkDist04   = new THnSparseF("fSparseRecK0s_PCAcut95_TrkDist04", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut97_TrkDist04   = new THnSparseF("fSparseRecK0s_PCAcut97_TrkDist04", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut997_TrkDist04  = new THnSparseF("fSparseRecK0s_PCAcut997_TrkDist04", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fSparseRecK0s_PCAcut999_TrkDist04  = new THnSparseF("fSparseRecK0s_PCAcut999_TrkDist04", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fOutputList -> Add(fSparseRecK0s_PCAcut95_TrkDist04);
     fOutputList -> Add(fSparseRecK0s_PCAcut97_TrkDist04);
     fOutputList -> Add(fSparseRecK0s_PCAcut997_TrkDist04);
     fOutputList -> Add(fSparseRecK0s_PCAcut999_TrkDist04);

     fSparseRecK0s_KstarSel   = new THnSparseF("fSparseRecK0s_KstarSel", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);
     fOutputList -> Add(fSparseRecK0s_KstarSel);

     fSparseTrueK0s           = new THnSparseF("fSparseTrueK0s", "", 3, bins_K0s_MC,min_bins_K0s_MC, max_bins_K0s_MC);     
     fOutputList    -> Add(fSparseTrueK0s);

     int bins_K0sParams[4]        = {10, 400, 400, 600};
     double min_bins_K0sParams[4] = {0,    0,   0,  -1};
     double max_bins_K0sParams[4] = {10,   2,   2,   2};
     
     fHistResoK0sMassPt_PCAcut95_TrkDist04 = new TH2F("fHistResoK0sMassPt_PCAcut95_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut95_TrkDist04   = new TH2F("fHistResoK0sPtPt_PCAcut95_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut95_TrkDist04  = new TH2F("fHistResoK0sEtaPt_PCAcut95_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut95_TrkDist04  = new TH2F("fHistResoK0sPhiPt_PCAcut95_TrkDist04","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut97_TrkDist04 = new TH2F("fHistResoK0sMassPt_PCAcut97_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut97_TrkDist04   = new TH2F("fHistResoK0sPtPt_PCAcut97_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut97_TrkDist04  = new TH2F("fHistResoK0sEtaPt_PCAcut97_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut97_TrkDist04  = new TH2F("fHistResoK0sPhiPt_PCAcut97_TrkDist04","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut997_TrkDist04 = new TH2F("fHistResoK0sMassPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut997_TrkDist04   = new TH2F("fHistResoK0sPtPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut997_TrkDist04  = new TH2F("fHistResoK0sEtaPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut997_TrkDist04  = new TH2F("fHistResoK0sPhiPt_PCAcut997_TrkDist04","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut997_TrkDist04 = new TH2F("fHistResoK0sMassPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut997_TrkDist04   = new TH2F("fHistResoK0sPtPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut997_TrkDist04  = new TH2F("fHistResoK0sEtaPt_PCAcut997_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut997_TrkDist04  = new TH2F("fHistResoK0sPhiPt_PCAcut997_TrkDist04","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut999_TrkDist04 = new TH2F("fHistResoK0sMassPt_PCAcut999_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut999_TrkDist04   = new TH2F("fHistResoK0sPtPt_PCAcut999_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut999_TrkDist04  = new TH2F("fHistResoK0sEtaPt_PCAcut999_TrkDist04","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut999_TrkDist04  = new TH2F("fHistResoK0sPhiPt_PCAcut999_TrkDist04","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut95_TrkDist1 = new TH2F("fHistResoK0sMassPt_PCAcut95_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut95_TrkDist1   = new TH2F("fHistResoK0sPtPt_PCAcut95_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut95_TrkDist1  = new TH2F("fHistResoK0sEtaPt_PCAcut95_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut95_TrkDist1  = new TH2F("fHistResoK0sPhiPt_PCAcut95_TrkDist1","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut97_TrkDist1 = new TH2F("fHistResoK0sMassPt_PCAcut97_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut97_TrkDist1   = new TH2F("fHistResoK0sPtPt_PCAcut97_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut97_TrkDist1  = new TH2F("fHistResoK0sEtaPt_PCAcut97_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut97_TrkDist1  = new TH2F("fHistResoK0sPhiPt_PCAcut97_TrkDist1","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut997_TrkDist1 = new TH2F("fHistResoK0sMassPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut997_TrkDist1   = new TH2F("fHistResoK0sPtPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut997_TrkDist1  = new TH2F("fHistResoK0sEtaPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut997_TrkDist1  = new TH2F("fHistResoK0sPhiPt_PCAcut997_TrkDist1","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut997_TrkDist1 = new TH2F("fHistResoK0sMassPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut997_TrkDist1   = new TH2F("fHistResoK0sPtPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut997_TrkDist1  = new TH2F("fHistResoK0sEtaPt_PCAcut997_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut997_TrkDist1  = new TH2F("fHistResoK0sPhiPt_PCAcut997_TrkDist1","",200,0,20,5000,0.75,1.25);

     fHistResoK0sMassPt_PCAcut999_TrkDist1 = new TH2F("fHistResoK0sMassPt_PCAcut999_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_PCAcut999_TrkDist1   = new TH2F("fHistResoK0sPtPt_PCAcut999_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_PCAcut999_TrkDist1  = new TH2F("fHistResoK0sEtaPt_PCAcut999_TrkDist1","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_PCAcut999_TrkDist1  = new TH2F("fHistResoK0sPhiPt_PCAcut999_TrkDist1","",200,0,20,5000,0.75,1.25);
     //_KstarSel
     fHistResoK0sMassPt_KstarSel = new TH2F("fHistResoK0sMassPt_KstarSel","",200,0,20,2000,0,2);
     fHistResoK0sPtPt_KstarSel   = new TH2F("fHistResoK0sPtPt_KstarSel","",200,0,20,2000,0,2);
     fHistResoK0sEtaPt_KstarSel  = new TH2F("fHistResoK0sEtaPt_KstarSel","",200,0,20,2000,0,2);
     fHistResoK0sPhiPt_KstarSel  = new TH2F("fHistResoK0sPhiPt_KstarSel","",200,0,20,5000,0.75,1.25);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut95_TrkDist1);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut95_TrkDist1);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut95_TrkDist1);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut95_TrkDist1);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut97_TrkDist1);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut97_TrkDist1);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut97_TrkDist1);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut97_TrkDist1);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut997_TrkDist1);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut997_TrkDist1);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut997_TrkDist1);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut997_TrkDist1);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut999_TrkDist1);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut999_TrkDist1);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut999_TrkDist1);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut999_TrkDist1);


     fOutputList -> Add(fHistResoK0sMassPt_PCAcut95_TrkDist04);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut95_TrkDist04);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut95_TrkDist04);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut95_TrkDist04);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut97_TrkDist04);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut97_TrkDist04);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut97_TrkDist04);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut97_TrkDist04);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut997_TrkDist04);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut997_TrkDist04);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut997_TrkDist04);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut997_TrkDist04);

     fOutputList -> Add(fHistResoK0sMassPt_PCAcut999_TrkDist04);
     fOutputList -> Add(fHistResoK0sPtPt_PCAcut999_TrkDist04);
     fOutputList -> Add(fHistResoK0sEtaPt_PCAcut999_TrkDist04);
     fOutputList -> Add(fHistResoK0sPhiPt_PCAcut999_TrkDist04);

     fOutputList -> Add(fHistResoK0sMassPt_KstarSel);
     fOutputList -> Add(fHistResoK0sPtPt_KstarSel);
     fOutputList -> Add(fHistResoK0sEtaPt_KstarSel);
     fOutputList -> Add(fHistResoK0sPhiPt_KstarSel);
   }   
   
   fHistNumK0sCandidates = new TH1F("fHistNumK0sCandidates","",100,0,100);
   fOutputList -> Add(fHistNumK0sCandidates);

   width_pt = 0.1;
   
   fHistV0PV0DecayLengthXY   = new TH2F("fHistV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
   fHistV0PV0DecayLength     = new TH2F("fHistV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
   fHistV0PV0PointingAngleXY = new TH2F("fHistV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
   fHistV0PV0PointingAngle   = new TH2F("fHistV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
   fHistV0PV0DCAXY           = new TH2F("fHistV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
   fHistV0PV0DCA             = new TH2F("fHistV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
   fHistV0PV0TrackDistanceXY = new TH2F("fHistV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistV0PV0TrackDistance   = new TH2F("fHistV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistV0PV0Chi2perNDF      = new TH2F("fHistV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistV0PV0PropLifeTime    = new TH2F("fHistV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0, 50.);
   fHistV0TrackDCAXY         = new TH2F("fHistV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
   fHistV0TrackDCA           = new TH2F("fHistV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
   fHistV0Armenteros         = new TH2F("fHistV0Armenteros", "",200,-1,1,300,0,0.3);
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
   fOutputList -> Add(fHistV0Armenteros);

   fHistSelV0PV0DecayLengthXY   = new TH2F("fHistSelV0PV0DecayLengthXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
   fHistSelV0PV0DecayLength     = new TH2F("fHistSelV0PV0DecayLength", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
   fHistSelV0PV0PointingAngleXY = new TH2F("fHistSelV0PV0PointingAngleXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
   fHistSelV0PV0PointingAngle   = new TH2F("fHistSelV0PV0PointingAngle", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
   fHistSelV0PV0DCAXY           = new TH2F("fHistSelV0PV0DCAXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
   fHistSelV0PV0DCA             = new TH2F("fHistSelV0PV0DCA", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
   fHistSelV0PV0TrackDistanceXY = new TH2F("fHistSelV0PV0TrackDistanceXY", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistSelV0PV0TrackDistance   = new TH2F("fHistSelV0PV0TrackDistance", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
   fHistSelV0PV0Chi2perNDF      = new TH2F("fHistSelV0PV0Chi2perNDF", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
   fHistSelV0PV0PropLifeTime    = new TH2F("fHistSelV0PV0PropLifeTime", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0, 50.);
   fHistSelV0TrackDCAXY         = new TH2F("fHistSelV0TrackDCAXY", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
   fHistSelV0TrackDCA           = new TH2F("fHistSelV0TrackDCA", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
   fHistSelV0Armenteros         = new TH2F("fHistSelV0Armenteros", "",200,-1,1,300,0,0.3);
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
   fOutputList -> Add(fHistSelV0Armenteros);

   if (fIsMC) {
     fHistV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
     fHistV0PV0DecayLength_TrueK0s     = new TH2F("fHistV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
     fHistV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
     fHistV0PV0PointingAngle_TrueK0s   = new TH2F("fHistV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
     fHistV0PV0DCAXY_TrueK0s           = new TH2F("fHistV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
     fHistV0PV0DCA_TrueK0s             = new TH2F("fHistV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
     fHistV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistV0PV0TrackDistance_TrueK0s   = new TH2F("fHistV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0, 50.);
     fHistV0TrackDCAXY_TrueK0s         = new TH2F("fHistV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
     fHistV0TrackDCA_TrueK0s           = new TH2F("fHistV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
     fHistV0Armenteros_TrueK0s         = new TH2F("fHistV0Armenteros_TrueK0s", "",200,-1,1,300,0,0.3);
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
     fOutputList -> Add(fHistV0Armenteros_TrueK0s);

     fHistSelV0PV0DecayLengthXY_TrueK0s   = new TH2F("fHistSelV0PV0DecayLengthXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
     fHistSelV0PV0DecayLength_TrueK0s     = new TH2F("fHistSelV0PV0DecayLength_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0., 100.);
     fHistSelV0PV0PointingAngleXY_TrueK0s = new TH2F("fHistSelV0PV0PointingAngleXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
     fHistSelV0PV0PointingAngle_TrueK0s   = new TH2F("fHistSelV0PV0PointingAngle_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,10000, 0.99, 1.0);
     fHistSelV0PV0DCAXY_TrueK0s           = new TH2F("fHistSelV0PV0DCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
     fHistSelV0PV0DCA_TrueK0s             = new TH2F("fHistSelV0PV0DCA_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 1000, 0., 10.);
     fHistSelV0PV0TrackDistanceXY_TrueK0s = new TH2F("fHistSelV0PV0TrackDistanceXY_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistSelV0PV0TrackDistance_TrueK0s   = new TH2F("fHistSelV0PV0TrackDistance_TrueK0s", "",int((max_pt - min_pt) / width_pt),min_pt, max_pt, 10000, 0, 10);
     fHistSelV0PV0Chi2perNDF_TrueK0s      = new TH2F("fHistSelV0PV0Chi2perNDF_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,1000, 0., 10.);
     fHistSelV0PV0PropLifeTime_TrueK0s    = new TH2F("fHistSelV0PV0PropLifeTime_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt,500, 0, 50.);
     fHistSelV0TrackDCAXY_TrueK0s         = new TH2F("fHistSelV0TrackDCAXY_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
     fHistSelV0TrackDCA_TrueK0s           = new TH2F("fHistSelV0TrackDCA_TrueK0s", "",int((max_pt - min_pt) / width_pt), min_pt, max_pt, 5000, 0., 50.);
     fHistSelV0Armenteros_TrueK0s         = new TH2F("fHistSelV0Armenteros_TrueK0s", "",200,-1,1,300,0,0.3);
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
     fOutputList -> Add(fHistSelV0Armenteros_TrueK0s);
   }
   
   fHistTrackP                    = new TH1F("fHistTrackP", "", 100, 0, 10);
   fHistTrackPt                   = new TH1F("fHistTrackPt", "", 100, 0, 10);
   fHistTrackEta                  = new TH1F("fHistTrackEta", "", 20, -1, 1);   
   fHistTPCNClusts                = new TH2F("fHistTPCNClusts", "", 160, 0, 160, 100, 0, 10);
   fHistSPDNClusts                = new TH2F("fHistSPDNClusts", "", 6, 0, 6, 100, 0, 10);
   fHistTPCCrossRowsFindableRatio = new TH2F("fHistTPCCrossRowsFindableRatio", "", 200, 0, 2, 100, 0, 10);
   fHistReducedChi2TPC            = new TH2F("fHistReducedChi2TPC", "", 100, 0, 10, 100, 0, 10);
   fHistReducedChi2ITS            = new TH2F("fHistReducedChi2ITS", "", 250, 0, 50, 100, 0, 10);
   fHistDCAz                      = new TH2F("fHistDCAz", "", 1200, -60, 60, 100, 0, 10);
   fHistDCAxy                     = new TH2F("fHistDCAxy", "", 2000, 0, 20, 100, 0, 10);
   
   fOutputList->Add(fHistTrackP);
   fOutputList->Add(fHistTrackPt);
   fOutputList->Add(fHistTrackEta);
   fOutputList->Add(fHistTPCNClusts);
   fOutputList->Add(fHistSPDNClusts);
   fOutputList->Add(fHistTPCCrossRowsFindableRatio);
   fOutputList->Add(fHistReducedChi2TPC);
   fOutputList->Add(fHistReducedChi2ITS);
   fOutputList->Add(fHistDCAz);
   fOutputList->Add(fHistDCAxy);
   
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
     
   fArrayK0s              = new TObjArray();      
   fArrayK0sDaughterTrack = new TObjArray();
   fArrayPrimaryPionTrack = new TObjArray();
   
   if (!Initialize())            return ;
   if (!fUtils->isAcceptEvent()) return ;
   
   //checkPrimaryTrack();

   AddK0sArray(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));
   fHistNumK0sCandidates->Fill(fNK0s);

   if (fNK0s>1){
     K0sPairAnalysis(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));     
     K0sPairAnalysisEventMixing(fUtils->getPairTargetPIDs(0),fUtils->getPairTargetPIDs(1));
   }

   if(fIsMC) ProcessMC();

   delete fArrayK0s;
   delete fArrayK0sDaughterTrack;
   delete fArrayPrimaryPionTrack;
   
   if (nEvtAnalyzed == 0){
     start = clock();
   }else if (nEvtAnalyzed % 1000 == 0){
     end = clock();
     const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;   
     cout<<Form("(%d) time %lf[s]\n", nEvtAnalyzed, time)<<"   AOD file ===>  "<<fInputHandler->GetInputFileName()<<endl;     
     start = clock();
   } 
   
   nEvtAnalyzed++;

 }

 bool AliAnalysisTaskAODTrackPair::Initialize() {
   
   fNK0s = 0;

   fEvent = dynamic_cast<AliAODEvent *>(InputEvent());
   
   if (!fUtils->setEvent(fEvent, fInputHandler)) return false;
   
   fPIDResponse = fInputHandler->GetPIDResponse();
   if ( !fPIDResponse ) return false; 

   fPrimVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();
   
   if (!fPrimVtx)                                       return false;   
   if (fPrimVtx->GetNContributors() < fMinNContPrimVtx) return false;
   
   fPrimVtx->GetXYZ(fPrimVtxPos);
   fPrimVtx->GetCovarianceMatrix(fPrimVtxCov);
   
   if (fMinVertexCutZ > fPrimVtxPos[2] || fMaxVertexCutZ < fPrimVtxPos[2]) return false;

   fMultSelection = (AliMultSelection *)fEvent->FindListObject("MultSelection");
   
   if (!fMultSelection) return false;

   fCent = fMultSelection->GetMultiplicityPercentile(fMethodCent.c_str(), false);
   
   if (0>fCent || fCent>100) return false;
   
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
  double beta   = sigTOF > 0 ? length / (2.99792457999974e-02 * sigTOF) : -999;
  double dEdx   = track->GetTPCsignal();
    
  if (isSel) {
    fHistSelTPCdEdxP         -> Fill(p, dEdx);
    fHistSelTPCSigmaElectron -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
    fHistSelTPCSigmaMuon     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
    fHistSelTPCSigmaPion     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    fHistSelTPCSigmaKaon     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
    fHistSelTPCSigmaProton   -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    if (track->GetTOFBunchCrossing() == 0 && beta > 0.) {
      fHistSelBetaP            -> Fill(p, beta);
      fHistSelTOFSigmaElectron -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
      fHistSelTOFSigmaMuon     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
      fHistSelTOFSigmaPion     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
      fHistSelTOFSigmaKaon     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
      fHistSelTOFSigmaProton   -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }
  } else {
    fHistTPCdEdxP         -> Fill(p, dEdx);
    fHistTPCSigmaElectron -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
    fHistTPCSigmaMuon     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
    fHistTPCSigmaPion     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    fHistTPCSigmaKaon     -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
    fHistTPCSigmaProton   -> Fill(p, fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    if (track->GetTOFBunchCrossing() == 0 && beta > 0.) {
      fHistBetaP            -> Fill(p, beta);
      fHistTOFSigmaElectron -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
      fHistTOFSigmaMuon     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
      fHistTOFSigmaPion     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
      fHistTOFSigmaKaon     -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
      fHistTOFSigmaProton   -> Fill(p, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::TrackQualityChecker(AliAODTrack *track) {

  fHistTPCNClusts->Fill(track->GetTPCNcls(), track->P());

  int nSPD = 0;

  if (track->HasPointOnITSLayer(0)) ++nSPD;
  if (track->HasPointOnITSLayer(1)) ++nSPD;

  fHistSPDNClusts->Fill(nSPD, track->P());
  
  double fr       = track->GetTPCCrossedRows() > 0 ? track->GetTPCNclsF() / track->GetTPCCrossedRows() : 0;
  double rchi2tpc = track->GetTPCNcls() > 0 ? track->GetTPCchi2() / track->GetTPCNcls() : 0;
  double rchi2its = track->GetITSNcls() > 0 ? track->GetITSchi2() / track->GetITSNcls() : 0;
  
  fHistTPCCrossRowsFindableRatio -> Fill(fr, track->P());
  fHistReducedChi2TPC            -> Fill(rchi2tpc, track->P());
  fHistReducedChi2ITS            -> Fill(rchi2its, track->P());
  
  double maxd   = 999.;
  double dz[2]  = {};
  double cov[3] = {};
  track->PropagateToDCA(fPrimVtx,fEvent->GetMagneticField(),maxd,dz,cov);

  KFParticle kfTrack;
  kfTrack = CreateKFParticle(track, 211);
  float dca, e_dca;
  float fPrimVtxPosF[] = {static_cast<float>(fPrimVtxPos[0]),static_cast<float>(fPrimVtxPos[1]),static_cast<float>(fPrimVtxPos[2])};
  float fPrimVtxCovF[] = {static_cast<float>(fPrimVtxCov[0]),static_cast<float>(fPrimVtxCov[1]),static_cast<float>(fPrimVtxCov[2]),
			  static_cast<float>(fPrimVtxCov[3]),static_cast<float>(fPrimVtxCov[4]),static_cast<float>(fPrimVtxCov[5])};
  kfTrack.GetDistanceFromVertexXY(fPrimVtxPosF,fPrimVtxCovF,dca,e_dca);
  double dca_xy = dca;
  double dca_z  = dz[1];
    
  fHistDCAz     -> Fill(dca_z, track->P());
  fHistDCAxy    -> Fill(dca_xy, track->P());  
  fHistTrackP   -> Fill(track->P());
  fHistTrackPt  -> Fill(track->Pt());
  fHistTrackEta -> Fill(track->Eta());

  return true;
}

bool AliAnalysisTaskAODTrackPair::V0QualityChecker(double mass, double p, int charge, double DCApairXY, double DCApair, double decay_length, 
						   double decay_length_xy, double cpv,
						   double cpv_xy, double dca, double dca_xy, double distance_daughter, double distance_daughter_xy,
						   double chi2, double lifetime, double armenteros_alpha, double armenteros_qt, int dlabel[2],bool isSel){  
  
  if ( !isAcceptK0sCandidateMassRange(mass) ) return false;
  
  if (isSel) {    
    
    fHistSelV0PV0DecayLengthXY   -> Fill(p, decay_length_xy);
    fHistSelV0PV0DecayLength     -> Fill(p, decay_length);
    fHistSelV0PV0PointingAngleXY -> Fill(p, cpv_xy);
    fHistSelV0PV0PointingAngle   -> Fill(p, cpv);
    fHistSelV0PV0DCAXY           -> Fill(p, dca_xy);
    fHistSelV0PV0DCA             -> Fill(p, dca);
    fHistSelV0PV0TrackDistanceXY -> Fill(p, distance_daughter_xy);
    fHistSelV0PV0TrackDistance   -> Fill(p, distance_daughter);
    fHistSelV0PV0Chi2perNDF      -> Fill(p, chi2);
    fHistSelV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * decay_length / p);    
    fHistSelV0TrackDCAXY         -> Fill(p, DCApairXY);
    fHistSelV0TrackDCA           -> Fill(p, DCApair);
    fHistSelV0Armenteros         -> Fill(armenteros_alpha,armenteros_qt);

    if(fIsMC && 
       fUtils->isSameMotherPair(dlabel[0],dlabel[1]) &&
       fUtils->getMotherPdgCode(dlabel[0]) == 310 &&
       fUtils->getGrandMotherLabel(dlabel[0]) == -1){      
    
      fHistSelV0PV0DecayLengthXY_TrueK0s   -> Fill(p, decay_length_xy);
      fHistSelV0PV0DecayLength_TrueK0s     -> Fill(p, decay_length);
      fHistSelV0PV0PointingAngleXY_TrueK0s -> Fill(p, cpv_xy);
      fHistSelV0PV0PointingAngle_TrueK0s   -> Fill(p, cpv);
      fHistSelV0PV0DCAXY_TrueK0s           -> Fill(p, dca_xy);
      fHistSelV0PV0DCA_TrueK0s             -> Fill(p, dca);
      fHistSelV0PV0TrackDistanceXY_TrueK0s -> Fill(p, distance_daughter_xy);
      fHistSelV0PV0TrackDistance_TrueK0s   -> Fill(p, distance_daughter);
      fHistSelV0PV0Chi2perNDF_TrueK0s      -> Fill(p, chi2);
      fHistSelV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * decay_length / p);    
      fHistSelV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
      fHistSelV0TrackDCA_TrueK0s           -> Fill(p, DCApair);      
      fHistSelV0Armenteros_TrueK0s         -> Fill(armenteros_alpha,armenteros_qt);

    }    
  } else {
    
    fHistV0PV0DecayLengthXY   -> Fill(p, decay_length_xy);
    fHistV0PV0DecayLength     -> Fill(p, decay_length);
    fHistV0PV0PointingAngleXY -> Fill(p, cpv_xy);
    fHistV0PV0PointingAngle   -> Fill(p, cpv);
    fHistV0PV0DCAXY           -> Fill(p, dca_xy);
    fHistV0PV0DCA             -> Fill(p, dca);
    fHistV0PV0TrackDistanceXY -> Fill(p, distance_daughter_xy);
    fHistV0PV0TrackDistance   -> Fill(p, distance_daughter);
    fHistV0PV0Chi2perNDF      -> Fill(p, chi2);
    fHistV0PV0PropLifeTime    -> Fill(p, fUtils->fPdgK0sMass * decay_length / p);    
    fHistV0TrackDCAXY         -> Fill(p, DCApairXY);
    fHistV0TrackDCA           -> Fill(p, DCApair);
    fHistV0Armenteros         -> Fill(armenteros_alpha,armenteros_qt);

    if(fIsMC && 
       fUtils->isSameMotherPair(dlabel[0],dlabel[1]) &&
       fUtils->getMotherPdgCode(dlabel[0]) == 310 &&
       fUtils->getGrandMotherLabel(dlabel[0]) == -1){      
    
      fHistV0PV0DecayLengthXY_TrueK0s   -> Fill(p, decay_length_xy);
      fHistV0PV0DecayLength_TrueK0s     -> Fill(p, decay_length);
      fHistV0PV0PointingAngleXY_TrueK0s -> Fill(p, cpv_xy);
      fHistV0PV0PointingAngle_TrueK0s   -> Fill(p, cpv);
      fHistV0PV0DCAXY_TrueK0s           -> Fill(p, dca_xy);
      fHistV0PV0DCA_TrueK0s             -> Fill(p, dca);
      fHistV0PV0TrackDistanceXY_TrueK0s -> Fill(p, distance_daughter_xy);
      fHistV0PV0TrackDistance_TrueK0s   -> Fill(p, distance_daughter);
      fHistV0PV0Chi2perNDF_TrueK0s      -> Fill(p, chi2);
      fHistV0PV0PropLifeTime_TrueK0s    -> Fill(p, fUtils->fPdgK0sMass * decay_length / p);    
      fHistV0TrackDCAXY_TrueK0s         -> Fill(p, DCApairXY);
      fHistV0TrackDCA_TrueK0s           -> Fill(p, DCApair);      
      fHistV0Armenteros_TrueK0s         -> Fill(armenteros_alpha,armenteros_qt);

    }    
  }
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::V0QualityChecker(AliAODv0* v0, bool isSel) {
  
  double p  = v0->P();
  
  double DCApairXY = v0->DcaPosToPrimVertex()*v0->DcaPosToPrimVertex() + v0->DcaNegToPrimVertex()*v0->DcaNegToPrimVertex() ? 
    sqrt(v0->DcaPosToPrimVertex()*v0->DcaPosToPrimVertex() + v0->DcaNegToPrimVertex()*v0->DcaNegToPrimVertex()) : 0;
  double DCApair   = DCApairXY;
  
  int charge                  = v0->Charge();
  double mass                 = v0->MassK0Short();
  double decay_length         = v0->DecayLength(fPrimVtxPos);
  double decay_length_xy      = v0->DecayLengthXY(fPrimVtxPos);
  double cpv                  = v0->CosPointingAngle(fPrimVtxPos);
  double cpv_xy               = v0->CosPointingAngleXY(fPrimVtxPos);
  double dca                  = v0->DcaV0ToPrimVertex ();
  double dca_xy               = 0;
  double distance_daughter    = v0->DcaV0Daughters();
  double distance_daughter_xy = 0;
  double chi2                 = v0->Chi2V0();
  double lifetime             = fUtils->fPdgK0sMass * decay_length / p;
  double armenteros_alpha     = v0->AlphaV0();
  double armenteros_qt        = v0->PtArmV0();
    
  int dlabel[2] = {};
  
  AliAODTrack *aodTrack1 = (AliAODTrack *)v0->GetDaughter(0);
  AliAODTrack *aodTrack2 = (AliAODTrack *)v0->GetDaughter(1);

  if (!aodTrack1 || !aodTrack2) {
    dlabel[0] = -1;
    dlabel[1] = -1;
  } else {
    dlabel[0] = aodTrack1->GetLabel();
    dlabel[1] = aodTrack2->GetLabel();
  }

  return V0QualityChecker(mass, p, charge, DCApairXY, DCApair, decay_length, decay_length_xy, cpv,
			  cpv_xy, dca, dca_xy, distance_daughter, distance_daughter_xy,
			  chi2, lifetime, armenteros_alpha, armenteros_qt,dlabel,isSel);
}

bool AliAnalysisTaskAODTrackPair::K0sPairAnalysis(AliPID::EParticleType pid1,AliPID::EParticleType pid2) {  
  
  bool existLeading = false;
  
  unique_ptr<TVector3> leading(new TVector3(0.,0.,0.));
  if (fLeadingTrack) {
    leading->SetXYZ(fLeadingTrack->Px(),fLeadingTrack->Py(),fLeadingTrack->Pz()); 
    existLeading = true;
  }

  for (int iV0_1 = 0; iV0_1 < fNK0s; ++iV0_1) {
    
    AliAODv0* v0_1 = static_cast<AliAODv0*>(fArrayK0s->At(iV0_1));
    
    unique_ptr<TLorentzVector> lv1(new TLorentzVector());    
    lv1->SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(), massK0s);
    
    for (int iV0_2 = iV0_1 + 1; iV0_2 < fNK0s; ++iV0_2) {
      
      AliAODv0* v0_2 = static_cast<AliAODv0*>(fArrayK0s->At(iV0_2));
      
      unique_ptr<TLorentzVector> lv2(new TLorentzVector());    
      lv2->SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(), massK0s);
      
      if (!isAcceptK0sPair(v0_1,v0_2)) continue; 
      
      unique_ptr<TLorentzVector> lv12(new TLorentzVector(lv1->Vect()+lv2->Vect(),lv1->E()+lv2->E()));
      
      double angle = angle = leading->Angle(lv12->Vect()) * 180. / TMath::Pi();
      
      double fill[] = {lv12->M(), lv12->Pt(), fCent};      
      
      double cpa_1      = v0_1->CosPointingAngle(fPrimVtxPos);
      double trkDist_1  = v0_1->DcaV0Daughters();
      bool isKstarSel_1 = false;
      double *dummy_1 = v0_1->GetPID();      
      if (dummy_1[0]>1.0) {
	isKstarSel_1 = true;
      }      
      double cpa_2      = v0_2->CosPointingAngle(fPrimVtxPos);
      double trkDist_2  = v0_2->DcaV0Daughters();
      bool isKstarSel_2 = false;
      double *dummy_2 = v0_2->GetPID();      
      if (dummy_2[0]>1.0) {
	isKstarSel_2 = true;
      } 
      if (isKstarSel_1 && isKstarSel_2) {
	fSparseNeutralK0sPair_KstarSel->Fill(fill);
      }
      if (trkDist_1 < 1.0 && trkDist_2 < 1.0) {
	if (cpa_1>0.997 && cpa_2>0.997)
	  fSparseNeutralK0sPair_PCAcut997_TrkDist1   -> Fill(fill);
	if (cpa_1>0.999 && cpa_2>0.999)
	  fSparseNeutralK0sPair_PCAcut999_TrkDist1   -> Fill(fill);
	if (cpa_1>0.97 && cpa_2>0.97)
	  fSparseNeutralK0sPair_PCAcut97_TrkDist1    -> Fill(fill);
	if (cpa_1>0.95 && cpa_2>0.95)
	  fSparseNeutralK0sPair_PCAcut95_TrkDist1    -> Fill(fill);
      }       
      if (trkDist_1 < 0.4 && trkDist_2 < 0.4) {
	if (cpa_1>0.997 && cpa_2>0.997)
	  fSparseNeutralK0sPair_PCAcut997_TrkDist04   -> Fill(fill);
	if (cpa_1>0.999 && cpa_2>0.999)
	  fSparseNeutralK0sPair_PCAcut999_TrkDist04   -> Fill(fill);
	if (cpa_1>0.97 && cpa_2>0.97)
	  fSparseNeutralK0sPair_PCAcut97_TrkDist04    -> Fill(fill);
	if (cpa_1>0.95 && cpa_2>0.95)
	  fSparseNeutralK0sPair_PCAcut95_TrkDist04    -> Fill(fill);
      } 
	
      //GetOnFlyStatus is dummy flag indicating K*+ decay K0s	
      if ( v0_1->GetOnFlyStatus() == false && v0_2->GetOnFlyStatus() == false )
	fSparseNeutralK0sPairRejectKstar->Fill(fill);

      if (existLeading) {
	if ( 0.<angle && angle<60. ) {
	  fSparseNeutralK0sPairToward -> Fill(fill);
	  fSparseNeutralK0sPairInJet  -> Fill(fill);
	} else if ( 60.<angle && angle<120. ){
	  fSparseNeutralK0sPairTransverse -> Fill(fill);
	} else {
	  fSparseNeutralK0sPairAway  -> Fill(fill);
	  fSparseNeutralK0sPairInJet -> Fill(fill);
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
  double poolPsi  = 0.;
  
  if (onEvtMixingPoolVtxZ) poolVtxZ = fPrimVtxPos[2];
  if (onEvtMixingPoolCent) poolCent = fCent;
  if (onEvtMixingPoolPsi)  poolPsi  = 0;

  poolCent = 0.;
  poolVtxZ = 0.;
  poolPsi  = 0.;

  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgrK0s->GetEventPool(poolCent, poolVtxZ, poolPsi);
  
  for (int iTrack1 = 0; iTrack1 < fNK0s; ++iTrack1) {
    
    AliAODv0* v0_1 = static_cast<AliAODv0*>(fArrayK0s->At(iTrack1));
    if ( !v0_1 )  continue;
    
    unique_ptr<TLorentzVector> lv1(new TLorentzVector());    
    lv1->SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(), massK0s);
    
    if (pool->IsReady()) {
      
      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {
	
	TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

	for (int iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast(); ++iTrack2) {
	  
	  AliAODv0* v0_2 = dynamic_cast<AliAODv0*>(poolTracks->At(iTrack2));
	  if ( !v0_2 ) continue;
	  
	  unique_ptr<TLorentzVector> lv2(new TLorentzVector());    
	  lv2->SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(), massK0s);

	  unique_ptr<TLorentzVector> lv12(new TLorentzVector(lv1->Vect()+lv2->Vect(),lv1->E()+lv2->E()));

	  double fill[] = {lv12->M(), lv12->Pt(), fCent};
	  
	  double cpa_1 = v0_1->CosPointingAngle(fPrimVtxPos);
	  double trkDist_1  = v0_1->DcaV0Daughters();
	  bool isKstarSel_1 = false;	  
	  double *dummy_1 = v0_1->GetPID();      
	  if (dummy_1[0]>1.0) {
	    isKstarSel_1 = true;
	  }
	  double cpa_2 = 0;	  	  	 
	  if (v0_2->GetCharge() == 0)      cpa_2 = 0.99991; 
	  else if (v0_2->GetCharge() == 1) cpa_2 = 0.9991; 
	  else if (v0_2->GetCharge() == 2) cpa_2 = 0.9971; 
	  else if (v0_2->GetCharge() == 3) cpa_2 = 0.971; 
	  else if (v0_2->GetCharge() == 4) cpa_2 = 0.951; 
	  double trkDist_2  = v0_2->DcaV0Daughters();	  	  
	  bool isKstarSel_2 = false;
	  double *dummy_2 = v0_1->GetPID();      
	  if (dummy_2[0]>1.0) {
	    isKstarSel_2 = true;
	  }
	  if (isKstarSel_1 && isKstarSel_2) {
	    fSparseMixNeutralK0sPair_KstarSel->Fill(fill);
	  }
	  if (trkDist_1 < 1.0 && trkDist_2 < 1.0) {
	    if (cpa_1>0.997 && cpa_2>0.997)
	      fSparseMixNeutralK0sPair_PCAcut997_TrkDist1   -> Fill(fill);
	    if (cpa_1>0.999 && cpa_2>0.999)
	      fSparseMixNeutralK0sPair_PCAcut999_TrkDist1   -> Fill(fill);
	    if (cpa_1>0.97 && cpa_2>0.97)
	      fSparseMixNeutralK0sPair_PCAcut97_TrkDist1    -> Fill(fill);
	    if (cpa_1>0.95 && cpa_2>0.95)
	      fSparseMixNeutralK0sPair_PCAcut95_TrkDist1    -> Fill(fill);
	  } 
	  if (trkDist_1 < 0.4 && trkDist_2 < 0.4) {
	    if (cpa_1>0.997 && cpa_2>0.997)
	      fSparseMixNeutralK0sPair_PCAcut997_TrkDist04   -> Fill(fill);
	    if (cpa_1>0.999 && cpa_2>0.999)
	      fSparseMixNeutralK0sPair_PCAcut999_TrkDist04   -> Fill(fill);
	    if (cpa_1>0.97 && cpa_2>0.97)
	      fSparseMixNeutralK0sPair_PCAcut97_TrkDist04 -> Fill(fill);
	    if (cpa_1>0.95 && cpa_2>0.95)
	      fSparseMixNeutralK0sPair_PCAcut95_TrkDist04 -> Fill(fill);
	  } 

	  if ( v0_1->GetOnFlyStatus() == false && v0_2->GetOnFlyStatus() == false )
	    fSparseMixNeutralK0sPairRejectKstar->Fill(fill);
	    
	}//end of loop iTrack2
      }//end of loop iMixEvt	  	
    }//end of pool->IsReady()  
    
    fTrackArray->Add(v0_1);
  }//end of loop iTrack1
  
  //TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  TObjArray *fTrackArrayClone = dynamic_cast<TObjArray*>(fTrackArray->Clone());  
  fTrackArrayClone->SetOwner(true);
  pool->UpdatePool(fTrackArrayClone);
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::AddK0sArray(AliPID::EParticleType pid1,AliPID::EParticleType pid2) {  
  
  TObjArray* fTrackArray = new TObjArray();
  fTrackArray->SetOwner(true);
    
  double poolCent = 0.;
  double poolVtxZ = 0.;
  double poolPsi  = 0.;
  
  if (onEvtMixingPoolVtxZ) poolVtxZ = fPrimVtxPos[2];
  if (onEvtMixingPoolCent) poolCent = fCent;
  if (onEvtMixingPoolPsi)  poolPsi  = fPsi;
  
  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgrPion->GetEventPool(poolCent, poolVtxZ, poolPsi);
  
  int nV0 = fEvent->GetNumberOfV0s();

  for (int iv0=0; iv0<nV0; ++ iv0) {

    AliAODv0 *aodv0 = (AliAODv0 *)fEvent->GetV0(iv0);      

    if ( aodv0->GetOnFlyStatus() )                             continue;      
    if ( !isAcceptK0sKinematics(aodv0) )                       continue;
    if ( 0.4>aodv0->MassK0Short() || aodv0->MassK0Short()>0.6) continue;

    AliAODTrack *aodTrack1 = (AliAODTrack *)aodv0->GetDaughter(0);
    AliAODTrack *aodTrack2 = (AliAODTrack *)aodv0->GetDaughter(1);
    
    if ( !isAcceptDaughterTrack(aodTrack1,pid1) || !isAcceptDaughterTrack(aodTrack2,pid2)) continue;           
    if ( !isAcceptDaughterTrackPairAngle(aodTrack1,aodTrack2) )                            continue;      
    
    AliAODv0 *v0 = aodv0;//updateAODv0(aodTrack1,aodTrack2);
    if (!v0) continue;

    int charge = v0->Charge();

    V0QualityChecker(v0,false);

    double fill[] = {v0->MassK0Short(), v0->Pt(), fCent};

    if ( !isAcceptK0sKinematics(v0) ) continue;

    //Set K* -> K0s + pion candidate      
    bool isFromKstar = isKstarCandidate(v0);      
    v0->SetOnFlyStatus(isFromKstar);
    fSparseULSPionPairBeforeCuts  -> Fill(fill);

    bool isAcceptCuts[7];          
    
    V0QualityChecker(v0,true);

    double chi2     = v0->Chi2V0();    
    double lengthXY = v0->DecayLengthXY(fPrimVtxPos);
    double DCAV0    = v0->DcaV0ToPrimVertex();   
    double cpa      = v0->CosPointingAngle(fPrimVtxPos);
    double trkDist  = v0->DcaV0Daughters();
    double radius   = v0->RadiusV0();
    double lifetime = fUtils->fPdgK0sMass * v0->DecayLength(fPrimVtxPos) / v0->P();
    
    double isKstarSel = 0.;
    if (lifetime<20 && lengthXY>0.5 && cpa>0.97 && DCAV0<0.3 && fabs(v0->RapK0Short())<0.8 && 
	fabs(v0->MassLambda()-1.115683)>0.0043 && fabs(v0->MassAntiLambda()-1.115683)>0.0043 && trkDist<1.0 ) {
      isKstarSel = 2.;
    }    
    if (isKstarSel>1.){
      fSparseULSPionPair_KstarSel -> Fill(fill);
    }
    if (trkDist < 1.0) {
      if (cpa>0.999)
	fSparseULSPionPair_PCAcut999_TrkDist1   -> Fill(fill);
      if (cpa>0.997)
	fSparseULSPionPair_PCAcut997_TrkDist1   -> Fill(fill);
      if (cpa>0.97)
	fSparseULSPionPair_PCAcut97_TrkDist1    -> Fill(fill);
      if (cpa>0.95)
	fSparseULSPionPair_PCAcut95_TrkDist1    -> Fill(fill);      
    }       
    if (trkDist < 0.4) {
      if (cpa>0.999)
	fSparseULSPionPair_PCAcut999_TrkDist04   -> Fill(fill);
      if (cpa>0.997)
	fSparseULSPionPair_PCAcut997_TrkDist04   -> Fill(fill);
      if (cpa>0.97)
	fSparseULSPionPair_PCAcut97_TrkDist04    -> Fill(fill);
      if (cpa>0.95)
	fSparseULSPionPair_PCAcut95_TrkDist04    -> Fill(fill);
    } 
    
    if (!isFromKstar) fSparseULSPionPairRejectKstar -> Fill(fill);      
    
    if(isAcceptK0sCandidateMassRange(v0->MassK0Short())) {			
      short int _cpa_ = 0;      
      if (cpa>0.9999) _cpa_ = 0;
      else if (cpa>0.999)  _cpa_ = 1;
      else if (cpa>0.997)  _cpa_ = 2;
      else if (cpa>0.97)   _cpa_ = 3;
      else if (cpa>0.95)   _cpa_ = 4;      
      v0->SetCharge(_cpa_);
      double dummy[]={isKstarSel,isKstarSel};
      v0->SetPID(2,dummy);
      fArrayK0s->Add(v0);
    }
    
  }
  
  fNK0s = fArrayK0s->GetEntriesFast();
  
  return true;  
}

AliAODv0 *AliAnalysisTaskAODTrackPair::updateAODv0(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2){
  
  bool isSignalFlag = false;
  
  double MagF = fEvent->GetMagneticField();
  
  KFPVertex kfVtxPrimeVertex;
  double fPrimVtxCov[6]={};
  fPrimVtx->GetCovarianceMatrix(fPrimVtxCov);

  float fPrimVtxPosF[3],fPrimVtxCovF[6];    
  for(int iEl = 0; iEl < 3; iEl++) fPrimVtxPosF[iEl] = (float)fPrimVtxPos[iEl];
  for(int iEl = 0; iEl < 6; iEl++) fPrimVtxCovF[iEl] = (float)fPrimVtxCov[iEl];  

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
  /*
  unique_ptr<KFParticleK0s> kfCheckIfSafePair(new KFParticleK0s());
  kfCheckIfSafePair->AddDaughter(kfTrack[1]);
  
  if (!kfCheckIfSafePair->CheckIfSafePair(kfTrack[0])) {
    return NULL;
  }
  */
  KFParticle kfK0s;
  kfK0s.Initialize();
  kfK0s.SetConstructMethod(0); 
  
  const KFParticle *kfKs0Daughters[2] = {&kfTrack[0],&kfTrack[1]};  
  kfK0s.Construct(kfKs0Daughters, 2);
  
  float chi2K0s   = kfK0s.GetChi2() / kfK0s.GetNDF();  
  if (chi2K0s > 1e+2) return NULL;

  float momK0s[3] = {kfK0s.GetPx(),kfK0s.GetPy(),kfK0s.GetPz()};
  
  float posK0s[3] = {};
  float covK0s[36]= {};
  posK0s[0] = kfK0s.GetX();
  posK0s[1] = kfK0s.GetY();
  posK0s[2] = kfK0s.GetZ();
  
  //Propagate to daughetrs to the secondary (their production) vertex
  kfTrack[0].TransportToPoint(posK0s);
  kfTrack[1].TransportToPoint(posK0s);
  
  double rDcaV0Daughters    = kfTrack[0].GetDistanceFromParticleXY(kfTrack[1]);
  double rDcaV0ToPrimVertex = DecayLengthFromKF(kfK0s,kfParticlePrimeVertex);
  const double rMomPos[3]={kfTrack[0].Px(),kfTrack[0].Py(),kfTrack[0].Pz()};
  const double rMomNeg[3]={kfTrack[1].Px(),kfTrack[1].Py(),kfTrack[1].Pz()};
  
  double rDcaDaughterToPrimVertex[2] = {};   
  float dcaDaughter_xy1, e_dcaDaughter_xy1;
  float dcaDaughter_xy2, e_dcaDaughter_xy2;  
  kfTrack[0].GetDistanceFromVertexXY(kfParticlePrimeVertex,dcaDaughter_xy1,e_dcaDaughter_xy1);
  kfTrack[1].GetDistanceFromVertexXY(kfParticlePrimeVertex,dcaDaughter_xy2,e_dcaDaughter_xy2);  
  rDcaDaughterToPrimVertex[0] = static_cast<double>(dcaDaughter_xy1);
  rDcaDaughterToPrimVertex[1] = static_cast<double>(dcaDaughter_xy2);
  
  AliAODVertex *vertex = new AliAODVertex();
  vertex->SetX( kfK0s.GetX() );
  vertex->SetY( kfK0s.GetY() );
  vertex->SetZ( kfK0s.GetZ() );
  vertex->SetChi2perNDF(static_cast<double>(kfK0s.GetChi2() / kfK0s.GetNDF()));
  vertex->AddDaughter(aodTrack1);
  vertex->AddDaughter(aodTrack2);
  
  AliAODv0* updated_v0 = new AliAODv0(vertex,rDcaV0Daughters,rDcaV0ToPrimVertex,rMomPos,rMomNeg,rDcaDaughterToPrimVertex);
  updated_v0->SetCharge(aodTrack1->Charge()+aodTrack2->Charge());
  
  if (0) {

    AliAODMCParticle* particle1 = (AliAODMCParticle *)fMCTrackArray->At(aodTrack1->GetLabel());
    AliAODMCParticle* particle2 = (AliAODMCParticle *)fMCTrackArray->At(aodTrack2->GetLabel());
    
    if(particle1 && particle2) {
      AliAODMCParticle* mother1 = (AliAODMCParticle *)fMCTrackArray->At(particle1->GetMother());
      AliAODMCParticle* mother2 = (AliAODMCParticle *)fMCTrackArray->At(particle2->GetMother());
      if (mother1 && particle1->GetMother()==particle2->GetMother() && mother1->GetPdgCode() == 310) {
	//if (mother->GetMother()<0) {
	  if (isSignalFlag) cout<<"=================================================================================================================="<<endl;
	  cout<<"MC Dau Mom1 : "<<particle1->Px()<<"   "<<particle1->Py()<<"   "<<particle1->Pz()<<endl;
	  cout<<"KF Dau Mom1 : "<<kfTrack[0].Px()<<"   "<<kfTrack[0].Py()<<"   "<<kfTrack[0].Pz()<<endl;
	  cout<<"AOD Dau Mom1: "<<aodTrack1->Px()<<"   "<<aodTrack1->Py()<<"   "<<aodTrack1->Pz()<<endl;
	  cout<<"KF Dau Pos1 : "<<kfTrack[0].GetX()<<"   "<<kfTrack[0].GetY()<<"   "<<kfTrack[0].GetZ()<<endl;
	  cout<<"AOD Dau Pos1: "<<aodTrack1->Xv()<<"   "<<aodTrack1->Yv()<<"   "<<aodTrack1->Zv()<<endl;
	  
	  cout<<"MC VTX     : "<<particle1->Xv()<<"   "<<particle1->Yv()<<"   "<<particle1->Zv()<<endl;
	  cout<<"KF VTX     : "<<kfK0s.GetX()<<"   "<<kfK0s.GetY()<<"   "<<kfK0s.GetZ()<<endl;
	  cout<<"MC MOM     : "<<mother1->Px()<<"   "<<mother1->Py()<<"   "<<mother1->Pz()<<endl;
	  cout<<"KF MOM     : "<<kfK0s.Px()<<"   "<<kfK0s.Py()<<"   "<<kfK0s.Pz()<<endl;
	  	  
	  cout<<endl;	  
	  //}
      }      
    }
  }
  


  return updated_v0;
}

void AliAnalysisTaskAODTrackPair::checkPrimaryTrack(){
  
  int nTrack = fEvent->GetNumberOfTracks();
  
  AliAODTrack* leading_track=NULL;
  
  double highest_pt = 0.;

  for (int iTrack = 0; iTrack < nTrack; ++iTrack){    
    AliAODTrack* track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));    
    
    if ( !isAcceptPrimaryTrack(track) )
      continue;

    if ( fUtils->isAcceptMidPid(track,AliPID::kPion) )
      fArrayPrimaryPionTrack->Add(track);     

    if ( track->Pt() < fMinLeadingTrackPt )
      continue; 

    if ( track->Pt() > highest_pt ) 
      leading_track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrack));          
  }
  
  fLeadingTrack = leading_track;
}

bool AliAnalysisTaskAODTrackPair::isKstarCandidate(AliAODv0* v0){
  
  unique_ptr<TLorentzVector> K0s(new TLorentzVector());
  K0s->SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), massK0s);  

  unique_ptr<TLorentzVector> pion(new TLorentzVector());

  AliAODTrack* daughter1 = (AliAODTrack*)v0->GetDaughter(0);
  AliAODTrack* daughter2 = (AliAODTrack*)v0->GetDaughter(1);

  int id1 = daughter1->GetID();
  int id2 = daughter2->GetID();
  
  for (int iTrack=0; iTrack<fArrayPrimaryPionTrack->GetEntriesFast(); ++iTrack) {
    
    AliAODTrack* track = static_cast<AliAODTrack*>(fArrayPrimaryPionTrack->At(iTrack));
    
    pion->SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), massPion);
    
    unique_ptr<TLorentzVector> Kstar(new TLorentzVector(K0s->Vect()+pion->Vect(),K0s->E()+pion->E()));
    
    int id = track->GetID();
    
    double fill[] = {Kstar->M(),Kstar->Pt()};
    
    if ( id!=id1 && id!=id2 ) fSparseK0sPion->Fill(fill);
    
    if ((0.85<Kstar->M() && Kstar->M()<0.95) && id!=id1 && id!=id2) {
      return true;
    }
  }
  
  return false;
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
  
  KFParticle kfp;
  kfp.Initialize();  
  kfp.Create(fkP,fkC,track->Charge(),massPion);

  return kfp;
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
  
  vector<int> plabels1;
  vector<int> plabels2;
  vector<int> plabels12;
  vector<int> v0label1;
  vector<int> v0label2;
  vector<int> v0label12;
  
  for (int iV0 = 0; iV0 < fNK0s; ++iV0) {

    AliAODv0* v0 = static_cast<AliAODv0*>(fArrayK0s->At(iV0));
    AliAODTrack* track1 = static_cast<AliAODTrack*>(v0->GetDaughter(0));
    AliAODTrack* track2 = static_cast<AliAODTrack*>(v0->GetDaughter(1));
    
    plabels12.push_back(track1->GetLabel());
    plabels12.push_back(track2->GetLabel());
    v0label12.push_back(iV0);
    v0label12.push_back(iV0);
    
    plabels1.push_back(track1->GetLabel());
    plabels2.push_back(track2->GetLabel());
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
    
    double fill_k0s[] = {particle1->Pt(), particle1->Eta(), particle1->Phi()};
    fSparseTrueK0s->Fill(fill_k0s);
    
    bool det = false;
    int idv0 = 0;
    
    for (int iV0 = 0; iV0 < fNK0s; ++iV0) {
      if ((plabels1[iV0]==label1 && plabels2[iV0]==label2) || (plabels1[iV0]==label2 && plabels2[iV0]==label1)) {
	det  = true;
	idv0 = iV0;
      }
    }
    
    if (det){      
      
      AliAODv0* v0 = static_cast<AliAODv0*>(fArrayK0s->At(idv0));
      
      double trkDist1 = v0->DcaV0Daughters();
      double cpa1     = v0->CosPointingAngle(fPrimVtxPos);
      
      if (trkDist1 < 1.0) {
	if (cpa1>0.997) {
	  fHistResoK0sMassPt_PCAcut997_TrkDist1 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	  fHistResoK0sPtPt_PCAcut997_TrkDist1   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	  fHistResoK0sEtaPt_PCAcut997_TrkDist1  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	  fHistResoK0sPhiPt_PCAcut997_TrkDist1  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	  fSparseRecK0s_PCAcut997_TrkDist1      -> Fill(fill_k0s);
	}
	if (cpa1>0.999) {
	  fHistResoK0sMassPt_PCAcut999_TrkDist1 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	  fHistResoK0sPtPt_PCAcut999_TrkDist1   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	  fHistResoK0sEtaPt_PCAcut999_TrkDist1  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	  fHistResoK0sPhiPt_PCAcut999_TrkDist1  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	  fSparseRecK0s_PCAcut999_TrkDist1      -> Fill(fill_k0s);
	}
	if (cpa1>0.97) {
	  fHistResoK0sMassPt_PCAcut97_TrkDist1 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	  fHistResoK0sPtPt_PCAcut97_TrkDist1   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	  fHistResoK0sEtaPt_PCAcut97_TrkDist1  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	  fHistResoK0sPhiPt_PCAcut97_TrkDist1  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	  fSparseRecK0s_PCAcut97_TrkDist1->Fill(fill_k0s);
	}
	if (cpa1>0.95) {
	  fHistResoK0sMassPt_PCAcut95_TrkDist1 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	  fHistResoK0sPtPt_PCAcut95_TrkDist1   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	  fHistResoK0sEtaPt_PCAcut95_TrkDist1  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	  fHistResoK0sPhiPt_PCAcut95_TrkDist1  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	  fSparseRecK0s_PCAcut95_TrkDist1->Fill(fill_k0s);
	}

	if (trkDist1 < 0.4) {
	  if (cpa1>0.997) {
	    fHistResoK0sMassPt_PCAcut997_TrkDist04 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	    fHistResoK0sPtPt_PCAcut997_TrkDist04   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	    fHistResoK0sEtaPt_PCAcut997_TrkDist04  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	    fHistResoK0sPhiPt_PCAcut997_TrkDist04  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	    fSparseRecK0s_PCAcut997_TrkDist04->Fill(fill_k0s);
	  }
	  if (cpa1>0.999) {
	    fHistResoK0sMassPt_PCAcut999_TrkDist04 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	    fHistResoK0sPtPt_PCAcut999_TrkDist04   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	    fHistResoK0sEtaPt_PCAcut999_TrkDist04  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	    fHistResoK0sPhiPt_PCAcut999_TrkDist04  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	    fSparseRecK0s_PCAcut999_TrkDist04->Fill(fill_k0s);
	  }
	  if (cpa1>0.97) {
	    fHistResoK0sMassPt_PCAcut97_TrkDist04 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	    fHistResoK0sPtPt_PCAcut97_TrkDist04   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	    fHistResoK0sEtaPt_PCAcut97_TrkDist04  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	    fHistResoK0sPhiPt_PCAcut97_TrkDist04  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	    fSparseRecK0s_PCAcut97_TrkDist04->Fill(fill_k0s);
	  }
	  if (cpa1>0.95) {
	    fHistResoK0sMassPt_PCAcut95_TrkDist04 -> Fill(particle1->Pt(),v0->MassK0Short()/particle1->M());
	    fHistResoK0sPtPt_PCAcut95_TrkDist04   -> Fill(particle1->Pt(),v0->Pt()/particle1->Pt());
	    fHistResoK0sEtaPt_PCAcut95_TrkDist04  -> Fill(particle1->Pt(),v0->Eta()/particle1->Eta());
	    fHistResoK0sPhiPt_PCAcut95_TrkDist04  -> Fill(particle1->Pt(),v0->Phi()/particle1->Phi());
	    fSparseRecK0s_PCAcut95_TrkDist04->Fill(fill_k0s);
	  }
	}

      }                       
    }    
  }
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sPair(AliAODv0* v0_1, AliAODv0* v0_2){
  
  AliAODTrack* aodTrack1[2] = {static_cast<AliAODTrack*>(v0_1->GetDaughter(0)),static_cast<AliAODTrack*>(v0_1->GetDaughter(1))};
  AliAODTrack* aodTrack2[2] = {static_cast<AliAODTrack*>(v0_2->GetDaughter(0)),static_cast<AliAODTrack*>(v0_2->GetDaughter(1))};

  int iTrack1[2];  
  iTrack1[0] = aodTrack1[0]->GetID();
  iTrack1[1] = aodTrack1[1]->GetID();

  int iTrack2[2];  
  iTrack2[0] = aodTrack2[0]->GetID();
  iTrack2[1] = aodTrack2[1]->GetID();
  
  return isAcceptK0sPair(iTrack1,iTrack2);
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sPair(int iTrack1[2], int iTrack2[2]){
  if (iTrack1[0] == iTrack2[0] || iTrack1[0] == iTrack2[1] || iTrack1[1] == iTrack2[0] || iTrack1[1] == iTrack2[1]) {
    return false;
  } else {
    return true;
  }
}

bool AliAnalysisTaskAODTrackPair::isAcceptV0QualityCuts(AliAODv0* v0, bool isAcceptCuts[7]){  
  
  double p  = v0->P();
  
  double chi2            = v0->Chi2V0();
  double cpv             = v0->CosPointingAngle(fPrimVtxPos);
  double lengthXY        = v0->DecayLengthXY(fPrimVtxPos);
  double trackDistanceXY = v0->DcaV0Daughters();
  double DCAxy           = v0->DcaV0ToPrimVertex();
  double lifetime        = fUtils->fPdgK0sMass * v0->DecayLength(fPrimVtxPos) / p;
  
  return isAcceptV0QualityCuts(p,chi2,cpv,lengthXY,trackDistanceXY,DCAxy,lifetime,isAcceptCuts);
}

bool AliAnalysisTaskAODTrackPair::isAcceptV0QualityCuts(double p, double chi2, double cpv, double lengthXY,
							double trackDistanceXY, double DCAxy, double lifetime, bool isAcceptCuts[7]){     
  isAcceptCuts[0] = false;
  isAcceptCuts[1] = false;
  isAcceptCuts[2] = false;
  isAcceptCuts[3] = false;
  isAcceptCuts[4] = false;
  isAcceptCuts[5] = false;
  isAcceptCuts[6] = false;

  double max_chi2            = 0;
  double min_cpv             = 0;
  double max_lengthXY        = 0;
  double max_trackDistanceXY = 0;
  double max_DCAxy           = 0;
  double max_lifetime        = 0;
  
  max_chi2 = 9999.;
  min_cpv = 0.90;
  max_lengthXY = 0;
  max_trackDistanceXY = 0;
  max_DCAxy = 0;
  max_lifetime = 20;
  
  return true;
  /*

  if (onChi2perNDFCut) {
    if (!fHistChi2perNDFCut) cout<<"Missing "<<endl;
    max_chi2            = fHistChi2perNDFCut->GetBinContent(fHistChi2perNDFCut->GetXaxis()->FindBin(p));
  }
  if (onPointingAngleCut){
    if (!fHistPointingAngleCut) cout<<"Missing fHistPointingAngleCut"<<endl;
    min_cpv             = fHistPointingAngleCut->GetBinContent(fHistPointingAngleCut->GetXaxis()->FindBin(p));
    
  }
  if (onDecayLengthCut){
    if (!fHistDecayLengthCut) cout<<"Missing fHistDecayLengthCut"<<endl;
    max_lengthXY        = fHistDecayLengthCut->GetBinContent(fHistDecayLengthCut->GetXaxis()->FindBin(p));
  }
  if (onDaughterTrackDistanceXYCut){
    if (!fHistDaughterTrackDistanceXYCut) cout<<"Missing fHistDaughterTrackDistanceXYCut"<<endl;
    max_trackDistanceXY = fHistDaughterTrackDistanceXYCut->GetBinContent(fHistDaughterTrackDistanceXYCut->GetXaxis()->FindBin(p));
    max_trackDistanceXY = 1.0;
  }
  if (onDCACut) {
    if (!fHistDCACut) cout<<"Missing fHistDCACut"<<endl;
    max_DCAxy           = fHistDCACut->GetBinContent(fHistDCACut->GetXaxis()->FindBin(p));
  }
  if (onLifetimeCut) {
    if (!fHistProperLifeTimeCut) cout<<"Missing fHistProperLifeTimeCut"<<endl;
    max_lifetime        = fHistProperLifeTimeCut->GetBinContent(fHistProperLifeTimeCut->GetXaxis()->FindBin(p));
  }

  if ( onChi2perNDFCut ) {
    if ( chi2 < max_chi2 ) isAcceptCuts[0] = true;
    else                   isAcceptCuts[0] = false;
  } else {
    isAcceptCuts[0] = true;
  }
  
  if ( onPointingAngleCut ) {
    if ( cpv > min_cpv ) isAcceptCuts[1] = true;
    else                 isAcceptCuts[1] = false;
  } else {
    isAcceptCuts[1] = true;
  }    
  
  if ( onDecayLengthCut ) {
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
  
  if ( onDCACut ) {
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
  
  if ( onArmenterosCut ) {
  
  //if ( isAcceptK0sArmenteros(v0) ) isAcceptCuts[6] = true;
//    else                             isAcceptCuts[6] = false;

    isAcceptCuts[6] = true;
  } else {
    isAcceptCuts[6] = true;
  }
  
  if (isAcceptCuts[0]==true && isAcceptCuts[1]==true && 
      isAcceptCuts[2]==true && isAcceptCuts[3]==true && 
      isAcceptCuts[4]==true && isAcceptCuts[5]==true && 
      isAcceptCuts[6]==true) return true;
  else return false;


  */

}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sKinematics(AliAODv0* v0) {  
  unique_ptr<TLorentzVector> lv(new TLorentzVector());
  lv->SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), v0->MassK0Short());  
  return isAcceptK0sKinematics(lv->Rapidity(), lv->Pt());
}
bool AliAnalysisTaskAODTrackPair::isAcceptK0sKinematics(AliAODMCParticle* v0) {
  unique_ptr<TLorentzVector> lv(new TLorentzVector());
  lv->SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), v0->M());  
  return isAcceptK0sKinematics(lv->Rapidity(), lv->Pt());
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sKinematics(double rap, double pt) {
  if (fMinK0sRap > rap || fMaxK0sRap < rap) return false;
  if (fMinK0sPt > pt || fMaxK0sPt < pt)     return false; 
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptPrimaryTrack(AliAODTrack* track){
  /*
  if ( (track->GetStatus() & AliVTrack::kTPCrefit) == 0 )               return false; 
  if ( (track->GetStatus() & AliVTrack::kITSrefit) == 0 )               return false; 
  if ( track->GetTPCNcls() < 70 )                                       return false;
  if ( track->GetITSNcls() < 1 )                                        return false;  
  if ( !track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) ) return false;
  if ( track->GetKinkIndex(0) != 0)                                     return false;
  if ( track->GetTPCNclsF() / track->GetTPCCrossedRows() < 0.8 )        return false;
  if ( track->GetTPCchi2() / track->GetTPCNcls() > 4 )                  return false;
  if ( track->GetITSchi2() / track->GetITSNcls() > 36 )                 return false;
  
  float dca_xy = 9999;
  float dca_z  = 9999;  
  track->GetImpactParameters(dca_xy, dca_z);
  
  double cut_dca_xy = 0.0105 + 0.035/pow(track->Pt(),1.1);
  
  if ( fabs(dca_xy) > cut_dca_xy ) return false;
  if ( fabs(dca_z) > 2.0 )         return false;
  */
  if (!track->TestFilterMask(AliAODTrack::kTrkGlobal)) return false;
  
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptK0sCandidateMassRange(double mass){
  if ( fMinK0sMassRange>mass || mass > fMaxK0sMassRange) return false;
  else                                                   return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptDaughterTrackPairAngle(AliAODTrack* track1, AliAODTrack* track2){  
  TVector3 v1(track1->Px(),track1->Py(),track1->Pz());
  TVector3 v2(track2->Px(),track2->Py(),track2->Pz());
  if (v1.Angle(v2)<1E-4)  return false;
  else                    return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptDaughterTrack(AliAODTrack *track, AliPID::EParticleType pid1) {  

  if ( !isAcceptDaughterTrackKinematics(track) ) return false;
  if ( !isAcceptDaughterTrackQuality(track) )    return false;
  if ( !fUtils->isAcceptMidPid(track,pid1) )     return false;  
  if ( !fIsMC || (fIsMC && fUtils->getMotherPdgCode(track) == 310 && fUtils->getGrandMotherLabel(track->GetLabel()) == -1) ) {
    TrackQualityChecker(track);    
    TrackPIDChecker(track,pid1,true);
  }
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptDaughterTrackKinematics(AliAODTrack *track) {
  if (track->P() < fMinTrackP     || fMaxTrackP < track->P())     return false;
  if (track->Eta() < fMinTrackEta || fMaxTrackEta < track->Eta()) return false;
  if (track->Pt() < fMinTrackPt   || fMaxTrackPt < track->Pt())   return false;
  return true;
}

bool AliAnalysisTaskAODTrackPair::isAcceptDaughterTrackQuality(AliAODTrack *track) {
  
  if ( !track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) && 
       !track->HasPointOnITSLayer(4) && !track->HasPointOnITSLayer(5) && 
       track->GetTOFBunchCrossing() != 0 )                                                                   return false;
  if ( track->GetTOFBunchCrossing() == 0 && track->GetTOFsignal() > 99998)                                   return false; 
  if ( (track->GetStatus() & AliVTrack::kTPCrefit) == 0 )                                                    return false;
  if ( fMinTrackTPCNClusts > track->GetTPCNcls() )                                                           return false;  
  if ( track->GetKinkIndex(0) != 0 )                                                                         return false;
  if ( fMinCrossRowsFindableRatio > track->GetTPCNclsF() / track->GetTPCCrossedRows() )                      return false;
  if ( fMaxReducedChi2TPC < track->GetTPCchi2() / track->GetTPCNcls() )                                      return false;
  
  double its_chi2 = track->GetITSNcls() > 0 ? track->GetITSchi2() / track->GetITSNcls() : 0;
  if ( fMaxReducedChi2ITS < its_chi2 )                                                                       return false;

  double maxd   = 999.;
  double dz[2]  = {};
  double c[3]   = {};

  if (track->PropagateToDCA(fPrimVtx,fEvent->GetMagneticField(),maxd,dz,c) ){ 
    if ( fabs(dz[1]) > 10 )   return false;
  } 
  

  /*
  KFParticle kfTrack;
  kfTrack = CreateKFParticle(track, 211);
  float dca, e_dca;
  float fPrimVtxPosF[] = {static_cast<float>(fPrimVtxPos[0]),static_cast<float>(fPrimVtxPos[1]),static_cast<float>(fPrimVtxPos[2])};
  float fPrimVtxCovF[] = {static_cast<float>(fPrimVtxCov[0]),static_cast<float>(fPrimVtxCov[1]),static_cast<float>(fPrimVtxCov[2]),
			  static_cast<float>(fPrimVtxCov[3]),static_cast<float>(fPrimVtxCov[4]),static_cast<float>(fPrimVtxCov[5])};
  kfTrack.GetDistanceFromVertexXY(fPrimVtxPosF,fPrimVtxCovF,dca,e_dca);  
  //if ( fabs(dca) < 0.05 ) return false;
  */
  return true;
}

