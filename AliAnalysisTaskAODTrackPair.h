#ifndef AliAnalysisTaskAODTrackPair_cxx
#define AliAnalysisTaskAODTrackPair_cxx

#include "AliAnalysisTaskAODTrackPairUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "TH3F.h"
#include "THnSparse.h"
//#include "TObjectArray.h"

#include "KFParticleBase.h"
#include "KFParticle.h"

#include <time.h>

class AliAnalysisTaskAODTrackPair : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskAODTrackPair();
  AliAnalysisTaskAODTrackPair(const char *name);
  virtual ~AliAnalysisTaskAODTrackPair();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void setMC(bool isMC) { fIsMC = isMC; }
  void setManualV0Analysis(bool isManual) {fIsManualV0Analysis=isManual;}
  
  void setMixingAnalysis(bool isMix) { fIsMixingAnalysis = isMix; }
  void setCentralityMethod(std::string method){fMethodCent = method;}


  void setUtils(AliAnalysisTaskAODTrackPairUtils *utils) { fUtils = utils; }
  
  void setEvtMixingTrackDepth(int depth) { fTrackDepth = depth; }
  void setEvtMixingPoolSize(int size) { fPoolSize = size; }
  void setEvtMixingReadyFraction(double frac) { fReadyFraction = frac; }
  void setEvtMixingPoolVtxZ(bool flag) { onEvtMixingPoolVtxZ = flag; }
  void setEvtMixingPoolCent(bool flag) { onEvtMixingPoolCent = flag; }
  void setEvtMixingPoolPsi(bool flag) { onEvtMixingPoolPsi = flag; }
  
  void setVtxZRange(double min, double max) {
    fMinVertexCutZ = min;
    fMaxVertexCutZ = max;
  }

  void setK0sPtRange(double min, double max){
    fMinK0sPt = min;
    fMaxK0sPt = max;
  }
  void setK0sRapRange(double min, double max){
    fMinK0sRap = min;
    fMaxK0sRap = max;
  }
  void setK0sDaughterTrackPRange(double min, double max){
    fMinTrackP = min;
    fMaxTrackP = max;
  }
  void setK0sDaughterTrackPtRange(double min, double max){
    fMinTrackPt = min;
    fMaxTrackPt = max;
  }
  void setK0sDaughterTrackEtaRange(double min, double max){
    fMinTrackEta = min;
    fMaxTrackEta = max;
  }
  void setK0sDaughterTrackChi2TPC(double max){
    fMaxReducedChi2TPC = max;
  }
  void setK0sDaughterTrackChi2ITS(double max){
    fMaxReducedChi2ITS = max;
  }
  void setK0sDaughterTrackNClustTPC(double min){
    fMinTrackTPCNClusts = min;
  }
  void setK0sDaughterTrackNClustSPD(double min){    
    fMinTrackSPDNClusts = min;
  }
  void setK0sDaughterTrackFindableTPC(double min){
    fMinCrossRowsFindableRatio = min;
  }

  void setK0sCuts(TFile *inFile,bool decaylength, bool pointangle, bool chi2, bool paiDCA, bool trackdistance, bool dcaxy, bool lifetime, bool armenteros){
    onDecayLengthCut             = decaylength;
    onPointingAngleCut           = pointangle;
    onChi2perNDFCut              = chi2;
    onDaughterPairDCAXYCut       = paiDCA;
    onDaughterTrackDistanceXYCut = trackdistance;
    onDCACut                     = dcaxy;
    onLifetimeCut                = lifetime;
    onArmenterosCut              = armenteros;
    
    if (onDecayLengthCut){
      fHistDecayLengthCut             = (TH1F*)inFile -> Get("fHistDecayLengthCut");
      if (!fHistDecayLengthCut) cout<<" Missing set fHistDecayLengthCut"<<endl;
      else                      cout<<" Find set fHistDecayLengthCut"<<endl;
    }
    if (onPointingAngleCut){
      fHistPointingAngleCut           = (TH1F*)inFile -> Get("fHistPointingAngleCut");
      if (!fHistPointingAngleCut) cout<<" Missing set fHistPointingAngleCut"<<endl;      
    }
    if (onChi2perNDFCut){
      fHistChi2perNDFCut              = (TH1F*)inFile -> Get("fHistChi2perNDFCut");
      if (!fHistChi2perNDFCut) cout<<" Missing set fHistChi2perNDFCut"<<endl;
    }
    if (onDaughterPairDCAXYCut){
      fHistDaughterPairDCAXYCut       = (TH1F*)inFile -> Get("fHistDaughterPairDCAXYCut");
      if (!fHistDaughterPairDCAXYCut) cout<<" Missing set fHistDaughterPairDCAXYCut"<<endl;
    }
    if (onDaughterTrackDistanceXYCut){
      fHistDaughterTrackDistanceXYCut = (TH1F*)inFile -> Get("fHistTrackDistanceXYCut");
      if (!fHistDaughterTrackDistanceXYCut) cout<<" Missing set fHistDaughterTrackDistanceXYCut"<<endl;
    }
    if (onDCACut){
      fHistDCACut                     = (TH1F*)inFile -> Get("fHistDCACut");
      if (!fHistDCACut) cout<<" Missing set fHistDCACut"<<endl;
    }
    if (onLifetimeCut){
      fHistProperLifeTimeCut          = (TH1F*)inFile -> Get("fHistPropLifeTimeCut");
      if (!fHistProperLifeTimeCut) cout<<" Missing set fHistProperLifeTimeCut"<<endl;
    }
  }

private:
  AliAnalysisTaskAODTrackPair(
      const AliAnalysisTaskAODTrackPair &); // not implemented
  AliAnalysisTaskAODTrackPair &
  operator=(const AliAnalysisTaskAODTrackPair &); // not implemented

  bool Initialize();
      
  bool TrackPIDChecker(AliAODTrack *track, AliPID::EParticleType pid, bool isSel);  
  bool TrackQualityChecker(AliAODTrack *track);

  bool V0QualityChecker(AliAODv0* v0, bool isSel);
  bool V0QualityChecker(double mass, double p, int charge, double DCApairXY, double DCApair, double decay_length, double decay_length_xy, double cpv,
			double cpv_xy, double dca, double dca_xy, double distance_daughter, double distance_daughter_xy,
			double chi2, double lifetime, double armenteros_alpha, double armenteros_qt, int dlabel[2], bool isSel);

  bool K0sPairAnalysis(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  bool K0sPairAnalysisEventMixing(AliPID::EParticleType pid1,AliPID::EParticleType pid2);

  bool AddK0sArray(AliPID::EParticleType pid1, AliPID::EParticleType pid2);

  AliAODv0 *updateAODv0(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2);

  void checkPrimaryTrack();
  
  bool isKstarCandidate(AliAODv0* v0);

  KFParticle CreateKFParticle(AliAODTrack* track,int pdg);
  
  double DecayLengthFromKF(KFParticle kfpParticle, KFParticle PV);
  double DecayLengthXYFromKF(KFParticle kfpParticle, KFParticle PV);

  bool ProcessMC();

  bool isAcceptK0sPair(AliAODv0* v0_1, AliAODv0* v0_2);
  bool isAcceptK0sPair(int label1[2], int label2[2]);

  bool isAcceptV0QualityCuts(AliAODv0* v0, bool isAcceptCuts[7]);
  bool isAcceptV0QualityCuts(double p, double chi2, double cpv, double lengthXY, double trackDistanceXY, double DCAxy, double lifetime, bool isAcceptCuts[7]);
  
  bool isAcceptK0sKinematics(AliAODv0* v0);
  bool isAcceptK0sKinematics(AliAODMCParticle* v0);
  bool isAcceptK0sKinematics(double rap, double pt);
  
  bool isAcceptPrimaryTrack(AliAODTrack* track);
  bool isAcceptK0sCandidateMassRange(double mass);
  
  bool isAcceptDaughterTrackPairAngle(AliAODTrack* track1, AliAODTrack* track2);
  bool isAcceptDaughterTrack(AliAODTrack *track, AliPID::EParticleType pid1);
  bool isAcceptDaughterTrackKinematics(AliAODTrack *track);
  bool isAcceptDaughterTrackQuality(AliAODTrack *track);
  

  AliAODEvent *fEvent;
  AliAODVertex *fPrimVtx;
  AliAODTrack *fLeadingTrack;
 
  AliEventPoolManager *fPoolMuonTrackMgrK0s;  
  AliEventPoolManager *fPoolMuonTrackMgrPion;  
  
  AliMultSelection *fMultSelection;
  AliPIDResponse *fPIDResponse;

  AliAnalysisTaskAODTrackPairUtils *fUtils;
  
  TClonesArray *fMCTrackArray;

  TObjArray *fArrayK0s;
  TObjArray *fArrayK0sDaughterTrack;
  TObjArray *fArrayPrimaryPionTrack;
  
  TH1F* fHistDecayLengthCut;
  TH1F* fHistPointingAngleCut;
  TH1F* fHistChi2perNDFCut;
  TH1F* fHistDaughterPairDCAXYCut;
  TH1F* fHistDaughterTrackDistanceXYCut;
  TH1F* fHistDCACut;
  TH1F* fHistProperLifeTimeCut;
  
  int fTrackDepth;
  int fPoolSize;

  double fReadyFraction;
  
  bool onVerbose;

  bool onEvtMixingPoolVtxZ;
  bool onEvtMixingPoolCent;
  bool onEvtMixingPoolPsi;

  bool onDecayLengthCut;
  bool onPointingAngleCut;
  bool onChi2perNDFCut;
  bool onDaughterPairDCAXYCut;
  bool onDaughterTrackDistanceXYCut;
  bool onDCACut;
  bool onLifetimeCut;
  bool onArmenterosCut;

  bool fIsMC;
  bool fIsMixingAnalysis;
  bool fIsManualV0Analysis;
  
  double fArmenterosAlpha;
  
  double fMinVertexCutZ;
  double fMaxVertexCutZ;

  double fMinK0sRap;
  double fMaxK0sRap;
  double fMinK0sPt;
  double fMaxK0sPt;
  double fMinK0sMassRange;
  double fMaxK0sMassRange;
  
  int fMinTrackTPCNClusts;
  int fMinTrackSPDNClusts;
  double fMaxReducedChi2TPC;
  double fMaxReducedChi2ITS;
  double fMinCrossRowsFindableRatio;
  
  double fMinTrackP;
  double fMaxTrackP;
  double fMinTrackPt;
  double fMaxTrackPt;
  double fMinTrackEta;
  double fMaxTrackEta;
  
  double fMinLeadingTrackPt;
  
  std::string fMethodCent;
  
  double fPrimVtxPos[3];
  double fPrimVtxCov[6];
  double fCent;
  double fPsi;
  
  int fNK0s;

  int fMinNContPrimVtx;

  const double massPion = 0.13957000;
  const double massK0s  = 0.49761400;
  
  ////////////////////////////////////////////////
  // Output histos
  ////////////////////////////////////////////////

  TList *fOutputList;
  
  TH1F *fHistEventVtxZ;
  TH1F *fHistEventCent;
  TH1F *fHistEventMulti;
  TH1F *fHistEventVtxCont;

  THnSparse *fSparseULSPionPairBeforeCuts;

  THnSparse *fSparseULSPionPair_PCAcut997_TrkDist1;
  THnSparse *fSparseULSPionPair_PCAcut999_TrkDist1;
  THnSparse *fSparseULSPionPair_PCAcut97_TrkDist1;
  THnSparse *fSparseULSPionPair_PCAcut95_TrkDist1;
  THnSparse *fSparseULSPionPair_PCAcut997_TrkDist04;
  THnSparse *fSparseULSPionPair_PCAcut999_TrkDist04;
  THnSparse *fSparseULSPionPair_PCAcut97_TrkDist04;
  THnSparse *fSparseULSPionPair_PCAcut95_TrkDist04;
  THnSparse *fSparseULSPionPair_KstarSel;
  

  THnSparse *fSparseNeutralK0sPair_PCAcut997_TrkDist1;
  THnSparse *fSparseNeutralK0sPair_PCAcut999_TrkDist1;
  THnSparse *fSparseNeutralK0sPair_PCAcut97_TrkDist1;
  THnSparse *fSparseNeutralK0sPair_PCAcut95_TrkDist1;
  THnSparse *fSparseNeutralK0sPair_PCAcut997_TrkDist04;
  THnSparse *fSparseNeutralK0sPair_PCAcut999_TrkDist04;
  THnSparse *fSparseNeutralK0sPair_PCAcut97_TrkDist04;
  THnSparse *fSparseNeutralK0sPair_PCAcut95_TrkDist04;
  THnSparse *fSparseNeutralK0sPair_KstarSel;

  THnSparse *fSparseMixNeutralK0sPair_PCAcut997_TrkDist1;
  THnSparse *fSparseMixNeutralK0sPair_PCAcut999_TrkDist1;
  THnSparse *fSparseMixNeutralK0sPair_PCAcut97_TrkDist1;  
  THnSparse *fSparseMixNeutralK0sPair_PCAcut95_TrkDist1;  
  THnSparse *fSparseMixNeutralK0sPair_PCAcut997_TrkDist04;
  THnSparse *fSparseMixNeutralK0sPair_PCAcut999_TrkDist04;
  THnSparse *fSparseMixNeutralK0sPair_PCAcut97_TrkDist04;
  THnSparse *fSparseMixNeutralK0sPair_PCAcut95_TrkDist04;
  THnSparse *fSparseMixNeutralK0sPair_KstarSel;
  
  THnSparse *fSparseNeutralK0sPairRejectKstar;
  THnSparse *fSparseMixNeutralK0sPairRejectKstar;

  THnSparse *fSparseNeutralK0sPairTransverse;
  THnSparse *fSparseNeutralK0sPairToward;
  THnSparse *fSparseNeutralK0sPairAway;
  THnSparse *fSparseNeutralK0sPairInJet;
  
  THnSparse *fSparseK0sPion;
  THnSparse *fSparseULSPionPairRejectKstar;
  
  THnSparse *fSparseRecK0s_PCAcut95_TrkDist1;
  THnSparse *fSparseRecK0s_PCAcut97_TrkDist1;
  THnSparse *fSparseRecK0s_PCAcut997_TrkDist1;
  THnSparse *fSparseRecK0s_PCAcut999_TrkDist1;
  THnSparse *fSparseRecK0s_PCAcut95_TrkDist04;
  THnSparse *fSparseRecK0s_PCAcut97_TrkDist04;
  THnSparse *fSparseRecK0s_PCAcut997_TrkDist04;
  THnSparse *fSparseRecK0s_PCAcut999_TrkDist04;  
  THnSparse *fSparseRecK0s_KstarSel;

  THnSparse *fSparseTrueK0s;

  TH2F* fHistResoK0sMassPt_PCAcut95_TrkDist04;
  TH2F* fHistResoK0sPtPt_PCAcut95_TrkDist04;
  TH2F* fHistResoK0sEtaPt_PCAcut95_TrkDist04;
  TH2F* fHistResoK0sPhiPt_PCAcut95_TrkDist04;

  TH2F* fHistResoK0sMassPt_PCAcut97_TrkDist04;
  TH2F* fHistResoK0sPtPt_PCAcut97_TrkDist04;
  TH2F* fHistResoK0sEtaPt_PCAcut97_TrkDist04;
  TH2F* fHistResoK0sPhiPt_PCAcut97_TrkDist04;

  TH2F* fHistResoK0sMassPt_PCAcut997_TrkDist04;
  TH2F* fHistResoK0sPtPt_PCAcut997_TrkDist04;
  TH2F* fHistResoK0sEtaPt_PCAcut997_TrkDist04;
  TH2F* fHistResoK0sPhiPt_PCAcut997_TrkDist04;

  TH2F* fHistResoK0sMassPt_PCAcut999_TrkDist04;
  TH2F* fHistResoK0sPtPt_PCAcut999_TrkDist04;
  TH2F* fHistResoK0sEtaPt_PCAcut999_TrkDist04;
  TH2F* fHistResoK0sPhiPt_PCAcut999_TrkDist04;

  TH2F* fHistResoK0sMassPt_PCAcut95_TrkDist1;
  TH2F* fHistResoK0sPtPt_PCAcut95_TrkDist1;
  TH2F* fHistResoK0sEtaPt_PCAcut95_TrkDist1;
  TH2F* fHistResoK0sPhiPt_PCAcut95_TrkDist1;

  TH2F* fHistResoK0sMassPt_PCAcut97_TrkDist1;
  TH2F* fHistResoK0sPtPt_PCAcut97_TrkDist1;
  TH2F* fHistResoK0sEtaPt_PCAcut97_TrkDist1;
  TH2F* fHistResoK0sPhiPt_PCAcut97_TrkDist1;

  TH2F* fHistResoK0sMassPt_PCAcut997_TrkDist1;
  TH2F* fHistResoK0sPtPt_PCAcut997_TrkDist1;
  TH2F* fHistResoK0sEtaPt_PCAcut997_TrkDist1;
  TH2F* fHistResoK0sPhiPt_PCAcut997_TrkDist1;

  TH2F* fHistResoK0sMassPt_PCAcut999_TrkDist1;
  TH2F* fHistResoK0sPtPt_PCAcut999_TrkDist1;
  TH2F* fHistResoK0sEtaPt_PCAcut999_TrkDist1;
  TH2F* fHistResoK0sPhiPt_PCAcut999_TrkDist1;

  TH2F* fHistResoK0sMassPt_KstarSel;
  TH2F* fHistResoK0sPtPt_KstarSel;
  TH2F* fHistResoK0sEtaPt_KstarSel;
  TH2F* fHistResoK0sPhiPt_KstarSel;

  TH1F* fHistNumK0sCandidates;

  TH2F *fHistTPCdEdxP;
  TH2F *fHistTPCSigmaElectron;
  TH2F *fHistTPCSigmaMuon;
  TH2F *fHistTPCSigmaPion;
  TH2F *fHistTPCSigmaKaon;
  TH2F *fHistTPCSigmaProton;  
  TH2F *fHistBetaP;  
  TH2F *fHistTOFSigmaElectron;  
  TH2F *fHistTOFSigmaMuon;  
  TH2F *fHistTOFSigmaPion;
  TH2F *fHistTOFSigmaKaon;
  TH2F *fHistTOFSigmaProton;  
  
  TH2F *fHistSelTPCdEdxP;
  TH2F *fHistSelTPCSigmaElectron;
  TH2F *fHistSelTPCSigmaMuon;
  TH2F *fHistSelTPCSigmaPion;
  TH2F *fHistSelTPCSigmaKaon;
  TH2F *fHistSelTPCSigmaProton;  
  TH2F *fHistSelBetaP;  
  TH2F *fHistSelTOFSigmaElectron;  
  TH2F *fHistSelTOFSigmaMuon;  
  TH2F *fHistSelTOFSigmaPion;
  TH2F *fHistSelTOFSigmaKaon;
  TH2F *fHistSelTOFSigmaProton;
  
  TH2F *fHistTPCNClusts;
  TH2F *fHistSPDNClusts;
  TH2F *fHistTPCCrossRowsFindableRatio;
  TH2F *fHistReducedChi2TPC;
  TH2F *fHistReducedChi2ITS;
  TH2F *fHistDCAz;
  TH2F *fHistDCAxy;
  TH1F *fHistTrackP;
  TH1F *fHistTrackPt;
  TH1F *fHistTrackEta;    
  
  TH2F *fHistV0PV0DecayLengthXY;
  TH2F *fHistV0PV0DecayLength;
  TH2F *fHistV0PV0PointingAngleXY;
  TH2F *fHistV0PV0PointingAngle;
  TH2F *fHistV0PV0DCAXY;
  TH2F *fHistV0PV0DCA;
  TH2F *fHistV0PV0TrackDistanceXY;
  TH2F *fHistV0PV0TrackDistance;
  TH2F *fHistV0PV0Chi2perNDF;
  TH2F *fHistV0PV0PropLifeTime;
  TH2F *fHistV0TrackDCAXY; 
  TH2F *fHistV0TrackDCA; 
  TH2F *fHistV0Armenteros;

  TH2F *fHistSelV0PV0DecayLengthXY;
  TH2F *fHistSelV0PV0DecayLength;
  TH2F *fHistSelV0PV0PointingAngleXY;
  TH2F *fHistSelV0PV0PointingAngle;
  TH2F *fHistSelV0PV0DCAXY;
  TH2F *fHistSelV0PV0DCA;
  TH2F *fHistSelV0PV0TrackDistanceXY;
  TH2F *fHistSelV0PV0TrackDistance;
  TH2F *fHistSelV0PV0Chi2perNDF;
  TH2F *fHistSelV0PV0PropLifeTime;
  TH2F *fHistSelV0TrackDCAXY; 
  TH2F *fHistSelV0TrackDCA; 
  TH2F *fHistSelV0Armenteros;

  TH2F *fHistV0PV0DecayLengthXY_TrueK0s;
  TH2F *fHistV0PV0DecayLength_TrueK0s;
  TH2F *fHistV0PV0PointingAngleXY_TrueK0s;
  TH2F *fHistV0PV0PointingAngle_TrueK0s;
  TH2F *fHistV0PV0DCAXY_TrueK0s;
  TH2F *fHistV0PV0DCA_TrueK0s;
  TH2F *fHistV0PV0TrackDistanceXY_TrueK0s;
  TH2F *fHistV0PV0TrackDistance_TrueK0s;
  TH2F *fHistV0PV0Chi2perNDF_TrueK0s;
  TH2F *fHistV0PV0PropLifeTime_TrueK0s;
  TH2F *fHistV0TrackDCAXY_TrueK0s; 
  TH2F *fHistV0TrackDCA_TrueK0s; 
  TH2F *fHistV0Armenteros_TrueK0s;

  TH2F *fHistSelV0PV0DecayLengthXY_TrueK0s;
  TH2F *fHistSelV0PV0DecayLength_TrueK0s;
  TH2F *fHistSelV0PV0PointingAngleXY_TrueK0s;
  TH2F *fHistSelV0PV0PointingAngle_TrueK0s;
  TH2F *fHistSelV0PV0DCAXY_TrueK0s;
  TH2F *fHistSelV0PV0DCA_TrueK0s;
  TH2F *fHistSelV0PV0TrackDistanceXY_TrueK0s;
  TH2F *fHistSelV0PV0TrackDistance_TrueK0s;
  TH2F *fHistSelV0PV0Chi2perNDF_TrueK0s;
  TH2F *fHistSelV0PV0PropLifeTime_TrueK0s;
  TH2F *fHistSelV0TrackDCAXY_TrueK0s; 
  TH2F *fHistSelV0TrackDCA_TrueK0s; 
  TH2F *fHistSelV0Armenteros_TrueK0s;

  int nEvtAnalyzed;

  clock_t start;
  clock_t end;
  
  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
