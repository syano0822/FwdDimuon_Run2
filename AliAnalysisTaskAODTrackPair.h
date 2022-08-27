#ifndef AliAnalysisTaskAODTrackPair_cxx
#define AliAnalysisTaskAODTrackPair_cxx

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "TH3F.h"
#include "AliAnalysisTaskAODTrackPairUtils.h"

class AliAnalysisTaskAODTrackPair : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskAODTrackPair();
  AliAnalysisTaskAODTrackPair(const char* name);
  virtual ~AliAnalysisTaskAODTrackPair();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  void setMC(bool isMC){
    fIsMC = isMC;
  }
  
  void setMidTrackAna(bool isMidTrack){
    fIsMidTrackAna = isMidTrack;
  }  
  void setK0sAna(bool isK0s) {
    fIsV0TrackPairAna = isK0s;
  }
  void setKaonTrackAna(bool isKaon) {
    fIsPrimTrackPairAna = isKaon;
  }
  
  void setMixingAnalysis(bool isMix)
  {
    fIsMixingAnalysis = isMix;
  }
  void setUtils(AliAnalysisTaskAODTrackPairUtils *utils)
  {
    fUtils = utils;
  }
  void setEvtMixingTrackDepth(int depth)
  {
    fTrackDepth=depth;
  }
  void setEvtMixingPoolSize(int size)
  {
    fPoolSize=size;
  }
  void setEvtMixingReadyFraction(float frac)
  {
    fReadyFraction=frac;
  }
  void setEvtMixingPoolVtxZ(bool flag)
  {
    onEvtMixingPoolVtxZ=flag;
  }
  void setEvtMixingPoolCent(bool flag)
  {
    onEvtMixingPoolCent=flag;
  }
  void setEvtMixingPoolPsi(bool flag)
  {
    onEvtMixingPoolPsi=flag;
  }
  void setMixingEventTrigger(unsigned int mask)
  {
    fTriggerMaskForMixing = mask;
  }
  void setSameEventTrigger(unsigned int mask)
  {
    fTriggerMaskForSame = mask;
  }

 private:

  AliAnalysisTaskAODTrackPair(const AliAnalysisTaskAODTrackPair&); // not implemented
  AliAnalysisTaskAODTrackPair& operator=(const AliAnalysisTaskAODTrackPair&); // not implemented

  bool EventQA();
  bool Initialize();
  bool FwdMuonPairAnalysis();
  bool FwdMuonPairAnalysisEveMixing();
  bool FwdMuonTrackQA(AliAODTrack* track);
  bool FwdMuonPairQA(AliAODDimuon* dimuon);
  
  bool MidTrackQA(AliAODTrack* track);
  bool MidTrackPIDChecker(AliAODTrack* track, AliPID::EParticleType pid, bool isSel);
  bool MidTrackQualityChecker(AliAODTrack* track);
  bool MidV0Checker(AliAODv0* v0, bool isSel);

  bool MidMuonPairQA(AliAODDimuon* dimuon);
  bool MidPairAnalysis(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  
  bool MidV0Analysis(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  bool MidV0AnalysisEventMixing(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  
  bool MidMuonPairAnalysisEveMixing();
  
  AliAODEvent    *fEvent;
  AliEventPoolManager *fPoolMuonTrackMgr;
  AliAnalysisTaskAODTrackPairUtils* fUtils;

  bool fIsMC;
  bool fIsMidTrackAna;
  bool fIsV0TrackPairAna;
  bool fIsPrimTrackPairAna;
  bool fIsMixingAnalysis;

  int fRunNumber;
  int fTrackDepth;
  int fPoolSize;

  unsigned int fTriggerMaskForSame;
  unsigned int fTriggerMaskForMixing;

  float fReadyFraction;

  bool onEvtMixingPoolVtxZ;
  bool onEvtMixingPoolCent;
  bool onEvtMixingPoolPsi;

  bool fIsCINT7;
  bool fIsCMSL7;
  bool fIsCMSH7;
  bool fIsCMUL7;
  bool fIsCMLL7;

  ////////////////////////////////////////////////
  //Output histos
  ////////////////////////////////////////////////

  TList* fOutputList;
  TH2F* fEventCounter;

  TH2F* fHistTrackPairPtBalance;
  TH2F* fHistTrackPairLocalBoardPair;


  TH2F* fHistTrackThetaAbs;
  TH2F* fHistTrackTriggerMatch;
  TH2F* fHistTrackPDCA;
  TH2F* fHistTrackChiSquare;
  TH2F* fHistTriggerChiSquare;

  TH1F* fHistEventVtxZ;
  TH1F* fHistEventCent;
  TH1F* fHistEventMulti;
  TH1F* fHistEventVtxCont;

  TTree* fTreeULSPair;
  TTree* fTreeLSppPair;
  TTree* fTreeLSmmPair;
  
  TTree* fTreeULSPair_ProngV0;
  TTree* fTreeLSppPair_ProngV0;
  TTree* fTreeLSmmPair_ProngV0;

  TTree* fTreeULSPair_TightCut;
  TTree* fTreeLSppPair_TightCut;
  TTree* fTreeLSmmPair_TightCut;
  
  TTree* fTreeMixULSPair;
  TTree* fTreeMixLSppPair;
  TTree* fTreeMixLSmmPair;

  TH2F* fHistULSPairMassPt;
  TH2F* fHistLSppPairMassPt;
  TH2F* fHistLSmmPairMassPt;
  
  TH2F* fHistULSPairMassPt_ProngV0;
  TH2F* fHistLSppPairMassPt_ProngV0;
  TH2F* fHistLSmmPairMassPt_ProngV0;

  TH2F* fHistULSPairMassPt_TightCut;
  TH2F* fHistLSppPairMassPt_TightCut;
  TH2F* fHistLSmmPairMassPt_TightCut;
  
  TH2F* fHistMixULSPairMassPt;
  TH2F* fHistMixLSppPairMassPt;
  TH2F* fHistMixLSmmPairMassPt;

  TH2F* fHistMassK0s1K0s2;
  
  float RecPairPt;
  float RecPairRap;
  float RecPairMass;
  float RecPairArmenterosArmPt;
  float RecPairArmenterosAlpha;
  float RecPairCent;
  float RecPairDS;
  
  TH2F* fHistTPCdEdxP;
  TH2F* fHistBetaP;    
  TH2F* fHistTPCSigmaElectron;
  TH2F* fHistTOFSigmaElectron;
  TH2F* fHistTPCSigmaMuon;
  TH2F* fHistTOFSigmaMuon;
  TH2F* fHistTPCSigmaPion;
  TH2F* fHistTOFSigmaPion;
  TH2F* fHistTPCSigmaKaon;
  TH2F* fHistTOFSigmaKaon;
  TH2F* fHistTPCSigmaProton;
  TH2F* fHistTOFSigmaProton;

  TH2F* fHistSelTPCdEdxP;
  TH2F* fHistSelBetaP;    
  TH2F* fHistSelTPCSigmaElectron;
  TH2F* fHistSelTOFSigmaElectron;
  TH2F* fHistSelTPCSigmaMuon;
  TH2F* fHistSelTOFSigmaMuon;
  TH2F* fHistSelTPCSigmaPion;
  TH2F* fHistSelTOFSigmaPion;
  TH2F* fHistSelTPCSigmaKaon;
  TH2F* fHistSelTOFSigmaKaon;
  TH2F* fHistSelTPCSigmaProton;
  TH2F* fHistSelTOFSigmaProton;
  
  TH1F* fHistTrackP;
  TH1F* fHistTrackPt;
  TH1F* fHistTrackEta;
  TH1F* fHistTPCNClusts;
  TH1F* fHistSPDNClusts;
  TH1F* fHistTPCCrossRowsFindableRatio;
  TH1F* fHistReducedChi2TPC;
  TH1F* fHistReducedChi2ITS;
  TH1F* fHistDCAz;
  TH2F* fHistDCAxyPt;

  TH2F* fHistArmenteros;
  TH2F* fHistSelArmenteros; 
  TH2F* fHistV0MassDecayLength;
  TH2F* fHistV0MassPointingAngle;
  TH2F* fHistV0MassV0DCA;
  TH2F* fHistV0MassV0TrackDCA;
  TH2F* fHistV0MassV0DecayRadius;
  TH2F* fHistV0MassV0PropLifeTime;
  TH2F* fHistSelV0MassDecayLength;
  TH2F* fHistSelV0MassPointingAngle;
  TH2F* fHistSelV0MassV0DCA;
  TH2F* fHistSelV0MassV0TrackDCA;
  TH2F* fHistSelV0MassV0DecayRadius;
  TH2F* fHistSelV0MassV0PropLifeTime;

  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
