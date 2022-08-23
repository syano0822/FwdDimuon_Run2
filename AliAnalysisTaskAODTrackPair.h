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
    fIsK0sAna = isK0s;
  }
  void setKaonTrackAna(bool isKaon) {
    fIsKaonTrackAna = isKaon;
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
  bool MidTrackPIDChecker(AliAODTrack* track, bool isSel);
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
  bool fIsK0sAna;
  bool fIsKaonTrackAna;
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

  TH2F* fHistTrackEta;
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
  TH2F* fHistTPCSigmaPKaon;
  TH2F* fHistTOFSigmaPKaon;
  TH3F* fHistTPCTOFSigmaKaon;
  TH2F* fHistSelTPCdEdxP;
  TH2F* fHistSelBetaP;      
  TH2F* fHistSelTPCSigmaPKaon;
  TH2F* fHistSelTOFSigmaPKaon;
  TH3F* fHistSelTPCTOFSigmaKaon;
  
  TH2F* fHistArmenteros;
  TH2F* fHistSelArmenteros;

  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
