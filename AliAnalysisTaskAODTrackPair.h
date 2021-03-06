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

  void setMidMuonAna(bool isMidMuon){
    fIsMidMuonAna = isMidMuon;
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
  
  bool MidMuonTrackQA(AliAODTrack* track);
  bool MidMuonPairQA(AliAODDimuon* dimuon);
  bool MidMuonPairAnalysis();
  bool MidMuonPairAnalysisEveMixing();
  
  AliAODEvent    *fEvent;
  AliEventPoolManager *fPoolMuonTrackMgr;
  AliAnalysisTaskAODTrackPairUtils* fUtils;

  bool fIsMC;
  bool fIsMidMuonAna;
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

  TTree* fTreeULSDimuon;
  TTree* fTreeLSppDimuon;
  TTree* fTreeLSmmDimuon;

  TTree* fTreeMixULSDimuon;
  TTree* fTreeMixLSppDimuon;
  TTree* fTreeMixLSmmDimuon;

  float RecDimuonPt;
  float RecDimuonRap;
  float RecDimuonMass;
  float RecDimuonCent;
  float RecDimuonDS;

  TTree *fTreeMidMuon;
  float fTrackPt;
  float fTrackP;
  float fTrackTheta;
  float fTrackPhi;
  float fTrackLength;
  float fTrackBeta;
  float fTrackTrackChi2perNDF;
  float fTrackTrackITSNcls;
  float fTrackTrackTPCNcls;
  float fTrackTrackTOFNcls;
  float fTrackTrackTPCChi2;
  float fTrackTrackITSChi2;
  float fTrackTPCCrossedRows;
  float fTrackTPCFindableNcls;
  float fTrackTOFBCTime;
  float fTrackTOFKinkIndex;
  float fTrackDCAxy;
  float fTrackDCAz;
  float fTrackTPCsigmaMuon;
  float fTrackTOFsigmaMuon;
  bool fTrackTrueMuonLabel;

  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
