#ifndef AliAnalysisTaskAODTrackPairMC_cxx
#define AliAnalysisTaskAODTrackPairMC_cxx

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "TH3F.h"
#include "AliAnalysisTaskAODTrackPairUtils.h"

class AliAnalysisTaskAODTrackPairMC : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskAODTrackPairMC();
  AliAnalysisTaskAODTrackPairMC(const char* name);
  virtual ~AliAnalysisTaskAODTrackPairMC();

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

  AliAnalysisTaskAODTrackPairMC(const AliAnalysisTaskAODTrackPairMC&); // not implemented
  AliAnalysisTaskAODTrackPairMC& operator=(const AliAnalysisTaskAODTrackPairMC&); // not implemented

  bool Initialize();
  bool FwdMuonPairAnalysis();
  bool FwdMuonPairAnalysisEveMixing();
  bool FwdMuonTrackQA(AliAODTrack* track);
  bool MidTrackQA();
  bool FillingRecMuonTree(AliAODTrack* track);
  bool EventQA();
  bool isPrimaryMuonTrack(AliAODMCParticle *particle1);
  bool processMC();

  AliAODEvent    *fEvent;
  AliEventPoolManager *fPoolMuonTrackMgr;
  AliAnalysisTaskAODTrackPairUtils* fUtils;
  TClonesArray   *fMCTrackArray;

  bool fIsMC;
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

  TTree *fTreeML;
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
  float fTrackTPCsigmaPion;
  float fTrackTOFsigmaPion;
  float fTrackTPCsigmaKaon;
  float fTrackTOFsigmaKaon;
  float fTrackTPCsigmaProton;
  float fTrackTOFsigmaProton;
  float fTrackTPCsigmaElectron;
  float fTrackTOFsigmaElectron;
  float fTrackTPCsigmaMuon;
  float fTrackTOFsigmaMuon;
  int fTrackGlobal;
  int fTrackGlobalNoDCA;
  int fTrackTPConly;
  int fTrackPdgCode;
  int fTrackMotherPdgCode;

  TTree* fTreeRecMuonP;
  TTree* fTreeRecMuonN;
  float RecMuonPt;
  float RecMuonEta;
  float RecMuonRap;
  float RecMuonPhi;
  float RecMuonThetaAbs;
  float RecMuonChiSquare;
  float RecMuonTriggerChiSquare;
  int RecMuonTriggerMatch;
  int RecMuonIsGoodTrack;
  
  float RecMCMuonPt;
  float RecMCMuonEta;
  float RecMCMuonRap;
  float RecMCMuonPhi;
  
  TTree* fTreeMCMuonP;
  TTree* fTreeMCMuonN;
  float MCMuonPt;
  float MCMuonEta;
  float MCMuonRap;
  float MCMuonPhi;
  int MCMuonDetect;

  TTree* fTreeULSDimuon;
  TTree* fTreeLSppDimuon;
  TTree* fTreeLSmmDimuon;

  TTree* fTreeMCULSDimuon;
  TTree* fTreeMCLSppDimuon;
  TTree* fTreeMCLSmmDimuon;

  TTree* fTreeMixULSDimuon;
  TTree* fTreeMixLSppDimuon;
  TTree* fTreeMixLSmmDimuon;
  
  float RecDimuonPt;
  float RecDimuonRap;
  float RecDimuonMass;
  float RecDimuonCent;
  float RecDimuonDS;

  float RecMCDimuonPt;
  float RecMCDimuonRap;
  float RecMCDimuonMass;
  float RecMCMomPdgCode;
  float RecMCMomMass;
  float RecMCMomEta;
  float RecMCMomPt;
  
  int RecMCMom2Body;
  int RecMCMomDalitz;

  float MCDimuonPt;
  float MCDimuonRap;
  float MCDimuonMass;
  float MCMomPt;
  float MCMomEta;
  float MCMomPdgCode;
  int MCMom2Body;
  int MCMomDalitz;
  int MCDimuonDetected;

  bool fIsMidMuonAna;
  
  ClassDef(AliAnalysisTaskAODTrackPairMC, 1); // example of analysis
};

#endif
