#ifndef AliAnalysisTaskAODTrackPair_cxx
#define AliAnalysisTaskAODTrackPair_cxx

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliAnalysisTaskAODTrackPairUtils.h"

class AliAnalysisTaskAODTrackPair : public AliAnalysisTaskSE {
  
 public:
  
  AliAnalysisTaskAODTrackPair();
  AliAnalysisTaskAODTrackPair(const char* name);
  virtual ~AliAnalysisTaskAODTrackPair();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
  void setMC(bool isMC){fIsMC = isMC;}
  
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
  
  bool Initialize();
  bool MuonPairAnalysis();
  bool MuonPairAnalysisEveMixing();  
  bool MuonTrackQA(AliAODTrack* track);
  bool ProcessMC();

  AliAODEvent    *fEvent;  
  AliEventPoolManager *fPoolMuonTrackMgr;
  AliAnalysisTaskAODTrackPairUtils* fUtils;
  TClonesArray   *fMCTrackArray;
  
  bool fIsMC;  
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

  THnSparse* fSparseULSDimuon;
  THnSparse* fSparseLSppDimuon;
  THnSparse* fSparseLSmmDimuon;
  THnSparse* fSparseMixULSDimuon;
  THnSparse* fSparseMixLSppDimuon;
  THnSparse* fSparseMixLSmmDimuon;

  THnSparse* fSparseEtaDalitz;
  THnSparse* fSparseEta2Body;
  THnSparse* fSparseRho2Body;
  THnSparse* fSparseOmegaDalitz;
  THnSparse* fSparseOmega2Body;
  THnSparse* fSparsePhi2Body;
  THnSparse* fSparseEtaPrimeDalitz;

  THnSparse* fSparseEtaDalitzMC;
  THnSparse* fSparseEta2BodyMC;
  THnSparse* fSparseRho2BodyMC;
  THnSparse* fSparseOmegaDalitzMC;
  THnSparse* fSparseOmega2BodyMC;
  THnSparse* fSparsePhi2BodyMC;
  THnSparse* fSparseEtaPrimeDalitzMC;


  TH2F* fHistTrackEta;
  TH2F* fHistTrackThetaAbs;
  TH2F* fHistTrackTriggerMatch;
  TH2F* fHistTrackPDCA;
  TH2F* fHistTrackChiSquare;
  TH2F* fHistTriggerChiSquare;

  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
