#ifndef AliAnalysisTaskAODTrackPairEventCounter_cxx
#define AliAnalysisTaskAODTrackPairEventCounter_cxx

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "memory"

class AliAnalysisTaskAODTrackPairEventCounter : public AliAnalysisTaskSE {
  
 public:
  
  AliAnalysisTaskAODTrackPairEventCounter();
  AliAnalysisTaskAODTrackPairEventCounter(const char *name);
  
  virtual ~AliAnalysisTaskAODTrackPairEventCounter();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
  void SetMC(bool isMC){fIsMC = isMC;}

 private:
  
  AliAnalysisTaskAODTrackPairEventCounter(const AliAnalysisTaskAODTrackPairEventCounter&); // not implemented
  AliAnalysisTaskAODTrackPairEventCounter& operator=(const AliAnalysisTaskAODTrackPairEventCounter&); // not implemented
  
  Bool_t InitEvent();  
  Bool_t SingleTrackAnalysis();
  
  Bool_t fIsMC;  
  TList* fOutputList;
  AliPIDResponse      *fPIDResponse;
  AliAODEvent    *fEvent;
  AliAODHandler  *fHandler;
  TClonesArray   *fMCTrackArray; 
  AliMultSelection         *fMultSelection;
  float fCent;
  float fVtxZ;
  float fPsi;

  ////////////////////////////////////////////////
  //Output histos
  ////////////////////////////////////////////////
  
  TH1F* fHistCentrality;
  TH1F* fHistPsi;
  TH1F* fHistAllPt;
  TH1F* fHistPionPt;
  TH1F* fHistKaonPt;


  ClassDef(AliAnalysisTaskAODTrackPairEventCounter, 1); // example of analysis
};

#endif
