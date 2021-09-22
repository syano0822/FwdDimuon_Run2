#ifndef AliAnalysisTaskAODTrackPairUtils_cxx
#define AliAnalysisTaskAODTrackPairUtils_cxx

#include "TNamed.h"
#include "TFile.h"
#include "AliMuonTrackCuts.h"
#include "AliVEventHandler.h"
#include "AliAnalysisMuonUtility.h"
#include <iostream>

class AliAnalysisTaskAODTrackPairUtils : public TNamed {
  
 public:
  
  AliAnalysisTaskAODTrackPairUtils();
  ~AliAnalysisTaskAODTrackPairUtils();

  bool setEvent(AliAODEvent* event, AliVEventHandler* handler);

  void setMC(bool isMC){fIsMC = isMC;}
  
  bool isAcceptEvent();
  bool isAcceptMuonTrack(AliAODTrack* track);
  bool isAcceptDimuon(AliAODDimuon* dimuon);
  
  bool isSameMotherPair(TClonesArray* fMCTrackArray,AliAODTrack* track1,AliAODTrack* track2);
  int getMotherPdgCode(TClonesArray* fMCTrackArray,AliAODTrack *track);

  bool isMC(){
    return fIsMC;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  //Set analysis cut flags
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  void setVertexCut(float min, float max)
  {
    fMinVertexCutZ = min;
    fMaxVertexCutZ = max;
    fIsVtxZcut = true;
  }  
  void setPairRapidityCut(float min, float max)
  {
    fMinPairRapCut = -4.0;
    fMaxPairRapCut = -2.5;
    fIsPairRapCut = true;
  }  
  void setPileupRejectionCut(bool flag)
  {
    fIsPUcut = flag;
  }
  void setLocalBoardCut(bool flag)
  {
    fIsLBCut = flag;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //Set analysis object
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  void setMultiEstimateMethod(std::string method)
  {
    fMultiMethod = method;
  }  
  void setMuonTrackCut(AliMuonTrackCuts* cut)
  {
    fMuonTrackCuts = cut;  
  }
  void setDownScalingHist(TFile* inFile)
  {    
    fHistDsCMLL7  = (TH1F*)inFile->Get("DS_MLL7")->Clone("fHistDsCMLL7");
    fHistDsCMSL7  = (TH1F*)inFile->Get("DS_MSL7")->Clone("fHistDsCMSL7");
    fHistDsCINT7  = (TH1F*)inFile->Get("DS_INT7")->Clone("fHistDsCINT7");
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////
  //Get the analysis variables
  //////////////////////////////////////////////////////////////////////////////////////////////

  float getVtxZ()
  {
    return fVtxZ;
  }
  float getCentClass(int spec)
  {
    if(spec==0)
      return fCentSPDTrk;
    else if(spec==1)
      return fCentV0M;
    else if(spec==2)
      return fCentV0A;
    else if(spec==3)
      return fCentV0C;
    else
      return -999;
  }
  float getCentClass()
  {
    return fCent;
  }
  float getPsi()
  {
    return fPsi;
  }
  float getDS()
  {
    return fDSfactor;
  }
  int getRunnumberIndex()
  {
    return fRunNumberIndex;
  }
  int getRunnumber(){
    return fRunNumber;
  }  
  int getNSPDTrkInfo(int spec){
    if(spec==0)
      return fNSPDTrk05;
    else if(spec==1)
      return fNSPDTrk10;
    else if(spec==2)
      return fNSPDTrk15;
    else if(spec==3)
      return fNSPDTrk20;
    else if(spec==4)
      return fNSPDTrkAll;
    else 
      return 0;
  }
  int getNSPDClustInfo(int spec){
    if(spec==0)
      return fNClustSPD1;
    else if(spec==1)
      return fNClustSPD2;
    else 
      return 0;
  }

  float getVzeroInfo(int spec)
  {
    if(spec==0)
      return fChV0A;
    else if(spec==1)
      return fChV0C;
    else if(spec==2)
      return fChV0M;
    else if(spec==3)
      return fTimeV0A-fTimeV0C;
    else if(spec==4)
      return fTimeV0A+fTimeV0C;
    else
      return 0;
  }

  std::string getPass(){
    return fPass;
  }
  void getPeriodInfo(std::string& period, std::string& system)
  {
    period = fPeriod;
    system = fCollSystem;
  }  
  void getTriggerInfo(bool &isCINT7, bool &isCMSL7, bool &isCMSH7, bool &isCMUL7, bool &isCMLL7, bool& is0MSL, bool& is0MSH, bool& is0MUL, bool& is0MLL)
  {
    is0MSL =  fIs0MSL;
    is0MSH =  fIs0MSH;
    is0MUL =  fIs0MUL;
    is0MLL =  fIs0MLL;

    isCINT7 =  fIsCINT7;
    isCMSL7 =  fIsCMSL7;
    isCMSH7 =  fIsCMSH7;
    isCMUL7 =  fIsCMUL7;
    isCMLL7 =  fIsCMLL7;    
  }
  void getTriggerInfo(bool &isCINT7, bool &isCMSL7, bool &isCMSH7, bool &isCMUL7, bool &isCMLL7)
  {
    isCINT7 =  fIsCINT7;
    isCMSL7 =  fIsCMSL7;
    isCMSH7 =  fIsCMSH7;
    isCMUL7 =  fIsCMUL7;
    isCMLL7 =  fIsCMLL7;    
  }
  
  AliMuonTrackCuts* getMuonTrackCuts()
  {
    return fMuonTrackCuts;
  }
  
  //private:
  
  void setInit();
  bool isSameRunnumber();
  bool setVtxZCentPsi();
  bool setDownScaleFactor();
  bool setPeriodInfo();  
  bool setRunnumberIndex();
  bool setTriggerInfo();
  bool setSPDTrk();
  bool setSPDClust();
  bool setVZERO();

  AliAODEvent* fEvent;
  AliMultSelection* fMultSelection;  
  AliMuonTrackCuts* fMuonTrackCuts;
  AliVEventHandler* fInputHandler;


  int fRunNumber;
  int fRunNumberIndex;

  std::string fPeriod;
  std::string fCollSystem;
  std::string fPass;

  bool fIsMC;
 
  bool fIsVtxZcut;
  float fMaxVertexCutZ;
  float fMinVertexCutZ;

  bool fIsPairRapCut;
  float fMinPairRapCut;
  float fMaxPairRapCut;
  
  bool fIsPUcut;  
  bool fIsLBCut;
  
  std::string fMultiMethod;
  
  TH1F* fHistDsCMSL7;
  TH1F* fHistDsCINT7;
  TH1F* fHistDsCMLL7;
  
  float fVtxZ;
  float fCent;
  float fPsi;

  float fCentSPDTrk;
  float fCentV0A;
  float fCentV0C;
  float fCentV0M;

  float fDSfactor;
  
  bool fIs0MSL;
  bool fIs0MSH;
  bool fIs0MUL;
  bool fIs0MLL;

  bool fIsCINT7;
  bool fIsCMSL7;
  bool fIsCMSH7;
  bool fIsCMUL7;
  bool fIsCMLL7;

  int fInput0MSH;
  int fInput0MLL;
  int fInput0MUL;
  int fInput0MSL;
    
  int fNSPDTrk05;
  int fNSPDTrk10;
  int fNSPDTrk15;
  int fNSPDTrk20;
  int fNSPDTrkAll;
  
  int fNClustSPD1;
  int fNClustSPD2;

  float fChV0A;
  float fChV0C;
  float fChV0M;
  float fTimeV0A;
  float fTimeV0C;

  ClassDef(AliAnalysisTaskAODTrackPairUtils, 1); // example of analysis
};

#endif
