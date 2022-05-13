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
  void setMCArray(TClonesArray* array){fMCArray = array;}
  void setMC(bool isMC){fIsMC = isMC;}
  void setEvtSelection(bool isEvtSel){fIsEvtSelect=isEvtSel;}
  bool setMCEventInfo();

  bool isAcceptEvent();
  bool isAcceptMuonTrack(AliAODTrack* track);
  bool isAcceptDimuon(AliAODDimuon* dimuon);
  
  bool isSameMotherPair(AliAODTrack* track1,AliAODTrack* track2);
  bool isSameMotherPair(AliAODMCParticle *part1, AliAODMCParticle *part2);
  bool isCharmQuarkOrigin(AliAODMCParticle* particle);
  bool isBeautyQuarkOrigin(AliAODMCParticle* particle);
  bool isHeavyFlavorOrigin(AliAODMCParticle* particle);
  bool isPrimary(AliAODMCParticle* particle);
  int getMotherPdgCode(AliAODTrack *track);
  int getMotherPdgCode(AliAODMCParticle *part);
  int getMotherLabel(AliAODTrack *track);
  int getMotherLabel(AliAODMCParticle *part);
  bool setTrueChPart();
  
  bool getTrueChPartInV0s(int &v0a, int& v0c)
  {
    v0a = fNChV0A;
    v0c = fNChV0C;
    return true;
  }

  int getNTrueChTrkInfo(int spec)
  {
    if (mode == 0) {
      trk = fNChEta05;
    } else if (mode == 1) {
      trk = fNChEta10;
    } else if (mode == 2){
      trk = fNChEta15;
    } else {
      trk = fNChEta20;
    }
    return true;
  }
  

  bool isMC(){
    return fIsMC;
  }
  
  void setDalitzProd(bool flag){
    fIsDalitzProd = flag;
  }
  void set2BodyProd(bool flag){
    fIs2BodyProd = flag;
  }

  bool isDalitzProd(){
    return fIsDalitzProd;
  }
  bool is2BodyProd(){
    return fIs2BodyProd;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  //Set analysis cut flags
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  void setVertexCut(double min, double max, int min_cont)
  {
    fMinVertexCutZ = min;
    fMaxVertexCutZ = max;
    fMinContVtx = min_cont;
    fIsVtxZcut = true;
  }  
  void setPairRapidityCut(double min, double max)
  {
    fMinPairRapCut = min;
    fMaxPairRapCut = max;
    fIsPairRapCut = true;
  }
  void setPairKinematicCut(int type=0,double min=0.){
    if(type==0){
      fIsPairPtCutForOneTrack = false;
      fIsPairPtCutForBothTracks = false;
      fMinPairPtCut = 0.;
    }
    else if(type==1){
      fIsPairPtCutForOneTrack = true;
      fIsPairPtCutForBothTracks = false;
      fMinPairPtCut = min;
    }
    else if(type==2){
      fIsPairPtCutForOneTrack = false;
      fIsPairPtCutForBothTracks = true;
      fMinPairPtCut = min;
    }
    else{
      std::cout<<"Pair pt cut is not correct..."<<std::endl;
      std::cout<<"You set type = "<<type<<"  but it should be between 0 - 2"<<std::endl;
    }
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

  double getTrueVtxZ(){
    return fTrueVtx[2];
  }

  double getVtxZ()
  {
    return fVtxZ;
  }
  double getCentClass(int spec)
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
  double getVtxCont() {
    return fNContVtx;
  }
  double getCentClass()
  {
    return fCent;
  }
  double getPsi()
  {
    return fPsi;
  }
  double getDS()
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

  double getVzeroInfo(int spec)
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
  

  const int fPdgCodeEta = 221;
  const int fPdgCodeRho = 113;
  const int fPdgCodeOmega = 223;
  const int fPdgCodeEtaPrime = 331;
  const int fPdgCodeK0star = 313;
  const int fPdgCodePhi = 333;
   
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
  TClonesArray* fMCArray;

  int fRunNumber;
  int fRunNumberIndex;

  std::string fPeriod;
  std::string fCollSystem;
  std::string fPass;

  bool fIsMC;
  bool fIsEvtSelect;

  bool fIsVtxZcut;
  double fMaxVertexCutZ;
  double fMinVertexCutZ;
  int fNContVtx;
  int fMinContVtx;

  bool fIsPairRapCut;
  double fMinPairRapCut;
  double fMaxPairRapCut;
  
  bool fIsPairPtCutForOneTrack;
  bool fIsPairPtCutForBothTracks;
  double fMinPairPtCut;

  bool fIsPUcut;  
  bool fIsLBCut;

  bool fIsDalitzProd;
  bool fIs2BodyProd;
  
  std::string fMultiMethod;
  
  TH1F* fHistDsCMSL7;
  TH1F* fHistDsCINT7;
  TH1F* fHistDsCMLL7;
  
  double fVtxZ;
  double fCent;
  double fPsi;

  double fTrueVtx[3];

  double fCentSPDTrk;
  double fCentV0A;
  double fCentV0C;
  double fCentV0M;

  double fDSfactor;
  
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

  double fChV0A;
  double fChV0C;
  double fChV0M;
  double fTimeV0A;
  double fTimeV0C;

  int fNChV0A;
  int fNChV0C;

  int fNChEta05;
  int fNChEta10;
  int fNChEta15;
  int fNChEta20;

  ClassDef(AliAnalysisTaskAODTrackPairUtils, 1); // example of analysis
};

#endif
