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

class KFParticleK0s : public KFParticle {
 
 public:

  bool CheckGetMeasurement( const KFParticleBase& daughter, float m[], float V[], float D[3][3] ) {
    
    float ds[2] = {0.f,0.f};
    float dsdr[4][6];
    float F1[36], F2[36], F3[36], F4[36];
    for(int i1=0; i1<36; i1++){
      F1[i1] = 0;
      F2[i1] = 0;
      F3[i1] = 0;
      F4[i1] = 0;
    }
    //cout<<"       [1.2.3.1]"<<endl;    
    GetDStoParticle( daughter, ds, dsdr );
    //cout<<"       [1.2.3.2]"<<endl;
    for (int index=0; index<2; ++index) if (fabs(ds[index])>safe_threshold) return false;
    /*
    for (int index1=0; index1<4; ++index1) {
      for (int index2=0; index2<6; ++index2) {
	if (fabs(dsdr[index1][index2])>safe_threshold) return false;
      }
    }
    */
    //cout<<"       [1.2.3.3]"<<endl;
    float fdP[8]={0.};
    float fdC[36]={0.};
    
    for(int index=0; index<36; ++index) {
      fdC[index] = daughter.GetCovariance(index);
      if (fabs(fdC[index])>safe_threshold) return false;
    }
    //cout<<"       [1.2.3.4]"<<endl;
    if( fabs(ds[0]*fP[5]) > 1000.f || fabs(ds[1]*daughter.GetParameter(5)) > 1000.f) return false;
    //cout<<"       [1.2.3.5]"<<endl;
    float V0Tmp[36] = {0.};
    float V1Tmp[36] = {0.};
    
    float C[36];
    for(int iC=0; iC<36; iC++) C[iC] = fC[iC];      
    //cout<<"       [1.2.3.6]"<<endl;
    Transport(ds[0], dsdr[0], fP, fC, dsdr[1], F1, F2);
    //cout<<"       [1.2.3.7]"<<endl;
    /*
    for(int i1=0; i1<36; i1++){
      if (fabs(F1[i1])>safe_threshold) return false;
      if (fabs(F2[i1])>safe_threshold) return false;
    }
    */
    //cout<<"       [1.2.3.8]"<<endl;
    daughter.Transport(ds[1], dsdr[3],  m,  V, dsdr[2], F4, F3);
    //cout<<"       [1.2.3.9]"<<endl;
    /*
    for(int i1=0; i1<36; i1++){
      if (fabs(F3[i1])>safe_threshold) return false;
      if (fabs(F4[i1])>safe_threshold) return false;
    }
    */
    //cout<<"       [1.2.3.10]"<<endl;
    MultQSQt(F2, fdC, V0Tmp, 6);   
    MultQSQt(F3, C, V1Tmp, 6);
    //cout<<"       [1.2.3.11]"<<endl;
    ///*
    for(int i1=0; i1<36; i1++){
      if (fabs(V0Tmp[i1])>safe_threshold) return false;
      if (fabs(V1Tmp[i1])>safe_threshold) return false;
    }
    //*/
    //cout<<"       [1.2.3.12]"<<endl;
    for(int iC=0; iC<21; iC++){
      fC[iC] += V0Tmp[iC];
      if (fabs(fC[iC])>safe_threshold) return false;
      V[iC]  += V1Tmp[iC];
      if (fabs(V[iC])>safe_threshold) return false;
    }
    //cout<<"       [1.2.3.13]"<<endl;
    float C1F1T[6][6];
    for(int i=0; i<6; i++){
      for(int j=0; j<6; j++){
	C1F1T[i][j] = 0;
	for(int k=0; k<6; k++){
	  C1F1T[i][j] +=  C[IJ(i,k)] * F1[j*6+k];
	  if (fabs(C1F1T[i][j])>safe_threshold) return false;
	}
      }
    }
    //cout<<"       [1.2.3.14]"<<endl;
    float F3C1F1T[6][6];
    for(int i=0; i<6; i++){
      for(int j=0; j<6; j++){
	F3C1F1T[i][j] = 0;
	for(int k=0; k<6; k++){
	  F3C1F1T[i][j] += F3[i*6+k] * C1F1T[k][j];
	  if (fabs(F3C1F1T[i][j])>safe_threshold) return false;
	}
      }
    }
    //cout<<"       [1.2.3.15]"<<endl;
    float C2F2T[6][6];
    for(int i=0; i<6; i++){
      for(int j=0; j<6; j++){
	C2F2T[i][j] = 0;
	for(int k=0; k<6; k++){
	  C2F2T[i][j] +=  fdC[IJ(i,k)] * F2[j*6+k];
	  if (fabs(C2F2T[i][j])>safe_threshold) return false;
	}
      }
    }
    //cout<<"       [1.2.3.16]"<<endl;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
	D[i][j] = F3C1F1T[i][j];
	for(int k=0; k<6; k++){
	  D[i][j] += F4[i*6+k] * C2F2T[k][j];
	  if (fabs(D[i][j])>safe_threshold) return false;
	}
      }    
    }
    //cout<<"       [1.2.3.17]"<<endl;
    return 1;
    
  }

  
  Bool_t CheckIfSafePair(KFParticle Daughter) {
    //cout<<"    [1.2.1.0]"<<endl;
    for ( int index =0; index<28; ++index) if (fabs(fC[index])>safe_threshold) return false;
    //cout<<"    [1.2.2.0]"<<endl;
    Int_t maxIter = 1;
    
    for( Int_t iter=0; iter<maxIter; iter++ ){
      
      float m[8], mV[36];
    
      float D[3][3];
      
      //cout<<"    [1.2.3.0]"<<endl;
      if(! CheckGetMeasurement(Daughter, m, mV, D) ) return false;
      //cout<<"    [1.2.4.0]"<<endl;
      //for ( int index =0; index<36; ++index) if (fabs(mV[index])>safe_threshold) return false;
      /*
      for ( int index1 =0; index1<3; ++index1){
	for ( int index2 =0; index2<3; ++index2){
	  if (fabs(D[index1][index2])>safe_threshold) return false;
	}
      }
      */
      //cout<<"    [1.2.5.0]"<<endl;
      float mS[6]= { fC[0]+mV[0], 
		     fC[1]+mV[1], fC[2]+mV[2], 
		     fC[3]+mV[3], fC[4]+mV[4], fC[5]+mV[5] };
      //cout<<"    [1.2.6.0]"<<endl;
      for ( int index =0; index<6; ++index) if (fabs(mS[index])>safe_threshold) return false;
      //cout<<"    [1.2.7.0]"<<endl;
      InvertCholetsky3(mS);
      //cout<<"    [1.2.8.0]"<<endl;
      //* Residual (measured - estimated)
      
      float zeta[3] = { m[0]-fP[0], m[1]-fP[1], m[2]-fP[2] };    
      //cout<<"    [1.2.9.0]"<<endl;
      for ( int index =0; index<3; ++index) if (fabs(zeta[index])>safe_threshold) return false;
      //cout<<"    [1.2.10.0]"<<endl;
      float dChi2 = (mS[0]*zeta[0] + mS[1]*zeta[1] + mS[3]*zeta[2])*zeta[0]
	+      (mS[1]*zeta[0] + mS[2]*zeta[1] + mS[4]*zeta[2])*zeta[1]
	+      (mS[3]*zeta[0] + mS[4]*zeta[1] + mS[5]*zeta[2])*zeta[2]; 
      if (dChi2 > 1e9) return false;
      //cout<<"    [1.2.11.0]"<<endl;
      float K[3][3];
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  K[i][j] = 0;
	  for(int k=0; k<3; k++){
	    K[i][j] += fC[IJ(i,k)] * mS[IJ(k,j)];
	    if (fabs(K[i][j])>safe_threshold) return false;
	  }
	}
      }
      //cout<<"    [1.2.12.0]"<<endl;
      //* CHt = CH' - D'
      float mCHt0[7], mCHt1[7], mCHt2[7];

      mCHt0[0]=fC[ 0] ;       mCHt1[0]=fC[ 1] ;       mCHt2[0]=fC[ 3] ;
      mCHt0[1]=fC[ 1] ;       mCHt1[1]=fC[ 2] ;       mCHt2[1]=fC[ 4] ;
      mCHt0[2]=fC[ 3] ;       mCHt1[2]=fC[ 4] ;       mCHt2[2]=fC[ 5] ;
      mCHt0[3]=fC[ 6]-mV[ 6]; mCHt1[3]=fC[ 7]-mV[ 7]; mCHt2[3]=fC[ 8]-mV[ 8];
      mCHt0[4]=fC[10]-mV[10]; mCHt1[4]=fC[11]-mV[11]; mCHt2[4]=fC[12]-mV[12];
      mCHt0[5]=fC[15]-mV[15]; mCHt1[5]=fC[16]-mV[16]; mCHt2[5]=fC[17]-mV[17];
      mCHt0[6]=fC[21]-mV[21]; mCHt1[6]=fC[22]-mV[22]; mCHt2[6]=fC[23]-mV[23];
      //cout<<"    [1.2.13.0]"<<endl;
      ///*
      for ( int index =0; index<7; ++index) {
	if (fabs(mCHt0[index])>safe_threshold) return false;
	if (fabs(mCHt1[index])>safe_threshold) return false;
	if (fabs(mCHt2[index])>safe_threshold) return false;
      }
      //*/
      //* Kalman gain K = mCH'*S
      
      float k0[7], k1[7], k2[7];
      //cout<<"    [1.2.14.0]"<<endl;
      for(Int_t i=0;i<7;++i){
	k0[i] = mCHt0[i]*mS[0] + mCHt1[i]*mS[1] + mCHt2[i]*mS[3];
	k1[i] = mCHt0[i]*mS[1] + mCHt1[i]*mS[2] + mCHt2[i]*mS[4];
	k2[i] = mCHt0[i]*mS[3] + mCHt1[i]*mS[4] + mCHt2[i]*mS[5];
	if (fabs(k0[i])>safe_threshold) return false;
	if (fabs(k1[i])>safe_threshold) return false;
	if (fabs(k2[i])>safe_threshold) return false;
      }
      //cout<<"    [1.2.15.0]"<<endl;
      //* Add the daughter momentum to the particle momentum
      
      fP[ 3] += m[ 3];
      fP[ 4] += m[ 4];
      fP[ 5] += m[ 5];
      fP[ 6] += m[ 6];
      //cout<<"    [1.2.16.0]"<<endl;
      /*
      for(Int_t i=0;i<7;++i) {
	if (fabs(fP[i])>safe_threshold) return false;
      }
      */
      fC[ 9] += mV[ 9];
      fC[13] += mV[13];
      fC[14] += mV[14];
      fC[18] += mV[18];
      fC[19] += mV[19];
      fC[20] += mV[20];
      fC[24] += mV[24];
      fC[25] += mV[25];
      fC[26] += mV[26];
      fC[27] += mV[27];
      //cout<<"    [1.2.17.0]"<<endl;
      /*
      for ( int index =0; index<36; ++index) {	
	if (fabs(fC[index])>safe_threshold) return false;      
      }
      */
      //cout<<"    [1.2.18.0]"<<endl;
      //* New estimation of the vertex position r += K*zeta
    
      for(Int_t i=0;i<7;++i) {
	fP[i] = fP[i] + k0[i]*zeta[0] + k1[i]*zeta[1] + k2[i]*zeta[2];
	if (fabs(fP[i])>safe_threshold) return false;
      }
      //cout<<"    [1.2.19.0]"<<endl;
      //* New covariance matrix C -= K*(mCH')'

      for(Int_t i=0, k=0;i<7;++i){
	for(Int_t j=0;j<=i;++j,++k){
	  fC[k] = fC[k] - (k0[i]*mCHt0[j] + k1[i]*mCHt1[j] + k2[i]*mCHt2[j] );
	  if (fabs(fC[k])>safe_threshold) return false;
	}
      }
      //cout<<"    [1.2.20.0]"<<endl;
      float K2[3][3];
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  K2[i][j] = -K[j][i];
	  if (fabs(K2[i][j])>safe_threshold) return false;
	}
	K2[i][i] += 1;
      }
      //cout<<"    [1.2.21.0]"<<endl;
      float A[3][3];
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  A[i][j] = 0;
	  for(int k=0; k<3; k++){
	    A[i][j] += D[i][k] * K2[k][j];
	    if (fabs(A[i][j])>safe_threshold) return false;
	  }
	}
      }
      //cout<<"    [1.2.21.0]"<<endl;
      double M[3][3];
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  M[i][j] = 0;
	  for(int k=0; k<3; k++){
	    M[i][j] += K[i][k] * A[k][j];
	    if (fabs(M[i][j])>safe_threshold) return false;
	  }
	}
      }
      //cout<<"    [1.2.22.0]"<<endl;
      fC[0] += 2*M[0][0];
      fC[1] += M[0][1] + M[1][0];
      fC[2] += 2*M[1][1];
      fC[3] += M[0][2] + M[2][0];
      fC[4] += M[1][2] + M[2][1];
      fC[5] += 2*M[2][2];
      //cout<<"    [1.2.23.0]"<<endl;
      for ( int index =0; index<36; ++index) if (fabs(fC[index])>safe_threshold) return false;

      //* Calculate Chi^2 
      //cout<<"    [1.2.24.0]"<<endl;
      fNDF  += 2;
      fQ    +=  Daughter.GetQ();
      fSFromDecay = 0;    
      fChi2 += dChi2;    
      //cout<<"    [1.2.25.0]"<<endl;
      if (fabs(fChi2)>safe_threshold) return false;
      //cout<<"    [1.2.26.0]"<<endl;
    }
  }

  double safe_threshold = 1.e2;
  
};

class K0sContainer : public TObject{ 
 public:
  
  K0sContainer(double mass, double eta, double mom[3], double svtx[3], int q, double decay_lengthxy, double decay_length, double d0xy, double d0, 
	       double pointxy, double point, double distance_daughtersxy, double distance_daughters, double chi2, double armenteros[2],
	       int did[2], int label[2], double dmom1[3], double dmom2[3], double dd0xy[2], double dd0[2]){
    
    Mass   = mass;
    Eta    = eta;
    P[0]   = mom[0];
    P[1]   = mom[1];
    P[2]   = mom[2];
    Vtx[0] = svtx[0];
    Vtx[1] = svtx[1];
    Vtx[2] = svtx[2];
    Charge = q;
    DcaXY  = d0xy;
    Dca    = d0;
    DecayLengthXY  = decay_lengthxy;
    DecayLength    = decay_length;
    CpvXY          = pointxy;
    Cpv            = point;
    DistanceDaughtersXY  = distance_daughtersxy;
    DistanceDaughters    = distance_daughters;
    rChi2V0              = chi2;
    Armenteros[0]        = armenteros[0];       
    Armenteros[1]        = armenteros[1];       
    dP1[0]      = dmom1[0];
    dP1[1]      = dmom1[1];
    dP1[2]      = dmom1[2];
    dP2[0]      = dmom2[0];
    dP2[1]      = dmom2[1];
    dP2[2]      = dmom2[2];
    dDcaXY[0]   = dd0xy[0];
    dDcaXY[1]   = dd0xy[1];
    dDca[0]     = dd0[0];
    dDca[1]     = dd0[1];
    dLabel[0]   = label[0];
    dLabel[1]   = label[1];
    dTrackid[0] = did[0];
    dTrackid[1] = did[1];
  }
  
  int Charge;
  int dLabel[2];
  int dTrackid[2];

  double Mass;
  double Eta;
  double P[3];  
  double Vtx[3];
  double DcaXY;
  double Dca;
  double DecayLengthXY;
  double DecayLength;
  double CpvXY;
  double Cpv;
  double DistanceDaughtersXY;
  double DistanceDaughters;
  double rChi2V0;
  double Armenteros[2];
  double dP1[3];
  double dP2[3];
  double dDcaXY[2];
  double dDca[2];
  
};

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
  void setUtils(AliAnalysisTaskAODTrackPairUtils *utils) { fUtils = utils; }
  void setEvtMixingTrackDepth(int depth) { fTrackDepth = depth; }
  void setEvtMixingPoolSize(int size) { fPoolSize = size; }
  void setEvtMixingReadyFraction(double frac) { fReadyFraction = frac; }
  void setEvtMixingPoolVtxZ(bool flag) { onEvtMixingPoolVtxZ = flag; }
  void setEvtMixingPoolCent(bool flag) { onEvtMixingPoolCent = flag; }
  void setEvtMixingPoolPsi(bool flag) { onEvtMixingPoolPsi = flag; }

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
  bool V0QualityChecker(K0sContainer* v0, bool isSel);
  bool V0QualityChecker(AliAODv0* v0, bool isSel);
  bool V0QualityChecker(double mass, double p, int charge, double DCApairXY, double DCApair, double decay_length, double decay_length_xy, double cpv,
			double cpv_xy, double dca, double dca_xy, double distance_daughter, double distance_daughter_xy,
			double chi2, double lifetime, double armenteros_alpha, double armenteros_qt, int dlabel[2], bool isSel);

  bool K0sPairAnalysis(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  bool K0sPairAnalysisEventMixing(AliPID::EParticleType pid1,AliPID::EParticleType pid2);

  bool AddK0sArray(AliPID::EParticleType pid1, AliPID::EParticleType pid2);

  unique_ptr<K0sContainer> calcK0sFromTracks(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2);
  unique_ptr<AliAODv0> calcAliAODv0FromTracks(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2);

  bool updateAODv0(AliAODv0* v0);

  int findLeadingTrack();
  
  KFParticle CreateKFParticle(AliAODTrack* track,int pdg);
  
  K0sContainer *CreateRawPointerK0sContainer(K0sContainer* unique);

  double DecayLengthFromKF(KFParticle kfpParticle, KFParticle PV);
  double DecayLengthXYFromKF(KFParticle kfpParticle, KFParticle PV);

  bool ProcessMC();

  bool isAcceptK0sPair(K0sContainer* v0_1, K0sContainer* v0_2);
  bool isAcceptV0QualityCuts(K0sContainer* v0, bool isAcceptCuts[7]);
  bool isAcceptV0QualityCuts(AliAODv0* v0, bool isAcceptCuts[7]);
  bool isAcceptV0QualityCuts(double p, double chi2, double cpv, double lengthXY, double trackDistanceXY, double DCAxy, double lifetime, bool isAcceptCuts[7]);
  bool isAcceptK0sArmenteros(K0sContainer* v0);
  bool isAcceptK0sCuts(K0sContainer* v0);
  bool isAcceptK0sKinematics(K0sContainer* v0);
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
  
  AliEventPoolManager *fPoolMuonTrackMgrK0s;  
  AliEventPoolManager *fPoolMuonTrackMgrPion;  
  
  AliMultSelection *fMultSelection;
  
  AliAnalysisTaskAODTrackPairUtils *fUtils;
  
  TClonesArray *fMCTrackArray;

  unique_ptr<TObjArray> fArrayK0s;
  unique_ptr<TObjArray> fArrayK0sDaughterTrack;

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
  
  THnSparse *fSparseULSPionPair;
  THnSparse *fSparseLSppPionPair;
  THnSparse *fSparseLSmmPionPair;

  THnSparse *fSparseMixULSPionPair;
  THnSparse *fSparseMixLSppPionPair;
  THnSparse *fSparseMixLSmmPionPair;

  THnSparse *fSparseULSPionPairBeforeCuts;
  THnSparse *fSparseLSppPionPairBeforeCuts;
  THnSparse *fSparseLSmmPionPairBeforeCuts;
  
  THnSparse* fSparseULSPionPair_PassChi2perNDFCut;
  THnSparse* fSparseULSPionPair_PassCPVXYCut;
  THnSparse* fSparseULSPionPair_PassDecayLengthXYCut;
  THnSparse* fSparseULSPionPair_PassDaughterDistanceXYCut;  
  THnSparse* fSparseULSPionPair_PassDCAXYCut;
  THnSparse* fSparseULSPionPair_PassLifetimeCut;
  THnSparse* fSparseULSPionPair_PassArmenterosCut;

  THnSparse *fSparseNeutralK0sPair;
  THnSparse *fSparseNeutralNegativeK0sPair;
  THnSparse *fSparseNeutralPositiveK0sPair;
  THnSparse *fSparsePositiveK0sPair;
  THnSparse *fSparseNegativeK0sPair;
  THnSparse *fSparsePositiveNegativeK0sPair;

  THnSparse *fSparseTrueK0s;
  THnSparse *fSparseTrueK0sRecK0s;

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

  TH2F *fHistLSV0PV0DecayLengthXY;
  TH2F *fHistLSV0PV0DecayLength;
  TH2F *fHistLSV0PV0PointingAngleXY;
  TH2F *fHistLSV0PV0PointingAngle;
  TH2F *fHistLSV0PV0DCAXY;
  TH2F *fHistLSV0PV0DCA;
  TH2F *fHistLSV0PV0TrackDistanceXY;
  TH2F *fHistLSV0PV0TrackDistance;
  TH2F *fHistLSV0PV0Chi2perNDF;
  TH2F *fHistLSV0PV0PropLifeTime;
  TH2F *fHistLSV0TrackDCAXY; 
  TH2F *fHistLSV0TrackDCA; 
  TH2F *fHistLSV0Armenteros;

  TH2F *fHistLSSelV0PV0DecayLengthXY;
  TH2F *fHistLSSelV0PV0DecayLength;
  TH2F *fHistLSSelV0PV0PointingAngleXY;
  TH2F *fHistLSSelV0PV0PointingAngle;
  TH2F *fHistLSSelV0PV0DCAXY;
  TH2F *fHistLSSelV0PV0DCA;
  TH2F *fHistLSSelV0PV0TrackDistanceXY;
  TH2F *fHistLSSelV0PV0TrackDistance;
  TH2F *fHistLSSelV0PV0Chi2perNDF;
  TH2F *fHistLSSelV0PV0PropLifeTime;
  TH2F *fHistLSSelV0TrackDCAXY; 
  TH2F *fHistLSSelV0TrackDCA; 
  TH2F *fHistLSSelV0Armenteros;

  TH2F *fHistLSV0PV0DecayLengthXY_TrueK0s;
  TH2F *fHistLSV0PV0DecayLength_TrueK0s;
  TH2F *fHistLSV0PV0PointingAngleXY_TrueK0s;
  TH2F *fHistLSV0PV0PointingAngle_TrueK0s;
  TH2F *fHistLSV0PV0DCAXY_TrueK0s;
  TH2F *fHistLSV0PV0DCA_TrueK0s;
  TH2F *fHistLSV0PV0TrackDistanceXY_TrueK0s;
  TH2F *fHistLSV0PV0TrackDistance_TrueK0s;
  TH2F *fHistLSV0PV0Chi2perNDF_TrueK0s;
  TH2F *fHistLSV0PV0PropLifeTime_TrueK0s;
  TH2F *fHistLSV0TrackDCAXY_TrueK0s; 
  TH2F *fHistLSV0TrackDCA_TrueK0s; 
  TH2F *fHistLSV0Armenteros_TrueK0s;

  TH2F *fHistLSSelV0PV0DecayLengthXY_TrueK0s;
  TH2F *fHistLSSelV0PV0DecayLength_TrueK0s;
  TH2F *fHistLSSelV0PV0PointingAngleXY_TrueK0s;
  TH2F *fHistLSSelV0PV0PointingAngle_TrueK0s;
  TH2F *fHistLSSelV0PV0DCAXY_TrueK0s;
  TH2F *fHistLSSelV0PV0DCA_TrueK0s;
  TH2F *fHistLSSelV0PV0TrackDistanceXY_TrueK0s;
  TH2F *fHistLSSelV0PV0TrackDistance_TrueK0s;
  TH2F *fHistLSSelV0PV0Chi2perNDF_TrueK0s;
  TH2F *fHistLSSelV0PV0PropLifeTime_TrueK0s;
  TH2F *fHistLSSelV0TrackDCAXY_TrueK0s; 
  TH2F *fHistLSSelV0TrackDCA_TrueK0s; 
  TH2F *fHistLSSelV0Armenteros_TrueK0s;

  int selectEvt;
  int all;
  int abord;

  clock_t start;
  clock_t end;
  
  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
