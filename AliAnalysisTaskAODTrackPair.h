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

class KFParticleK0s : public KFParticle {
 public:

  bool SetMassConstraintPrivate( float *mP, float *mC, float mJ[7][7], float mass ){
    
    //* Set nonlinear mass constraint (Mass) on the state vector mP with a covariance matrix mC.
    //cout<<"[1.13.1.0]"<<endl;
    const float energy2 = mP[6]*mP[6], p2 = mP[3]*mP[3]+mP[4]*mP[4]+mP[5]*mP[5], mass2 = mass*mass;
    //cout<<"[1.13.2.0]"<<endl;
    const float a = energy2 - p2 + 2.*mass2;
    //cout<<"[1.13.3.0]"<<endl;
    const float b = -2.*(energy2 + p2);
    //cout<<"[1.13.4.0]"<<endl;
    const float c = energy2 - p2 - mass2;
    //cout<<"[1.13.5.0]"<<endl;
    float lambda = 0;
    if(fabs(b) > 1.e-10) lambda = -c / b;
    
    cout<<"[1.13.6.0]"<<endl;
    cout<<"  energy2="<<energy2<<"  p2="<<p2<<"  mass2="<<mass2<<"  energy2-p2-2.*mass2="<<(energy2-p2-2.*mass2)<<endl;
    cout<<"[1.13.7.0]"<<endl;
    cout<<"  fabs(TMath::Log10(4.*energy2*p2 - mass2*(energy2-p2-2.*mass2)))="<<fabs(TMath::Log10(4.*energy2*p2 - mass2*(energy2-p2-2.*mass2)))<<endl;
    
    if (fabs(TMath::Log10(4.*energy2*p2 - mass2*(energy2-p2-2.*mass2))) > 30) {
      cout<<"    Dont pass fabs(TMath::Log10(4.*energy2*p2 - mass2*(energy2-p2-2.*mass2))) "<<fabs(TMath::Log10(4.*energy2*p2 - mass2*(energy2-p2-2.*mass2)))<<endl;
      return false;
    } 
    
    float d = 4.*energy2*p2 - mass2*(energy2-p2-2.*mass2);
    //cout<<"[1.13.8.0]"<<endl;
    if(d>=0 && fabs(a) > 1.e-10) lambda = (energy2 + p2 - sqrt(d))/a;
    //cout<<"[1.13.9.0]"<<endl;
    if(mP[6] < 0) //If energy < 0 we need a lambda < 0
      lambda = -1000000.; //Empirical, a better solution should be found
    //cout<<"[1.13.10.0]"<<endl;
    Int_t iIter=0;
    for(iIter=0; iIter<100; iIter++){
      float lambda2 = lambda*lambda;
      float lambda4 = lambda2*lambda2;
      
      float lambda0 = lambda;

      float f  = -mass2 * lambda4 + a*lambda2 + b*lambda + c;
      float df = -4.*mass2 * lambda2*lambda + 2.*a*lambda + b;
      if(fabs(df) < 1.e-10) break;
      lambda -= f/df;
      if(fabs(lambda0 - lambda) < 1.e-8) break;
    }    
    //cout<<"[1.13.11.0]  lambda = "<<lambda<<endl;
    const float lpi = 1./(1. + lambda);
    const float lmi = 1./(1. - lambda);
    const float lp2i = lpi*lpi;
    const float lm2i = lmi*lmi;
    //cout<<"[1.13.12.0]  lp2i="<<lp2i<<"  lm2i="<<lm2i<<endl;
    float lambda2 = lambda*lambda;
    //cout<<"[1.13.13.0]"<<endl;
    float dfl  = -4.*mass2 * lambda2*lambda + 2.*a*lambda + b;
    //cout<<"[1.13.14.0]  dfl="<<dfl<<endl;
    float dfx[7] = {0};//,0,0,0};
    dfx[0] = -2.*(1. + lambda)*(1. + lambda)*mP[3];
    dfx[1] = -2.*(1. + lambda)*(1. + lambda)*mP[4];
    dfx[2] = -2.*(1. + lambda)*(1. + lambda)*mP[5];
    dfx[3] = 2.*(1. - lambda)*(1. - lambda)*mP[6];
    //cout<<"[1.13.15.0]"<<endl;
    float dlx[4] = {1,1,1,1};
    
    if(fabs(dfl) > 1.e-10 ){
      for(int i=0; i<4; i++)
	dlx[i] = -dfx[i] / dfl;
    }
    //cout<<"[1.13.16.0]"<<endl;
    float dxx[4] = {mP[3]*lm2i, mP[4]*lm2i, mP[5]*lm2i, -mP[6]*lp2i};
    //cout<<"[1.13.17.0]"<<endl;
    for(Int_t i=0; i<7; i++)
      for(Int_t j=0; j<7; j++)
	mJ[i][j]=0;
    //cout<<"[1.13.18.0]"<<endl;
    mJ[0][0] = 1.;
    mJ[1][1] = 1.;
    mJ[2][2] = 1.;
    //cout<<"[1.13.19.0]"<<endl;
    for(Int_t i=3; i<7; i++)
      for(Int_t j=3; j<7; j++)
	mJ[i][j] = dlx[j-3]*dxx[i-3];

    for(Int_t i=3; i<6; i++)
      mJ[i][i] += lmi;
    mJ[6][6] += lpi;
    //cout<<"[1.13.20.0]"<<endl;
    float mCJ[7][7];

    for(Int_t i=0; i<7; i++) {
      for(Int_t j=0; j<7; j++) {
	mCJ[i][j] = 0;
	for(Int_t k=0; k<7; k++) {
	  mCJ[i][j] += mC[IJ(i,k)]*mJ[j][k];
	}
      }
    }
    //cout<<"[1.13.21.0]"<<endl;
    for(Int_t i=0; i<7; ++i){
      for(Int_t j=0; j<=i; ++j){
	mC[IJ(i,j)]=0;
	for(Int_t l=0; l<7; l++){
	  mC[IJ(i,j)] += mJ[i][l]*mCJ[l][j];
	}
      }
    }
    //cout<<"[1.13.22.0]"<<endl;
    //cout<<"  lmi="<<lmi<<"  lpi="<<lpi<<endl;
    mP[3] *= lmi;
    mP[4] *= lmi;
    mP[5] *= lmi;
    mP[6] *= lpi;
    //cout<<"[1.13.23.0]  mP[0]="<<mP[0]<<"  mP[1]="<<mP[1]<<"  mP[2]="<<mP[2]<<"  mP[3]="<<mP[3]<<"  mP[4]="<<mP[4]<<"  mP[5]="<<mP[5]<<"  mP[6]="<<mP[6]<<endl;

    return true;
  }
  
  Bool_t CheckIfSafePair(KFParticle Daughter) {
    
    Float_t m[8], mV[36], D[3][3];
    if (KFParticleBase::GetMeasurement(Daughter, m, mV, D)) {
      
      float mS[6]= { fC[0]+mV[0], 
		     fC[1]+mV[1], fC[2]+mV[2], 
		     fC[3]+mV[3], fC[4]+mV[4], fC[5]+mV[5] };
      
      cout<<"Original Covariance matrix"<<endl;
      for (int index=0; index<5; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=5; index<10; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=10; index<15; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=16; index<20; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=20; index<25; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=25; index<28; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;

      if (fC[0]>100 || fC[1]>100 || fC[2]>100 || fC[3]>100 || fC[4]>100 || fC[5]>100) {
	cout<<"Don't pass fCs"<<endl;
	return false;
      }
      
      InvertCholetsky3(mS);

      //* Residual (measured - estimated)
      
      float zeta[3] = { m[0]-fP[0], m[1]-fP[1], m[2]-fP[2] };    
      
      float K[3][6];
      for(int i=0; i<3; i++){
	for(int j=0; j<3; j++){
	  K[i][j] = 0;
	  for(int k=0; k<3; k++){
	    K[i][j] += fC[IJ(i,k)] * mS[IJ(k,j)];
	  }
	}
      }
    
      //* CHt = CH'
      
      float mCHt0[7], mCHt1[7], mCHt2[7];
    
      mCHt0[0]=fC[ 0] ; mCHt1[0]=fC[ 1] ; mCHt2[0]=fC[ 3] ;
      mCHt0[1]=fC[ 1] ; mCHt1[1]=fC[ 2] ; mCHt2[1]=fC[ 4] ;
      mCHt0[2]=fC[ 3] ; mCHt1[2]=fC[ 4] ; mCHt2[2]=fC[ 5] ;
      mCHt0[3]=fC[ 6] ; mCHt1[3]=fC[ 7] ; mCHt2[3]=fC[ 8] ;
      mCHt0[4]=fC[10] ; mCHt1[4]=fC[11] ; mCHt2[4]=fC[12] ;
      mCHt0[5]=fC[15] ; mCHt1[5]=fC[16] ; mCHt2[5]=fC[17] ;
      mCHt0[6]=fC[21] ; mCHt1[6]=fC[22] ; mCHt2[6]=fC[23] ;
  
      
      //* Kalman gain K = mCH'*S
      
      float k0[7], k1[7], k2[7];
    
      for(Int_t i=0;i<7;++i){
	k0[i] = mCHt0[i]*mS[0] + mCHt1[i]*mS[1] + mCHt2[i]*mS[3];
	k1[i] = mCHt0[i]*mS[1] + mCHt1[i]*mS[2] + mCHt2[i]*mS[4];
	k2[i] = mCHt0[i]*mS[3] + mCHt1[i]*mS[4] + mCHt2[i]*mS[5];
      }
      
      // last itearation -> update the particle

      //* VHt = VH'
    
      float mVHt0[7], mVHt1[7], mVHt2[7];
    
      mVHt0[0]=mV[ 0] ; mVHt1[0]=mV[ 1] ; mVHt2[0]=mV[ 3] ;
      mVHt0[1]=mV[ 1] ; mVHt1[1]=mV[ 2] ; mVHt2[1]=mV[ 4] ;
      mVHt0[2]=mV[ 3] ; mVHt1[2]=mV[ 4] ; mVHt2[2]=mV[ 5] ;
      mVHt0[3]=mV[ 6] ; mVHt1[3]=mV[ 7] ; mVHt2[3]=mV[ 8] ;
      mVHt0[4]=mV[10] ; mVHt1[4]=mV[11] ; mVHt2[4]=mV[12] ;
      mVHt0[5]=mV[15] ; mVHt1[5]=mV[16] ; mVHt2[5]=mV[17] ;
      mVHt0[6]=mV[21] ; mVHt1[6]=mV[22] ; mVHt2[6]=mV[23] ;
  
      //* Kalman gain Km = mCH'*S
      
      float km0[7], km1[7], km2[7];
      //cout<<"[1.2.0.0]"<<endl;
      for(Int_t i=0;i<7;++i){
	km0[i] = mVHt0[i]*mS[0] + mVHt1[i]*mS[1] + mVHt2[i]*mS[3];
	km1[i] = mVHt0[i]*mS[1] + mVHt1[i]*mS[2] + mVHt2[i]*mS[4];
	km2[i] = mVHt0[i]*mS[3] + mVHt1[i]*mS[4] + mVHt2[i]*mS[5];
      }
      //cout<<"[1.3.0.0]"<<endl;
      for(Int_t i=0;i<7;++i) {
	fP[i] = fP[i] + k0[i]*zeta[0] + k1[i]*zeta[1] + k2[i]*zeta[2];
      }
      //cout<<"[1.4.0.0]"<<endl;
      
      //cout<<"[1.5.0.0]"<<endl;
      for(Int_t i=0;i<7;++i) 
	m[i] = m[i] - km0[i]*zeta[0] - km1[i]*zeta[1] - km2[i]*zeta[2];
      //cout<<"[1.6.0.0]"<<endl;
      for(Int_t i=0, k=0;i<7;++i){
	for(Int_t j=0;j<=i;++j,++k){
	  fC[k] = fC[k] - (k0[i]*mCHt0[j] + k1[i]*mCHt1[j] + k2[i]*mCHt2[j] );
	}
      }
      //cout<<"[1.7.0.0]"<<endl;
      for(Int_t i=0, k=0;i<7;++i){
	for(Int_t j=0;j<=i;++j,++k){
	  mV[k] = mV[k] - (km0[i]*mVHt0[j] + km1[i]*mVHt1[j] + km2[i]*mVHt2[j] );
	}
      }
      //cout<<"[1.8.0.0]"<<endl;
      float mDf[7][7];

      for(Int_t i=0;i<7;++i){
	for(Int_t j=0;j<7;++j){
	  mDf[i][j] = (km0[i]*mCHt0[j] + km1[i]*mCHt1[j] + km2[i]*mCHt2[j] );
	}
      }
      //cout<<"[1.9.0.0]"<<endl;
      float mJ1[7][7], mJ2[7][7];
      for(Int_t iPar1=0; iPar1<7; iPar1++){
	for(Int_t iPar2=0; iPar2<7; iPar2++){
	  mJ1[iPar1][iPar2] = 0;
	  mJ2[iPar1][iPar2] = 0;
	}
      }
      //cout<<"[1.10.0.0]"<<endl;
      float mMassParticle  = fP[6]*fP[6] - (fP[3]*fP[3] + fP[4]*fP[4] + fP[5]*fP[5]);
      float mMassDaughter  = m[6]*m[6] - (m[3]*m[3] + m[4]*m[4] + m[5]*m[5]);
      
      cout<<"Original  fP[0]="<<fP[0]<<"  fP[1]="<<fP[1]<<"  fP[2]="<<fP[2]<<"  fP[3]="<<fP[3]<<
	"  fP[4]="<<fP[4]<<"  fP[5]="<<fP[5]<<"  fP[6]="<<fP[6]<<"  fP[7]="<<fP[7]<<endl;
      
      //cout<<"[1.11.0.0]"<<endl;
      if(mMassParticle > 0) mMassParticle = sqrt(mMassParticle);
      if(mMassDaughter > 0) mMassDaughter = sqrt(mMassDaughter);
      //cout<<"[1.12.0.0]"<<endl;
      if( fMassHypo > -0.5){
	//cout<<"[1.13.0.0   fMassHypo > -0.5]"<<endl;
	if (!SetMassConstraintPrivate(fP,fC,mJ1,SumDaughterMass)) return false;
	//cout<<"[1.14.0.0]"<<endl;
	//SetMassConstraint(fP,fC,mJ1,fMassHypo);
      } else if((mMassParticle < SumDaughterMass) || (fP[6]<0) ){
	//cout<<"[1.13.0.0   (mMassParticle < SumDaughterMass) || (fP[6]<0)]"<<endl;
	//SetMassConstraint(fP,fC,mJ1,SumDaughterMass);
	if (!SetMassConstraintPrivate(fP,fC,mJ1,SumDaughterMass)) return false;
	//cout<<"[1.14.0.0]"<<endl;
      }
      
      //cout<<"[1.15.0.0]"<<endl;
      
      float mDJ[7][7];

      for(Int_t i=0; i<7; i++) {
	for(Int_t j=0; j<7; j++) {
	  mDJ[i][j] = 0;
	  for(Int_t k=0; k<7; k++) {
	    mDJ[i][j] += mDf[i][k]*mJ1[j][k];
	  }
	}
      }
      //cout<<"[1.16.0.0]"<<endl;
      for(Int_t i=0; i<7; ++i){
	for(Int_t j=0; j<7; ++j){
	  mDf[i][j]=0;
	  for(Int_t l=0; l<7; l++){
	    mDf[i][j] += mJ2[i][l]*mDJ[l][j];
	  }
	}
      }
      //cout<<"[1.17.0.0]"<<endl;
      
      //* Add the daughter momentum to the particle momentum

      fP[ 3] += m[ 3];
      fP[ 4] += m[ 4];
      fP[ 5] += m[ 5];
      fP[ 6] += m[ 6];

      //cout<<"[1.18.0.0]"<<endl;
      cout<<"New  fP[0]="<<fP[0]<<"  fP[1]="<<fP[1]<<"  fP[2]="<<fP[2]<<"  fP[3]="<<fP[3]<<
	"  fP[4]="<<fP[4]<<"  fP[5]="<<fP[5]<<"  fP[6]="<<fP[6]<<"  fP[7]="<<fP[7]<<endl;

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
      //cout<<"[1.19.0.0]"<<endl;
      fC[6 ] += mDf[3][0]; fC[7 ] += mDf[3][1]; fC[8 ] += mDf[3][2];
      fC[10] += mDf[4][0]; fC[11] += mDf[4][1]; fC[12] += mDf[4][2];
      fC[15] += mDf[5][0]; fC[16] += mDf[5][1]; fC[17] += mDf[5][2];
      fC[21] += mDf[6][0]; fC[22] += mDf[6][1]; fC[23] += mDf[6][2];
      //cout<<"[1.20.0.0]"<<endl;
      fC[9 ] += mDf[3][3] + mDf[3][3];
      fC[13] += mDf[4][3] + mDf[3][4]; fC[14] += mDf[4][4] + mDf[4][4];
      fC[18] += mDf[5][3] + mDf[3][5]; fC[19] += mDf[5][4] + mDf[4][5]; fC[20] += mDf[5][5] + mDf[5][5];
      fC[24] += mDf[6][3] + mDf[3][6]; fC[25] += mDf[6][4] + mDf[4][6]; fC[26] += mDf[6][5] + mDf[5][6]; fC[27] += mDf[6][6] + mDf[6][6];
      //cout<<"[1.21.0.0]"<<endl;
      
      float K2[3][3];
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<3; j++)
	    K2[i][j] = -K[j][i];
	  K2[i][i] += 1;
	}
      //cout<<"[1.22.0.0]"<<endl;
      float A[3][3];
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  {
	    A[i][j] = 0;
	    for(int k=0; k<3; k++)
	      {
		A[i][j] += D[i][k] * K2[k][j];
	      }
	  }
      //cout<<"[1.23.0.0]"<<endl;
      double M[3][3];
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  {
	    M[i][j] = 0;
	    for(int k=0; k<3; k++)
	      {
		M[i][j] += K[i][k] * A[k][j];
	      }
	  }
      //cout<<"[1.24.0.0]"<<endl;
      fC[0] += 2*M[0][0];
      fC[1] += M[0][1] + M[1][0];
      fC[2] += 2*M[1][1];
      fC[3] += M[0][2] + M[2][0];
      fC[4] += M[1][2] + M[2][1];
      fC[5] += 2*M[2][2];
      //cout<<"[1.25.0.0]"<<endl;
      //* Calculate Chi^2 
      cout<<"New Covariance matrix"<<endl;
      for (int index=0; index<5; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=5; index<10; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=10; index<15; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=16; index<20; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=20; index<25; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;
      for (int index=25; index<28; ++index) cout<<Form("  fC[%d]=",index)<<fC[index];
      cout<<endl;

      cout<<"  mS[0]="<<mS[0]<<"  mS[1]="<<mS[1]<<"  mS[2]="<<mS[2]<<"  mS[3]="<<mS[3]<<"  mS[4]="<<mS[4]<<"  mS[5]="<<mS[5]<<endl;
      cout<<" zeta[0]="<<zeta[0]<<"  zeta[1]="<<zeta[1]<<"  zeta[2]="<<zeta[2]<<endl;

      
      fNDF  += 2;
      fQ    +=  Daughter.GetQ();
      fSFromDecay = 0;


      fChi2 += (mS[0]*zeta[0] + mS[1]*zeta[1] + mS[3]*zeta[2])*zeta[0]
	+      (mS[1]*zeta[0] + mS[2]*zeta[1] + mS[4]*zeta[2])*zeta[1]
	+      (mS[3]*zeta[0] + mS[4]*zeta[1] + mS[5]*zeta[2])*zeta[2];
      
      cout<<"fChi2="<<fChi2<<endl;
      
      if (fC[0]>100 || fC[1]>100 || fC[2]>100 || fC[3]>100 || fC[4]>100 || fC[5]>100) {
	cout<<"Don't pass fCs"<<endl;
	return false;
      }


      //cout<<"[1.26.0.0]"<<endl;
      return true;      
    }
  }
};

class K0sContainer : public TObject{ 
 public:
  
  K0sContainer(double m, double y, double mom[3], double svtx[3], int q, double dlxy, double dlxyz, double d0xy, double d0xyz, double pointxy, double pointxyz, 
	       double dcadauxy, double dcadauxyz, double chi2, int did[2], int label[2], double dmom1[3], double dmom2[3], double dd0xy[2], double dd0xyz[2]){
    Mass   = m;
    Rap    = y;
    P[0]   = mom[0];
    P[1]   = mom[1];
    P[2]   = mom[2];
    Vtx[0] = svtx[0];
    Vtx[1] = svtx[1];
    Vtx[2] = svtx[2];
    Charge = q;
    DcaXY  = d0xy;
    DcaXYZ = d0xyz;
    DecayLengthXY  = dlxy;
    DecayLengthXYZ = dlxyz;
    CpvXY  = pointxy;
    CpvXYZ = pointxyz;
    DcaDaughtersXY  = dcadauxy;
    DcaDaughtersXYZ = dcadauxyz;
    rChi2V0     = chi2;
    dP1[0]      = dmom1[0];
    dP1[1]      = dmom1[1];
    dP1[2]      = dmom1[2];
    dP2[0]      = dmom2[0];
    dP2[1]      = dmom2[1];
    dP2[2]      = dmom2[2];
    dDcaXY[0]   = dd0xy[0];
    dDcaXY[1]   = dd0xy[1];
    dDcaXYZ[0]  = dd0xyz[0];
    dDcaXYZ[1]  = dd0xyz[1];
    dLabel[0]   = label[0];
    dLabel[1]   = label[1];
    dTrackid[0] = did[0];
    dTrackid[1] = did[1];
  }
  
  int Charge;
  int dLabel[2];
  int dTrackid[2];

  double Mass;
  double Rap;
  double P[3];  
  double Vtx[3];
  double DcaXY;
  double DcaXYZ;
  double DecayLengthXY;
  double DecayLengthXYZ;
  double CpvXY;
  double CpvXYZ;
  double DcaDaughtersXY;
  double DcaDaughtersXYZ;
  double rChi2V0;    
  double dP1[3];
  double dP2[3];
  double dDcaXY[2];
  double dDcaXYZ[2];
  
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

  void setK0sCuts(bool decaylength, bool pointangle, bool chi2, bool paiDCA, bool trackdistance, bool dcaxy, bool lifetime){
    onDecayLengthXYCut           = decaylength;
    onPointingAngleXYCut         = pointangle;
    onChi2Cut                    = chi2;
    onDaughterPairDCAXYCut       = paiDCA;
    onDaughterTrackDistanceXYCut = trackdistance;
    onDCAXYCut                   = dcaxy;
    onLifetimeCut                = lifetime;
  }
  
  void setK0sCus(TFile *inFile){        
    fHistDecayLengthXYCut           = (TH1F*)inFile -> Get("fHistDecayLengthXYCut");
    fHistPointingAngleXYCut         = (TH1F*)inFile -> Get("fHistPointingAngleXYCut");
    fHistChi2Cut                    = (TH1F*)inFile -> Get("fHistChi2Cut");
    fHistDaughterPairDCAXYCut       = (TH1F*)inFile -> Get("fHistDaughterPairDCAXYCut");
    fHistDaughterTrackDistanceXYCut = (TH1F*)inFile -> Get("fHistDaughterTrackDistanceXYCut");
    fHistDCAxyCut                   = (TH1F*)inFile -> Get("fHistDCAxyCut");
    fHistLifeTimeCut                = (TH1F*)inFile -> Get("fHistLifeTimeCut");
  }

private:
  AliAnalysisTaskAODTrackPair(
      const AliAnalysisTaskAODTrackPair &); // not implemented
  AliAnalysisTaskAODTrackPair &
  operator=(const AliAnalysisTaskAODTrackPair &); // not implemented

  bool Initialize();
      
  bool TrackPIDChecker(AliAODTrack *track, AliPID::EParticleType pid, bool isSel);  
  bool TrackQualityChecker(AliAODTrack *track);
  bool V0QualityChecker(K0sContainer *v0, bool isSel);
  
  bool K0sPairAnalysis(AliPID::EParticleType pid1, AliPID::EParticleType pid2);
  bool K0sPairAnalysisEventMixing(AliPID::EParticleType pid1,AliPID::EParticleType pid2);

  bool AddK0sArray(AliPID::EParticleType pid1, AliPID::EParticleType pid2);

  K0sContainer* calcK0sFromTracks(AliAODTrack* aodTrack1,AliAODTrack* aodTrack2);
  
  bool updateAODv0(AliAODv0* v0);

  int findLeadingTrack();
  
  KFParticle CreateKFParticle(AliAODTrack* track,int pdg);
  
  bool CheckTrackCovariance(AliAODTrack* track);
  bool CheckTrackCovariance(KFParticle track);

  double DecayLengthFromKF(KFParticle kfpParticle, KFParticle PV);
  double DecayLengthXYFromKF(KFParticle kfpParticle, KFParticle PV);

  bool ProcessMC();

  bool  isAcceptK0sPair(K0sContainer* v0_1, K0sContainer* v0_2);
  bool  isAcceptV0QualityCuts(K0sContainer* v0, bool* isAcceptCuts);
  bool  isAcceptK0sCuts(K0sContainer* v0);
  bool  isAcceptK0sKinematics(K0sContainer* v0);
  bool  isAcceptK0sKinematics(AliAODMCParticle* v0);
  bool  isAcceptPrimaryTrack(AliAODTrack* track);

  AliAODEvent *fEvent;
  AliAODVertex *fPrimVtx;
  
  AliEventPoolManager *fPoolMuonTrackMgrK0s;  
  AliEventPoolManager *fPoolMuonTrackMgrPion;  
  
  AliMultSelection *fMultSelection;
  

  AliAnalysisTaskAODTrackPairUtils *fUtils;
  
  TClonesArray *fMCTrackArray;
  TObjArray* fArrayK0s;

  TH1F* fHistDecayLengthXYCut;
  TH1F* fHistPointingAngleXYCut;
  TH1F* fHistChi2Cut;
  TH1F* fHistDaughterPairDCAXYCut;
  TH1F* fHistDaughterTrackDistanceXYCut;
  TH1F* fHistDCAxyCut;
  TH1F* fHistLifeTimeCut;
  
  int fTrackDepth;
  int fPoolSize;

  double fReadyFraction;
  
  bool onVerbose;

  bool onEvtMixingPoolVtxZ;
  bool onEvtMixingPoolCent;
  bool onEvtMixingPoolPsi;

  bool onDecayLengthXYCut;
  bool onPointingAngleXYCut;
  bool onChi2Cut;
  bool onDaughterPairDCAXYCut;
  bool onDaughterTrackDistanceXYCut;
  bool onDCAXYCut;
  bool onLifetimeCut;
  

  bool fIsMC;
  bool fIsMixingAnalysis;
  bool fIsManualV0Analysis;
  
  double fMinK0sRap;
  double fMaxK0sRap;
  double fMinK0sPt;
  double fMaxK0sPt;

  std::string fMethodCent;
  
  double fPrimVtxPos[3];
  double fCent;
  double fPsi;

  int fMinNContPrimVtx;

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
  
  TH1F *fHistTPCNClusts;
  TH1F *fHistSPDNClusts;
  TH1F *fHistTPCCrossRowsFindableRatio;
  TH1F *fHistReducedChi2TPC;
  TH1F *fHistReducedChi2ITS;
  TH1F *fHistDCAz;
  TH2F *fHistDCAxyPt;
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




  int all;
  int abord;
  
  ClassDef(AliAnalysisTaskAODTrackPair, 1); // example of analysis
};

#endif
