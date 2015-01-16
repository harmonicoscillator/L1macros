#ifndef CalculateIsolation
#define CalculateIsolation

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TMath.h>

#include <vector>
#include <iostream>

class Isolation{
public:
  Isolation(){};
  ~Isolation(){};
  void FillRegions(int pt[396], int eta[396], int phi[396]);
  int iso3x3(int eta, int phi);
  int isoCross(int eta, int phi);
  int isoFourPoint(int eta, int phi);
  int isoSingle(int eta, int phi);
  int isoHolePunch(int eta, int phi);

private:
  int region_hwPt_[396];
  int region_hwEta_[396];
  int region_hwPhi_[396];

  int region_subPt_[396];
  
};

#endif

void Isolation::FillRegions(int pt[396], int eta[396], int phi[396])
{
  // fill local copy
  for(int i = 0; i < 396; ++i)
  {
    region_hwPt_[i] = pt[i];
    region_hwEta_[i] = eta[i];
    region_hwPhi_[i] = phi[i];
  }

  // compute bkg-subtracted pt for each region
  int puLevelHI[22];
  double r_puLevelHI[22];
  for(unsigned i = 0; i < 22; ++i)
  {
    puLevelHI[i] = 0;
    r_puLevelHI[i] = 0.0;
  }

  for(int i = 0; i < 396; ++i){
    r_puLevelHI[region_hwEta_[i]] += region_hwPt_[i];
  }

  for(unsigned i = 0; i < 22; ++i)
  {
    puLevelHI[i] = floor(r_puLevelHI[i]/18 + 0.5);
  }

  for(int i = 0; i < 396; ++i)
  {
    region_subPt_[i] = std::max(0, region_hwPt_[i] - puLevelHI[region_hwEta_[i]]);
  }
}

int Isolation::iso3x3(int eta, int phi)
{
  int iso = 0;
  for(int i = 0; i < 396; ++i){
    int diffeta = eta - region_hwEta_[i];
    int diffphi = phi - region_hwPhi_[i];

    if(TMath::Abs(diffphi) == 17)
      diffphi = 1;

    if((TMath::Abs(diffeta) <= 1) && (TMath::Abs(diffphi) <= 1))
      iso += region_subPt_[i];
  }

  return iso;
}

int Isolation::isoCross(int eta, int phi)
{
  int iso = 0;
  for(int i = 0; i < 396; ++i){
    int diffeta = eta - region_hwEta_[i];
    int diffphi = phi - region_hwPhi_[i];

    if(TMath::Abs(diffphi) == 17)
      diffphi = 1;

    if((TMath::Abs(diffeta) <= 1) && (TMath::Abs(diffphi) == 0)) {
      iso += region_subPt_[i];
    } else if ( (TMath::Abs(diffphi) == 1) && (TMath::Abs(diffeta) == 0)) {
      iso += region_subPt_[i];
    }
  }

  return iso;
}

int Isolation::isoFourPoint(int eta, int phi)
{
  int iso = 0;
  for(int i = 0; i < 396; ++i){
    int diffeta = eta - region_hwEta_[i];
    int diffphi = phi - region_hwPhi_[i];

    if(TMath::Abs(diffphi) == 17)
      diffphi = 1;
    
    if((TMath::Abs(diffeta) == 1) && (TMath::Abs(diffphi) == 1))
      iso += region_subPt_[i];
  }

  return iso;
}

int Isolation::isoSingle(int eta, int phi)
{
  int iso = 0;
  for(int i = 0; i < 396; ++i){
    int diffeta = eta - region_hwEta_[i];
    int diffphi = phi - region_hwPhi_[i];

    if(TMath::Abs(diffphi) == 17)
      diffphi = 1;
    
    if((TMath::Abs(diffeta) == 0) && (TMath::Abs(diffphi) == 0))
      iso += region_subPt_[i];
  }

  return iso;
}

int Isolation::isoHolePunch(int eta, int phi)
{
  int iso = 0;
  for(int i = 0; i < 396; ++i){
    int diffeta = eta - region_hwEta_[i];
    int diffphi = phi - region_hwPhi_[i];

    if(TMath::Abs(diffphi) == 17)
      diffphi = 1;

    if((TMath::Abs(diffeta) == 1) && (TMath::Abs(diffphi) == 0)) {
      iso += region_subPt_[i];
    } else if ( (TMath::Abs(diffphi) == 1) && (TMath::Abs(diffeta) == 0)) {
      iso += region_subPt_[i];
    }
  }

  return iso;
}
