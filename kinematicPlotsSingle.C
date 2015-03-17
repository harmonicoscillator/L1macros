#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

void kinematicPlotsSingle()
{
  TH1::SetDefaultSumw2();

  //TString l1_input = "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root";
  //TString l1_input = "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root";

  TString l1_input = "/mnt/hadoop/cms/store/user/men1/L1Data/HIL1DPG/MinBias/HIMinBiasUPC_Skim_HLT_HIMinBiasHfOrBSC_v2_CaloRegionEta516_CMSSW740pre7/L1NTupleMBHIFS.root";

  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_event, l1_run, l1_lumi;
  Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);

  const int NHISTS = 22;
  TH1D *ebeSigma[NHISTS];
  TH1D *range[NHISTS];
  range[0] = new TH1D("range0","range;max pT - min pT",1024,0,1024);
  TH1D *phiDelta[NHISTS];
  phiDelta[0] = new TH1D("phiDelta0","phiDelta;PhiMax - PhiMin",10,0,10);
  ebeSigma[0] = new TH1D("ebeSigma0","ebeSigma;#sigma;count",300,0,50);
  TH1D *distSigma = new TH1D("distSigma",";#phi index;#sigma", NHISTS,0,NHISTS);
  TH1D *distPt = new TH1D("distPt",";#phi index;<p_{T}>",NHISTS,0,NHISTS);

  TH1D *ptDists[NHISTS];
  ptDists[0] = new TH1D("ptDists0",";region_hwPt",1024,0,1024);
  for(int i = 1; i < NHISTS; i++)
  {
    ebeSigma[i] = (TH1D*)ebeSigma[0]->Clone(Form("ebeSigma%i",i));
    ptDists[i] = (TH1D*)ptDists[0]->Clone(Form("ptDists%i",i));
    range[i] = (TH1D*)range[0]->Clone(Form("range%i",i));
    phiDelta[i] = (TH1D*)phiDelta[0]->Clone(Form("phiDelta%i",i));
  }

  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);

    double sums[NHISTS];
    double sums2[NHISTS];
    double maxLocation[NHISTS];
    double maxValue[NHISTS];
    double minLocation[NHISTS];
    double minValue[NHISTS];
    for(int i = 0; i < NHISTS; i++)
    {
      sums[i] = 0.;
      sums2[i] = 0.;
      maxLocation[i] = -1;
      maxValue[i] = -1;
      minLocation[i] = -1;
      minValue[i] = -1;
    }
    for(int i = 0; i < 396; i++)
    {
      //if(region_hwEta[i] != 4) continue;
      ptDists[region_hwEta[i]]->Fill(region_hwPt[i]);
      sums[region_hwEta[i]] += region_hwPt[i];
      sums2[region_hwEta[i]] += (region_hwPt[i] * region_hwPt[i]);

      if(maxValue[region_hwEta[i]] < region_hwPt[i])
      {
	maxValue[region_hwEta[i]] = region_hwPt[i];
	maxLocation[region_hwEta[i]] = region_hwPhi[i];
      }
      if(minValue[region_hwEta[i]] > region_hwPt[i])
      {
	minValue[region_hwEta[i]] = region_hwPt[i];
	minLocation[region_hwEta[i]] = region_hwPhi[i];
      }
    }

    for(int i = 0; i < NHISTS; i++)
    {
      double ebesigma = TMath::Sqrt( (sums2[i]) - ((sums[i])*(sums[i])) );
      ebeSigma[i]->Fill(ebesigma);
      //std::cout << ebesigma << std::endl;
      range[i]->Fill(maxValue[i] - minValue[i]);
      double phidelta = TMath::Abs(maxLocation[i] - minLocation[i]);
      if(phidelta > 9)
	phidelta = 18 - phidelta;
      phiDelta[i]->Fill(phidelta);
    }
  }
  for(int i = 0; i < NHISTS; i++)
  {
    distSigma->SetBinContent(i+1, ptDists[i]->GetStdDev());
    distSigma->SetBinError(i+1, ptDists[i]->GetStdDevError());

    distPt->SetBinContent(i+1, ptDists[i]->GetMean());
    distPt->SetBinError(i+1, ptDists[i]->GetMeanError());
  }

  TFile *outFile = TFile::Open("kinematicDists_MBData.root","RECREATE");
  outFile->cd();
  for(int i = 0; i < NHISTS; i++)
  {
    ptDists[i]->Write();
    ebeSigma[i]->Write();
    range[i]->Write();
    phiDelta[i]->Write();
  }
  distSigma->Write();
  distPt->Write();
  outFile->Close();
}

int main()
{
  kinematicPlotsSingle();
  return 0;
}
