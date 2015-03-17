#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

void kinematicPlots()
{
  TH1::SetDefaultSumw2();

  TString l1_input = "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root";
  //TString l1_input = "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root";

  //TString l1_input = "/mnt/hadoop/cms/store/user/men1/L1Data/HIL1DPG/MinBias/HIMinBiasUPC_Skim_HLT_HIMinBiasHfOrBSC_v2_CaloRegionEta516_CMSSW740pre7/L1NTupleMBHI.root";

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

  TH1D *ebeSigma[22];
  ebeSigma[0] = new TH1D("ebeSigma0","ebeSigma;#sigma;count",300,0,50);
  TH1D *distSigma = new TH1D("distSigma",";#eta index;#sigma", 22,0,22);
  TH1D *distPt = new TH1D("distPt",";#eta index;<p_{T}>",22,0,22);

  TH1D *ptDists[22];
  ptDists[0] = new TH1D("ptDists0",";region_hwPt",1024,0,1024);
  for(int i = 1; i < 22; i++)
  {
    ebeSigma[i] = (TH1D*)ebeSigma[0]->Clone(Form("ebeSigma%i",i));
    ptDists[i] = (TH1D*)ptDists[0]->Clone(Form("ptDists%i",i));
  }

  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);

    double sums[22];
    double sums2[22];
    for(int i = 0; i < 22; i++)
    {
      sums[i] = 0.;
      sums2[i] = 0.;
    }
    for(int i = 0; i < 396; i++)
    {
      ptDists[region_hwEta[i]]->Fill(region_hwPt[i]);
      sums[region_hwEta[i]] += region_hwPt[i];
      sums2[region_hwEta[i]] += (region_hwPt[i] * region_hwPt[i]);
    }

    for(int i = 0; i < 22; i++)
    {
      double ebesigma = TMath::Sqrt( (sums2[i]/18.) - ((sums[i]/18.)*(sums[i]/18.)) );
      ebeSigma[i]->Fill(ebesigma);
      //std::cout << ebesigma << std::endl;
    }
  }
  for(int i = 0; i < 22; i++)
  {
    distSigma->SetBinContent(i+1, ptDists[i]->GetStdDev());
    distSigma->SetBinError(i+1, ptDists[i]->GetStdDevError());

    distPt->SetBinContent(i+1, ptDists[i]->GetMean());
    distPt->SetBinError(i+1, ptDists[i]->GetMeanError());
  }

  TFile *outFile = TFile::Open("kinematicDists_502Hydjet.root","RECREATE");
  outFile->cd();
  for(int i = 0; i < 22; i++)
  {
    ptDists[i]->Write();
    ebeSigma[i]->Write();
  }
  distSigma->Write();
  distPt->Write();
  outFile->Close();
}

int main()
{
  kinematicPlots();
  return 0;
}
