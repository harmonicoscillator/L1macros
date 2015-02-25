#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

void makeRateCurve_fromTuple()
{
  TH1::SetDefaultSumw2();

  const int MAXJETS = 10;
  const int nBins = 75;
  const double maxPt = 300;

  //const char *type = "hydjet_jets_gen";
  TFile *inFile = TFile::Open("/export/d00/scratch/luck/minbias_hydjet_l1ntuple_v2.root");
  TTree *inTree = (TTree*)inFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_num;
  inTree->SetBranchAddress("nJet",&l1_num);

  Double_t l1_pt[MAXJETS];
  inTree->SetBranchAddress("jet_pt",l1_pt);

  TH1D *counts = new TH1D("counts","counts",nBins,0,maxPt);

  long long entries = inTree->GetEntries();
  for(long long i = 0; i < entries; ++i)
  {
    inTree->GetEntry(i);

    double maxl1pt = -1;
    for(int j = 0; j < l1_num; ++j)
    {
      if(l1_pt[j] > maxl1pt)
	maxl1pt = l1_pt[j];
    }

    counts->Fill(maxl1pt);
  }

  TCanvas *c0 = new TCanvas();
  counts->Draw();

  TH1D *rate;
  rate = new TH1D("rate",";L1 p_{T};Rate (w.r.t. 5.02TeV Hydjet)",nBins,0,maxPt);
  double total_integral = counts->Integral();

  for(int i = 0; i < nBins; i++)
  {
    double j = (double)i*(double)maxPt/(double)nBins;
    //std::cout << j << std::endl;
    double integral = counts->Integral(i+1, nBins);
    //std::cout << integral << std::endl;
    rate->Fill(j, (double)integral/total_integral);
    //rate[n]->Fill(j,integral);
    std::cout << "Threshold: " << j << " Rate: " << (double)integral/total_integral << std::endl;
  }

  TCanvas *c1 = new TCanvas();
  rate->Draw("l");
  c1->SetLogy();

  TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  leg->SetTextSize(20);

  leg->AddEntry(rate,"Inclusive", "l");
  //leg->AddEntry(rate[1],"Isolated", "l");
  // leg->AddEntry(rate[2],"forward", "l");
  //leg->Draw();

  c1->SaveAs("hydjet_jets_rate.pdf");
}
