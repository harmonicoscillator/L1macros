#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>

void makeRateCurve()
{
  const char *type = "akPu3CaloJets";
  TFile *inFile = TFile::Open(Form("hist_out_%s_cleaned.root",type));
  TH1D *counts = (TH1D*)inFile->Get("l1Pt");

  const int nBins = 75;
  const double maxPt = 300;

  TH1D *rate = new TH1D("rate",";L1 p_{T};Rate (w.r.t. PbPb 2011)",nBins,0,maxPt);

  double total_integral = counts->Integral();

  for(int i = 0; i < nBins; i++)
  {
    double j = (double)i*(double)maxPt/(double)nBins;
    //std::cout << j << std::endl;
    double integral = counts->Integral(i+1, nBins);
    //std::cout << integral << std::endl;
    rate->Fill(j, (double)integral/total_integral);
    //std::cout << (double)integral/total_integral << std::endl;
  }

  rate->SetMaximum(2);
  rate->SetMinimum(0.00001);
  TCanvas *c1 = new TCanvas();
  rate->Draw("l");
  c1->SetLogy();

  c1->SaveAs(Form("minbiasHI_rate_%s.pdf",type));
}
