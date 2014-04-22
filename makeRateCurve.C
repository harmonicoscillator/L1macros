#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

void makeRateCurve()
{
  TFile *inFile = TFile::Open("hist_out_akPu3Calo.root");
  TH1D *counts = (TH1D*)inFile->Get("l1JetPt");

  const int nBins = 75;
  const double maxPt = 150;

  TH1D *rate = new TH1D("rate",";L1 Jet p_{T};Rate (w.r.t. PbPb 2011)",nBins,0,maxPt);

  double total_integral = counts->Integral();

  for(int i = 0; i < nBins; i++)
  {
    double j = (double)i*(double)maxPt/(double)nBins;
    double integral = counts->Integral(i, nBins);
    rate->Fill(j, (double)integral/total_integral);
  }

  TCanvas *c1 = new TCanvas();
  rate->Draw("l");
  c1->SetLogy();
}
