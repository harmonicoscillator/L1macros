#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

void makeRateCurve()
{
  const char *type = "photon_data";
  //TFile *inFile = TFile::Open(Form("hist_out_%s_cleaned.root",type));
  const int numfile = 2;
  TFile *inFile[numfile];
  inFile[0] = TFile::Open(Form("hist_%s.root",type));
  inFile[1] = TFile::Open("hist_photon_iso_data.root");
  // inFile[2] = TFile::Open("hist_minbias_forward_Pu.root");

  TH1D *counts[numfile];
  for(int i = 0; i < numfile; i++)
    counts[i] = (TH1D*)inFile[i]->Get("l1Pt");

  const int nBins = 100;
  const double maxPt = 100;

  TH1D *rate[numfile];
  rate[0] = new TH1D("rate",";L1 p_{T};Rate (w.r.t. 2011)",nBins,0,maxPt);
  rate[1] = (TH1D*)rate[0]->Clone("1");
  // rate[2] = (TH1D*)rate[0]->Clone("2");

  // double total_integral[2];
  // for(int i = 0; i < 2; i++)
  //   total_integral[i] = counts[i]->Integral();
  double total_integral = counts[0]->Integral();

  for(int n = 0; n < numfile; n++)
  {
    for(int i = 0; i < nBins; i++)
    {
      double j = (double)i*(double)maxPt/(double)nBins;
      //std::cout << j << std::endl;
      double integral = counts[n]->Integral(i+1, nBins);
      //std::cout << integral << std::endl;
      rate[n]->Fill(j, (double)integral/total_integral);
      //rate[n]->Fill(j,integral);
      std::cout << "Threshold: " << j << " Rate: " << (double)integral/total_integral << std::endl;
    }
  }

  //rate[0]->SetMaximum(2);
  //rate[0]->SetMinimum(0.0001);
  TCanvas *c1 = new TCanvas();
  rate[0]->Draw("l");
  c1->SetLogy();

  rate[1]->SetLineColor(kRed);
  rate[1]->Draw("l same");

  // rate[2]->SetLineColor(kBlue);
  // rate[2]->Draw("l same");

  TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  leg->SetTextSize(20);

  leg->AddEntry(rate[0],"inclusive", "l");
  leg->AddEntry(rate[1],"isolated", "l");
  // leg->AddEntry(rate[2],"forward", "l");
  leg->Draw();

  c1->SaveAs(Form("%s_rate.pdf",type));


  // for(int i = 1; i < rate[0]->GetNbinsX(); ++i)
  // {
  //   std::cout << "Threshold: " << (i-1)*300/75 << "GeV "
  // 	      << "Rate: " << rate[0]->GetBinContent(i) << std::endl;
  // }
}
