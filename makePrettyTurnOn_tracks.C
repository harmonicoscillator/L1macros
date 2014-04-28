#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>

void makePrettyTurnOn_tracks()
{
  TFile *inFile = TFile::Open("hist_out_tracks.root");
  TGraphAsymmErrors *asymm_pt_0_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_0_cen");
  TGraphAsymmErrors *asymm_pt_15_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_15_cen");
  TGraphAsymmErrors *asymm_pt_30_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_30_cen");
  TGraphAsymmErrors *asymm_pt_60_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_60_cen");
  TGraphAsymmErrors *asymm_pt_0_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_0_periph");
  TGraphAsymmErrors *asymm_pt_15_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_15_periph");
  TGraphAsymmErrors *asymm_pt_30_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_30_periph");
  TGraphAsymmErrors *asymm_pt_60_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_60_periph");


  asymm_pt_0_cen->SetMarkerColor(1);
  asymm_pt_0_cen->SetLineColor(1);
  asymm_pt_15_cen->SetMarkerColor(kBlue);
  asymm_pt_15_cen->SetLineColor(kBlue);
  asymm_pt_30_cen->SetMarkerColor(kRed);
  asymm_pt_30_cen->SetLineColor(kRed);
  asymm_pt_60_cen->SetMarkerColor(kOrange);
  asymm_pt_60_cen->SetLineColor(kOrange);

  asymm_pt_0_periph->SetMarkerColor(1);
  asymm_pt_0_periph->SetLineColor(1);
  asymm_pt_15_periph->SetMarkerColor(kBlue);
  asymm_pt_15_periph->SetLineColor(kBlue);
  asymm_pt_30_periph->SetMarkerColor(kRed);
  asymm_pt_30_periph->SetLineColor(kRed);
  asymm_pt_60_periph->SetMarkerColor(kOrange);
  asymm_pt_60_periph->SetLineColor(kOrange);

  asymm_pt_0_periph->SetMarkerStyle(24);
  asymm_pt_15_periph->SetMarkerStyle(24);
  asymm_pt_30_periph->SetMarkerStyle(24);
  asymm_pt_60_periph->SetMarkerStyle(24);


  const int nBins = 100;
  const double maxPt = 120;

  TH1D *hEmpty = new TH1D("hEmpty",";Offline max track p_{T} (GeV);Efficiency",nBins,10,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  asymm_pt_0_cen->Draw("p");
  asymm_pt_15_cen->Draw("p");
  asymm_pt_30_cen->Draw("p");
  asymm_pt_60_cen->Draw("p");
  asymm_pt_0_periph->Draw("p");
  asymm_pt_15_periph->Draw("p");
  asymm_pt_30_periph->Draw("p");
  asymm_pt_60_periph->Draw("p");

  TLegend *leg = new TLegend(0.6,0.2,0.9,0.5,"L1 Threshold (GeV)");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(20);

  leg->AddEntry(asymm_pt_0_cen,"0 cen","lp");
  leg->AddEntry(asymm_pt_15_cen,"15 cen","lp");
  leg->AddEntry(asymm_pt_30_cen,"30 cen","lp");
  leg->AddEntry(asymm_pt_60_cen,"60 cen","lp");

  leg->AddEntry(asymm_pt_0_periph,"0 periph","lp");
  leg->AddEntry(asymm_pt_15_periph,"15 periph","lp");
  leg->AddEntry(asymm_pt_30_periph,"30 periph","lp");
  leg->AddEntry(asymm_pt_60_periph,"60 periph","lp");

  leg->Draw();

  c1->SaveAs("minbiasHI_trkturnon_hisub.pdf");

}

int main()
{
  makePrettyTurnOn_tracks();
  return 0;
}
