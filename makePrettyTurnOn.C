#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>

void makePrettyTurnOn()
{
  const char *type = "dijet15_forward_gen";
  TFile *inFile = TFile::Open(Form("hist_%s.root",type));
  //TFile *inFile = TFile::Open("hist_hydjet_gen.root");
  //TGraphAsymmErrors *asymm_pt_0_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_0_0");
  //TGraphAsymmErrors *asymm_pt_15_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_15_0");
  //TGraphAsymmErrors *asymm_pt_30_cen = (TGraphAsymmErrors*)inFile->Get("asymm_pt_30_0");
  TGraphAsymmErrors *asymm_pt_12_0 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_12_0");
  TGraphAsymmErrors *asymm_pt_24_0 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_24_0");
  TGraphAsymmErrors *asymm_pt_40_0 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_40_0");
  //TGraphAsymmErrors *asymm_pt_0_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_0_1");
  //TGraphAsymmErrors *asymm_pt_15_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_15_1");
  //TGraphAsymmErrors *asymm_pt_30_periph = (TGraphAsymmErrors*)inFile->Get("asymm_pt_30_1");
  TGraphAsymmErrors *asymm_pt_12_1 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_12_1");
  TGraphAsymmErrors *asymm_pt_24_1 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_24_1");
  TGraphAsymmErrors *asymm_pt_40_1 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_40_1");

  // asymm_pt_0_cen->SetMarkerColor(1);
  // asymm_pt_0_cen->SetLineColor(1);
  // asymm_pt_15_cen->SetMarkerColor(kBlue);
  // asymm_pt_15_cen->SetLineColor(kBlue);
  // asymm_pt_30_cen->SetMarkerColor(kRed);
  // asymm_pt_30_cen->SetLineColor(kRed);
  //asymm_pt_60_cen->SetMarkerColor(kOrange);
  //asymm_pt_60_cen->SetLineColor(kOrange);

  // asymm_pt_0_periph->SetMarkerColor(1);
  // asymm_pt_0_periph->SetLineColor(1);

  // asymm_pt_15_periph->SetMarkerColor(kBlue);
  // asymm_pt_15_periph->SetLineColor(kBlue);
  // asymm_pt_30_periph->SetMarkerColor(kRed);
  // asymm_pt_30_periph->SetLineColor(kRed);
  //asymm_pt_60_periph->SetMarkerColor(kOrange);
  //asymm_pt_60_periph->SetLineColor(kOrange);

  // asymm_pt_0_periph->SetMarkerStyle(24);
  // asymm_pt_15_periph->SetMarkerStyle(24);
  // asymm_pt_30_periph->SetMarkerStyle(24);

  asymm_pt_12_1->SetMarkerStyle(25);
  asymm_pt_24_1->SetMarkerStyle(25);
  asymm_pt_40_1->SetMarkerStyle(25);

  asymm_pt_24_0->SetMarkerColor(kRed);
  asymm_pt_24_1->SetMarkerColor(kRed);
  asymm_pt_24_0->SetLineColor(kRed);
  asymm_pt_24_1->SetLineColor(kRed);

  asymm_pt_40_0->SetMarkerColor(kBlue);
  asymm_pt_40_1->SetMarkerColor(kBlue);
  asymm_pt_40_0->SetLineColor(kBlue);
  asymm_pt_40_1->SetLineColor(kBlue);

  const int nBins = 100;
  const double maxPt = 250;

  TH1D *hEmpty = new TH1D("hEmpty",Form(";Gen-Jet p_{T} (GeV);Efficiency"),nBins,0,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  TLine *line = new TLine(0,1,maxPt,1);
  line->Draw();

  // asymm_pt_0_cen->Draw("p");
  // asymm_pt_15_cen->Draw("p");
  // asymm_pt_30_cen->Draw("p");
  // asymm_pt_60_cen->Draw("p");
  // asymm_pt_0_periph->Draw("p");
  // asymm_pt_15_periph->Draw("p");
  // asymm_pt_30_periph->Draw("p");
  // asymm_pt_60_periph->Draw("p");

  asymm_pt_12_0->Draw("p");
  asymm_pt_12_1->Draw("p");
  asymm_pt_24_0->Draw("p");
  asymm_pt_24_1->Draw("p");
  asymm_pt_40_0->Draw("p");
  asymm_pt_40_1->Draw("p");

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.5);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(18);

  //leg->AddEntry(asymm_pt_0_cen,"0 cen","lp");
  //leg->AddEntry(asymm_pt_15_cen,"15 cen","lp");
  //leg->AddEntry(asymm_pt_30_cen,"30 cen","lp");
  leg->AddEntry(asymm_pt_12_0,"Thresh. 12, central","lp");

  // leg->AddEntry(asymm_pt_0_periph,"0 periph","lp");
  // leg->AddEntry(asymm_pt_15_periph,"15 periph","lp");
  // leg->AddEntry(asymm_pt_30_periph,"30 periph","lp");
  leg->AddEntry(asymm_pt_12_1,"Thresh. 12, peripheral","lp");
  leg->AddEntry(asymm_pt_24_0,"Thresh. 24, central","lp");
  leg->AddEntry(asymm_pt_24_1,"Thresh. 24, peripheral","lp");
  leg->AddEntry(asymm_pt_40_0,"Thresh. 40, central","lp");
  leg->AddEntry(asymm_pt_40_1,"Thresh. 40, peripheral","lp");



  leg->Draw();

  c1->SaveAs(Form("%s_turnon.pdf",type));

}

int main()
{
  makePrettyTurnOn();
  return 0;
}
