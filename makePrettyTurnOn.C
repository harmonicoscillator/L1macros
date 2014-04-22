#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>

void makePrettyTurnOn()
{
  TFile *inFile = TFile::Open("hist_out_akPu3Calo.root");
  TGraphAsymmErrors *asymm_pt_0 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_0");
  TGraphAsymmErrors *asymm_pt_15 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_15");
  TGraphAsymmErrors *asymm_pt_30 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_30");
  TGraphAsymmErrors *asymm_pt_60 = (TGraphAsymmErrors*)inFile->Get("asymm_pt_60");

  asymm_pt_0->SetMarkerColor(1);
  asymm_pt_15->SetMarkerColor(kBlue);
  asymm_pt_30->SetMarkerColor(kRed);
  asymm_pt_60->SetMarkerColor(kOrange);

  asymm_pt_0->SetLineColor(1);
  asymm_pt_15->SetLineColor(kBlue);
  asymm_pt_30->SetLineColor(kRed);
  asymm_pt_60->SetLineColor(kOrange);

  const int nBins = 75;
  const double maxPt = 150;

  TH1D *hEmpty = new TH1D("hEmpty",";akPu3Calo Jet p_{T} (GeV);Efficiency",nBins,10,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  asymm_pt_0->Draw("p");
  asymm_pt_15->Draw("p");
  asymm_pt_30->Draw("p");
  asymm_pt_60->Draw("p");

  TLegend *leg = new TLegend(0.7,0.2,0.9,0.4,"L1 Threshold (GeV)");
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  //leg->SetTextSize(20);

  leg->AddEntry(asymm_pt_0,"0","lp");
  leg->AddEntry(asymm_pt_15,"15","lp");
  leg->AddEntry(asymm_pt_30,"30","lp");
  leg->AddEntry(asymm_pt_60,"60","lp");

  leg->Draw();

  c1->SaveAs("minbiasHI_turnon_hisub.pdf");

}

int main()
{
  makePrettyTurnOn();
  return 0;
}
