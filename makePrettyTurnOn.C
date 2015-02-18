#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>

void makePrettyTurnOn()
{
  const char *type = "minbias_photon_data";
  TFile *inFile = TFile::Open(Form("hist_%s.root",type));

  const Int_t THRESHOLDS = 4;
  const Double_t L1_THRESHOLD[THRESHOLDS] = {2, 5, 8, 12};
  const Int_t COLORS[THRESHOLDS] = {kBlack, kRed, kBlue, kGreen+3};//, kMagenta+3};
  const char* LABELS[2] = {"central", "periph"};
  TGraphAsymmErrors *asymm[THRESHOLDS][2];

  for(int i = 0; i < THRESHOLDS; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      asymm[i][j] = (TGraphAsymmErrors*)inFile->Get(Form("asymm_pt_%d_%d",(int)L1_THRESHOLD[i],j+1));
      asymm[i][j]->SetMarkerColor(COLORS[i]);
      asymm[i][j]->SetLineColor(COLORS[i]);
    }
    asymm[i][1]->SetMarkerStyle(25);
  }

  const int nBins = 100;
  const double maxPt = 100;

  TH1D *hEmpty = new TH1D("hEmpty",Form(";Offline Photon p_{T} (GeV);Efficiency"),nBins,0,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  TLine *line = new TLine(0,1,maxPt,1);
  line->Draw();

  for(int i = 0; i < THRESHOLDS; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      asymm[i][j]->Draw("p");
    }
  }

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.5);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(18);

  for(int i = 0; i < THRESHOLDS; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      leg->AddEntry(asymm[i][j], Form("Upgrade RCT %d, %s", (int)L1_THRESHOLD[i], LABELS[j]), "lp");
    }
  }

  leg->Draw();

  c1->SaveAs(Form("%s_turnon.pdf",type));

}

int main()
{
  makePrettyTurnOn();
  return 0;
}
