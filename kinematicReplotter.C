#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

void kinematicReplotter()
{
  TH1::SetDefaultSumw2();

  TString l1_input[3];
  l1_input[0] = "kinematicDists_MBData.root";
  l1_input[1] = "kinematicDists_276Hydjet.root";
  l1_input[2] = "kinematicDists_502Hydjet.root";

  TString labels[3] = {"MB Data", "2.76 MB Hydjet", "5.02 MB Hydjet"};

  Int_t colors[3] = {kBlack, kBlue, kRed};

  TFile *files[3];
  TH1D *hists[3];
  TCanvas *c1 = new TCanvas();
  TLegend *leg = new TLegend(0.6, 0.55, 0.9, 0.85);
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  for(int i = 0; i < 3; i++)
  {
    files[i] = TFile::Open(l1_input[i]);
    hists[i] = (TH1D*)files[i]->Get("distPt");
    //hists[i]->Divide(hists[i]->GetEntries());
    hists[i]->SetLineColor(colors[i]);
    hists[i]->SetMarkerColor(colors[i]);
    //hists[i]->DrawNormalized("p e same");
    hists[i]->Draw("p e same");
    leg->AddEntry(hists[i],labels[i],"pl");
    c1->SetLogy();
  }
  leg->Draw();

  c1->SaveAs("sigmaComp_distPt.pdf");
}
