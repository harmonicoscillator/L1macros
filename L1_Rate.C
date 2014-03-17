#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <stdio.h>
#include <TLegendEntry.h>

// #include "TROOT.h"
// #include "TSystem.h"
// #include "TStyle.h"
// #include "FWCore/FWLite/interface/AutoLibraryLoader.h"

// #include "/net/hisrv0001/home/luck/.dotfiles/rootlogon.C"

void L1_Rate()
{
  //rootlogon();
  
  const int PT_MAX = 130;
  const int nPTBINS = PT_MAX/2; // LSB is 0.5

  TChain *c = new TChain("Events","chain");
  c->Add("/mnt/hadoop/cms/store/user/luck/L1Emulator/newLayer2_v2/L1Emulator_newLayer2_*.root");

  TH1D *pt = new TH1D("pt",";L1 Jet p_{T} (GeV);Count",nPTBINS,0,PT_MAX);
  // use [0][0] to get only the highest-pt object
  c->Project(pt->GetName(), "l1tJetBXVector_Layer2Phys__L1TEMULATION.obj.data_.l1t::L1Candidate.pt()[0][0]");

  TLegend *leg = new TLegend(0.65, 0.55, 0.9, 0.7);
  leg->SetFillColor(0);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->AddEntry(pt,"Layer2 HI Sub","l");
  TLatex *tex = new TLatex(0.45,0.9,"HI MinBias Data");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(20);
  
  TCanvas *c1 = new TCanvas();
  pt->Draw();
  leg->Draw();
  tex->Draw();
  c1->SaveAs("layer2_minbias_pt.pdf");

  TH1D *rate = new TH1D("rate",";L1 Jet p_{T} (GeV);Rate (Arb.)",
  			      nPTBINS,0,PT_MAX);
  double total_integral = pt->Integral();
  for(int i = 0; i < nPTBINS; i++)
  {
    double j = (double)i*(double)PT_MAX/(double)nPTBINS;
    double integral = pt->Integral(i, nPTBINS);
    rate->Fill(j, (double)integral/total_integral);
  }

  TCanvas *c2 = new TCanvas();
  rate->Draw();
  leg->Draw();
  tex->Draw();

  const TString trigNames[] = {"L1_SingleJet16",
			       "L1_SingleJet36_BptxAND",
			       "L1_SingleJet52_BptxAND",
			       "L1_SingleJet68_BptxAND",
			       "L1_SingleJet80_BptxAND",
			       "L1_SingleJet92",
			       "L1_SingleJet128"};
  const Double_t trigPt[] = {16,
			     36,
			     52,
			     68,
			     80,
			     92,
			     128};

  const int nTrigs = 7;
  
  Double_t newTrigRate[nTrigs];
  TChain *newhlt = new TChain("hltanalysis/HltTree","newhlt");
  newhlt->Add("~/scratch/openHLT_newLayer2.root");
  for(int i = 0; i < nTrigs; ++i)
  {
    int total = newhlt->GetEntries();
    TString newString = trigNames[i]+"==1";
    const char *selection = newString.Data();
    int pass = newhlt->GetEntries(selection);
    newTrigRate[i] = (double)pass/(double)total;
  }    
  TGraph *newhlt_rate = new TGraph(nTrigs,trigPt,newTrigRate);
  newhlt_rate->SetMarkerColor(kBlue);
  newhlt_rate->Draw("p");

  Double_t oldTrigRate[nTrigs];
  TChain *oldhlt = new TChain("hltanalysis/HltTree","oldhlt");
  oldhlt->Add("~/scratch/openHLT_oldGCT.root");
  for(int i = 0; i < nTrigs; ++i)
  {
    int total = oldhlt->GetEntries();
    TString newString = trigNames[i]+"==1";
    const char *selection = newString.Data();
    int pass = oldhlt->GetEntries(selection);
    oldTrigRate[i] = (double)pass/(double)total;
  }    
  TGraph *oldhlt_rate = new TGraph(nTrigs,trigPt,oldTrigRate);
  oldhlt_rate->SetMarkerColor(kRed);
  oldhlt_rate->Draw("p");
  TLegendEntry *le = leg->AddEntry(oldhlt_rate->GetName(),"oldGCT from HLT","p");
  le->SetMarkerColor(kRed);
  leg->AddEntry(newhlt_rate->GetName(),"Layer2 from HLT","p");

  leg->Draw();
  c2->SetLogy();
  c2->SaveAs("layer2_minbias_rate.pdf");
}
