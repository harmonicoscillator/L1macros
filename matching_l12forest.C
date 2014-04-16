#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <stdio.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>

#include <vector>
#include <map>

using namespace std;
#include "test.h"
#include "test2.h"

const Double_t L1_THRESHOLD = 30;

void matching_l12forest()
{
  const TString l1_input = "/export/d00/scratch/luck/L1Tree_minbias_chunk1.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *lTree = (TTree*)lFile->Get("l1ExtraTreeProducer/L1ExtraTree");
  TTree *lEvtTree = (TTree*)lFile->Get("l1NtupleProducer/L1Tree");

  test *l1extra = new test(lTree);
  test2 *l1 = new test2(lEvtTree);

  const TString forest_input = "/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root";
  TFile *fFile = TFile::Open(forest_input);
  TTree *fTree = (TTree*)fFile->Get("akVs3PFJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);

  Int_t f_evt, f_run, f_lumi;
  fTree->SetBranchAddress("evt",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);

  Int_t nref;
  Float_t jtpt[500], jteta[500], jtphi[500];
  fTree->SetBranchAddress("nref",&nref);
  fTree->SetBranchAddress("jtpt",jtpt);
  fTree->SetBranchAddress("jteta",jteta);
  fTree->SetBranchAddress("jtphi",jtphi);

  TFile *outFile = new TFile("hist_out.root","RECREATE");

  map<long,int> kmap;// = new unordered_map<long,int>();
  int f_entries = fTree->GetEntries();
  for(int j = 0; j < f_entries; ++j)
  {
    fTree->GetEntry(j);
    long key = 10000000000*f_run + 10000000*f_lumi + f_evt;
    pair<long,int> p(key,j);
    kmap.insert(p);
  }

  const int nBins = 75;
  const double maxPt = 150;

  TH1D *l1JetPt = new TH1D("l1JetPt",";L1 Jet p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fJetPt = new TH1D("fJetPt",";akVs3PF Jet p_{T}",nBins,0,maxPt);
  TH1D *accepted = new TH1D("accepted",";akVsPF Jet p_{T}",nBins,0,maxPt);
  TH2D *corr = new TH2D("corr","akVs pt;l1 pt",nBins,0,maxPt,nBins,0,maxPt);

  int count = 0;

  int entries = lTree->GetEntries();
  for(int j = 0; j < entries; ++j)
  {
    l1->GetEntry(j);

    int l1_evt = l1->event;
    int l1_run = l1->run;
    int l1_lumi = l1->lumi;
    long key = 10000000000*l1_run + 10000000*l1_lumi + l1_evt;

    map<long,int>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      fTree->GetEntry(got->second);
      // if(l1_evt != f_evt)
      // {
      // 	printf("ERROR: not actually the same event");
      // 	exit(1);
      // }
      l1extra->GetEntry(j);

      //double maxl1pt = TMath::Max(l1extra->cenJetEt[0],l1extra->fwdJetEt[0]);
      //double maxl1pt = l1extra->cenJetEt[0];
      double maxl1pt = l1extra->fwdJetEt[0];

      l1JetPt->Fill(maxl1pt);
      fJetPt->Fill(jtpt[0]);
      corr->Fill(jtpt[0],maxl1pt);

      count++;

      if(maxl1pt>L1_THRESHOLD)
	accepted->Fill(jtpt[0]);
    }

    // for(std::vector<double>::const_iterator it = l1extra->cenJetEt.begin(); it != l1extra->cenJetEt.end(); ++it)
    // {
    //   printf("%lf\n",*it);
    // }
    // //printf("%d\n",l1_run);
    //if(j > 10) break;
  }

  TGraphAsymmErrors *a = new TGraphAsymmErrors();
  a->SetName("asymm");
  a->BayesDivide(accepted,fJetPt);

  // TCanvas *c1 = new TCanvas();
  // l1JetPt->Draw();
  // TCanvas *c2 = new TCanvas();
  // fJetPt->Draw();
  // TCanvas *c3 = new TCanvas();
  // a->Draw("p A");

  l1JetPt->Write();
  fJetPt->Write();
  corr->Write();
  accepted->Write();
  a->Write();

  printf("Matching entries: %d\n",count);

  lFile->Close();
  fFile->Close();
  outFile->Close();
}

int main()
{
  matching_l12forest();
  return 0;
}
