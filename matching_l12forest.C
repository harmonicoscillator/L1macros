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
//#include <TMap.h>
#include <map>

using namespace std;
#include "l1ExtraTree.h"
#include "l1Tree.h"

const Double_t L1_THRESHOLD = 30;

long makeKey(long run, long lumi, long event){
  return (10000000000*run + 10000000*lumi + event);
}

void matching_l12forest()
{
  //const TString l1_input = "/export/d00/scratch/luck/L1Tree_minbias_chunk1.root";
  const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/L1Tree_MinBiasSkim_v1.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *lTree = (TTree*)lFile->Get("l1ExtraTreeProducer/L1ExtraTree");
  TTree *lEvtTree = (TTree*)lFile->Get("l1NtupleProducer/L1Tree");

  l1ExtraTree *l1extra = new l1ExtraTree(lTree);
  l1Tree *l1 = new l1Tree(lEvtTree);

  //const TString forest_input = "/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root";
  const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged/0.root";
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

  map<long, int> kmap;
  //TMap *kmap = new TMap();// = new unordered_map<long,int>();
  // choose loop over forest firest
  // int f_entries = fTree->GetEntries();
  // for(int j = 0; j < f_entries; ++j)
  // {
  //   fTree->GetEntry(j);
  //   long key = 10000000000*f_run + 10000000*f_lumi + f_evt;
  //   pair<long,int> p(key,j);
  //   kmap.insert(p);
  // }
  // choose loop over l1 tree first (smaller)
  int l_entries = lEvtTree->GetEntries();
  for(long j = 0; j < l_entries; ++j)
  {
    if(j % 1000 == 0)
      printf("%ld / %d\n",j,l_entries);

    l1->GetEntry(j);
    int l1_evt = l1->event;
    int l1_run = l1->run;
    int l1_lumi = l1->lumi;
    long key = makeKey(l1_run, l1_lumi, l1_evt);

    pair<long,int> p(key,j);
    kmap.insert(p);
    //kmap->Add(&key,&j);
    if(j > 5000) break;
  }

  //outFile->cd();

  const int nBins = 75;
  const double maxPt = 150;

  TH1D *l1JetPt = new TH1D("l1JetPt",";L1 Jet p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fJetPt = new TH1D("fJetPt",";akVs3PF Jet p_{T}",nBins,0,maxPt);
  TH1D *accepted = new TH1D("accepted",";akVsPF Jet p_{T}",nBins,0,maxPt);
  TH2D *corr = new TH2D("corr","akVs pt;l1 pt",nBins,0,maxPt,nBins,0,maxPt);

  int count = 0;

  int entries = lTree->GetEntries();
  for(long j = 0; j < entries; ++j)
  {
    if(j % 1000 == 0)
      printf("%ld / %d\n",j,entries);

    //l1->GetEntry(j);
    fTree->GetEntry(j);
    long key = makeKey(f_run, f_lumi, f_evt);

    // int l1_evt = l1->event;
    // int l1_run = l1->run;
    // int l1_lumi = l1->lumi;
    // long key = 10000000000*l1_run + 10000000*l1_lumi + l1_evt;

    // long *k = (long *)kmap->GetValue(&key);
    // if(k = 0){
    map<long,int>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      // fTree->GetEntry(got->second);
      // l1extra->GetEntry(j);

      l1->GetEntry(got->second);
      l1extra->GetEntry(got->second);
      //l1extra->GetEntry(*k);

      if(l1->event != f_evt)
      {
      	printf("ERROR: not actually the same event");
      	exit(1);
      }


      printf("%d \n",(int)l1extra->cenJetEt.size());
      double maxl1pt = TMath::Max(l1extra->cenJetEt[0],l1extra->fwdJetEt[0]);
      //double maxl1pt = l1extra->cenJetEt[0];
      //double maxl1pt = l1extra->fwdJetEt[0];

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

  //outFile->cd();
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
