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

const Double_t L1_THRESHOLD[4] = {0,15,30,60};

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
  TTree *fTree = (TTree*)fFile->Get("akVs3CaloJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);
  TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");
  fTree->AddFriend(fSkimTree);

  Int_t f_evt, f_run, f_lumi;
  fTree->SetBranchAddress("evt",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t nref;
  Float_t jtpt[500], jteta[500], jtphi[500];
  fTree->SetBranchAddress("nref",&nref);
  fTree->SetBranchAddress("jtpt",jtpt);
  fTree->SetBranchAddress("jteta",jteta);
  fTree->SetBranchAddress("jtphi",jtphi);

  TFile *outFile = new TFile("hist_out_akVs3Calo.root","RECREATE");

  map<long, int> kmap;
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
    // if(j % 1000 == 0)
    //   printf("%ld / %d\n",j,l_entries);

    l1->GetEntry(j);
    long key = makeKey(l1->run, l1->lumi, l1->event);

    pair<long,int> p(key,j);
    kmap.insert(p);
  }

  outFile->cd();

  const int nBins = 75;
  const double maxPt = 150;

  TH1D *l1JetPt = new TH1D("l1JetPt",";L1 Jet p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fJetPt = new TH1D("fJetPt",";offline Jet p_{T}",nBins,0,maxPt);
  TH1D *accepted[4];
  accepted[0] = new TH1D("accepted_pt0",";offline Jet p_{T}",nBins,0,maxPt);
  accepted[1] = (TH1D*)accepted[0]->Clone("accepted_pt15");
  accepted[2] = (TH1D*)accepted[0]->Clone("accepted_pt30");
  accepted[3] = (TH1D*)accepted[0]->Clone("accepted_pt60");
  TH2D *corr = new TH2D("corr","akVs pt;l1 pt",nBins,0,maxPt,nBins,0,maxPt);

  int count = 0;

  int entries = fTree->GetEntries();
  for(long j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%ld / %d\n",j,entries);

    //l1->GetEntry(j);
    fTree->GetEntry(j);
    long key = makeKey(f_run, f_lumi, f_evt);

    map<long,int>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      // fTree->GetEntry(got->second);
      // l1extra->GetEntry(j);

      //l1->GetEntry(got->second);
      l1extra->GetEntry(got->second);
      //l1extra->GetEntry(*k);

      // if(l1->event != f_evt)
      // {
      // 	printf("ERROR: not actually the same event");
      // 	exit(1);
      // }

      double maxCenJetEt = 0;
      if(l1extra->nCenJets != 0)
	maxCenJetEt = l1extra->cenJetEt[0];
      double maxFwdJetEt = 0;
      if(l1extra->nFwdJets != 0)
	maxFwdJetEt = l1extra->fwdJetEt[0];

      double maxl1pt = TMath::Max(maxCenJetEt, maxFwdJetEt);
      //double maxl1pt = l1extra->cenJetEt[0];
      //double maxl1pt = l1extra->fwdJetEt[0];

      double maxjtpt = 0;
      if(nref > 0)
	maxjtpt = jtpt[0];

      l1JetPt->Fill(maxl1pt);
      if((pcollisionEventSelection == 1) && (pHBHENoiseFilter == 1))
      {
	fJetPt->Fill(maxjtpt);
	corr->Fill(maxjtpt,maxl1pt);

	for(int k = 0; k < 4; ++k){
	  if(maxl1pt>L1_THRESHOLD[k])
	    accepted[k]->Fill(maxjtpt);
	}

      }
      count++;
    }

    // for(std::vector<double>::const_iterator it = l1extra->cenJetEt.begin(); it != l1extra->cenJetEt.end(); ++it)
    // {
    //   printf("%lf\n",*it);
    // }
    // //printf("%d\n",l1_run);
    //if(j > 10) break;
  }

  TGraphAsymmErrors *a[4];
  for(int k = 0; k < 4; ++k){
    a[k] = new TGraphAsymmErrors();
    a[k]->BayesDivide(accepted[k],fJetPt);
  }
  a[0]->SetName("asymm_pt_0");
  a[1]->SetName("asymm_pt_15");
  a[2]->SetName("asymm_pt_30");
  a[3]->SetName("asymm_pt_60");

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
  for(int k = 0; k < 4; ++k){
    accepted[k]->Write();
    a[k]->Write();
  }

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
