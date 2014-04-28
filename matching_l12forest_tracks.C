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
  TTree *fTree = (TTree*)fFile->Get("anaTrack/trackTree");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);
  TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");
  fTree->AddFriend(fSkimTree);

  Int_t f_evt, f_run, f_lumi, hiBin;
  fTree->SetBranchAddress("evt",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);
  fTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  const unsigned int MAXTRACKS = 20000;
  Int_t nTrk;
  Float_t trkPt[MAXTRACKS], trkPtError[MAXTRACKS];
  Bool_t highPurity[MAXTRACKS];
  Float_t trkDxy1[MAXTRACKS], trkDxyError1[MAXTRACKS], trkDz1[MAXTRACKS], trkDzError1[MAXTRACKS];

  fTree->SetBranchAddress("nTrk",&nTrk);
  fTree->SetBranchAddress("trkPt",trkPt);
  fTree->SetBranchAddress("trkPtError",trkPtError);
  fTree->SetBranchAddress("highPurity",highPurity);
  fTree->SetBranchAddress("trkDxy1",trkDxy1);
  fTree->SetBranchAddress("trkDxyError1",trkDxyError1);
  fTree->SetBranchAddress("trkDz1",trkDz1);
  fTree->SetBranchAddress("trkDzError1",trkDzError1);

  TFile *outFile = new TFile("hist_out_tracks.root","RECREATE");

  map<long, int> kmap;
  int l_entries = lEvtTree->GetEntries();
  for(long j = 0; j < l_entries; ++j)
  {
    l1->GetEntry(j);
    long key = makeKey(l1->run, l1->lumi, l1->event);

    pair<long,int> p(key,j);
    kmap.insert(p);
  }

  outFile->cd();

  const int nBins = 150;
  const double maxPt = 300;

  TH1D *l1TrkPt = new TH1D("l1TrkPt",";L1 Trk p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fTrkPt[2];
  fTrkPt[0] = new TH1D("fTrkPt_cen",";offline Trk p_{T}",nBins,0,maxPt);
  fTrkPt[1] = (TH1D*)fTrkPt[0]->Clone("fTrkPt_periph");
  TH1D *accepted[4][2];
  accepted[0][0] = new TH1D("accepted_pt0_cen",";offline Trk p_{T}",nBins,0,maxPt);
  accepted[1][0] = (TH1D*)accepted[0][0]->Clone("accepted_pt15_cen");
  accepted[2][0] = (TH1D*)accepted[0][0]->Clone("accepted_pt30_cen");
  accepted[3][0] = (TH1D*)accepted[0][0]->Clone("accepted_pt60_cen");
  accepted[0][1] = (TH1D*)accepted[0][0]->Clone("accepted_pt0_periph");
  accepted[1][1] = (TH1D*)accepted[0][0]->Clone("accepted_pt15_periph");
  accepted[2][1] = (TH1D*)accepted[0][0]->Clone("accepted_pt30_periph");
  accepted[3][1] = (TH1D*)accepted[0][0]->Clone("accepted_pt60_periph");
  TH2D *corr = new TH2D("corr","pt;l1 pt",nBins,0,maxPt,nBins,0,maxPt);

  int count = 0;

  int entries = fTree->GetEntries();
  for(long j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%ld / %d\n",j,entries);

    fTree->GetEntry(j);
    long key = makeKey(f_run, f_lumi, f_evt);

    map<long,int>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      l1extra->GetEntry(got->second);
      kmap.erase(key);

      Float_t maxl1pt = 0;
      if(l1extra->nTauJets > 0)
	maxl1pt = l1extra->tauJetEt[0];

      Float_t maxtrkpt = 0;
      for(int i = 0; i < nTrk; ++i)
      {
	bool goodTrack = (highPurity[i] &&
			  TMath::Abs(trkDz1[i]/trkDzError1[i]) < 3 &&
			  TMath::Abs(trkDxy1[i]/trkDxyError1[i]) < 3 &&
			  trkPtError[i]/trkPt[i] < 0.05);

	if(goodTrack)
	{
	  if(trkPt[i] > maxtrkpt)
	    maxtrkpt = trkPt[i];
	}
      }

      l1TrkPt->Fill(maxl1pt);
      if((pcollisionEventSelection == 1) && (pHBHENoiseFilter == 1))
      {
	if(hiBin < 60)
	  fTrkPt[0]->Fill(maxtrkpt);
	else if (hiBin >= 100)
	  fTrkPt[1]->Fill(maxtrkpt);

	corr->Fill(maxtrkpt,maxl1pt);

	for(int k = 0; k < 4; ++k){
	  if(maxl1pt>L1_THRESHOLD[k])
	  {
	    if(hiBin < 60)
	      accepted[k][0]->Fill(maxtrkpt);
	    else if (hiBin >= 100)
	      accepted[k][1]->Fill(maxtrkpt);
	  }
	}

      }
      count++;
    }

    // for(std::vector<double>::const_iterator it = l1extra->cenTrkEt.begin(); it != l1extra->cenTrkEt.end(); ++it)
    // {
    //   printf("%lf\n",*it);
    // }
    // //printf("%d\n",l1_run);
    //if(j > 10) break;
  }

  TGraphAsymmErrors *a[4][2];
  for(int k = 0; k < 4; ++k){
    for(int l = 0; l < 2; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fTrkPt[l]);
    }
  }
  a[0][0]->SetName("asymm_pt_0_cen");
  a[1][0]->SetName("asymm_pt_15_cen");
  a[2][0]->SetName("asymm_pt_30_cen");
  a[3][0]->SetName("asymm_pt_60_cen");
  a[0][1]->SetName("asymm_pt_0_periph");
  a[1][1]->SetName("asymm_pt_15_periph");
  a[2][1]->SetName("asymm_pt_30_periph");
  a[3][1]->SetName("asymm_pt_60_periph");


  // TCanvas *c1 = new TCanvas();
  // l1TrkPt->Draw();
  // TCanvas *c2 = new TCanvas();
  // fTrkPt->Draw();
  // TCanvas *c3 = new TCanvas();
  // a->Draw("p A");

  //outFile->cd();
  l1TrkPt->Write();
  fTrkPt[0]->Write();
  fTrkPt[1]->Write();
  corr->Write();
  for(int k = 0; k < 4; ++k){
    for(int l = 0; l < 2; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
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
