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

#include <vector>

using namespace std;
//#include "test2.h"

void matching_l12forest()
{
  const TString l1_input = "/export/d00/scratch/luck/L1Tree_minbias_chunk1.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *lTree = (TTree*)lFile->Get("l1ExtraTreeProducer/L1ExtraTree");
  //TTree *lEvtTree = (TTree*)lFile->Get("l1NtupleProducer/L1Tree");
  lTree->AddFriend("l1NtupleProducer/L1Tree","/export/d00/scratch/luck/L1Tree_minbias_chunk1.root");

  const TString forest_input = "/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root";
  TFile *fFile = TFile::Open(forest_input);
  TTree *fTree = (TTree*)fFile->Get("akVs3PFJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);

  Int_t l1_evt;
  Int_t l1_run;
  Int_t l1_lumi;
  TBranch        *b_Event_run;   //!
  TBranch        *b_Event_event;   //!
  TBranch        *b_Event_lumi;   //!
  lTree->SetBranchAddress("event",&l1_evt,&b_Event_event);
  lTree->SetBranchAddress("run",&l1_run,&b_Event_run);
  lTree->SetBranchAddress("lumi",&l1_lumi,&b_Event_lumi);

  lTree->SetBranchStatus("*",0);
  lTree->SetBranchStatus("event",1);
  lTree->SetBranchStatus("run",1);
  lTree->SetBranchStatus("lumi",1);

  //UInt_t          nCenJets;
  vector<double>  cenJetEt;
  vector<double>  cenJetEta;
  vector<double>  cenJetPhi;
  //UInt_t          nFwdJets;
  vector<double>  fwdJetEt;
  vector<double>  fwdJetEta;
  vector<double>  fwdJetPhi;

  lTree->SetBranchAddress("cenJetEta",&cenJetEta);
  lTree->SetBranchAddress("cenJetPhi",&cenJetPhi);
  lTree->SetBranchAddress("fwdJetEt",&fwdJetEt);
  lTree->SetBranchAddress("fwdJetEta",&fwdJetEta);
  lTree->SetBranchAddress("fwdJetPhi",&fwdJetPhi);

  // Int_t f_evt, f_run, f_lumi;
  // fTree->SetBranchAddress("evt",&f_evt);
  // fTree->SetBranchAddress("run",&f_run);
  // fTree->SetBranchAddress("lumi",&f_lumi);

  // Int_t nref;
  // Float_t jtpt[500], jteta[500], jtphi[500];
  // fTree->SetBranchAddress("nref",&nref);
  // fTree->SetBranchAddress("jtpt",jtpt);
  // fTree->SetBranchAddress("jteta",jteta);
  // fTree->SetBranchAddress("jtphi",jtphi);

  lTree->GetEntry(50);
  //lEvtTree->GetEntry(50);

  printf("%d\n",l1_evt);
  printf("%d\n",l1_run);
  printf("%d\n",l1_lumi);

  // int entries = lEvtTree->GetEntries();
  // for(int j = 0; j < entries; ++j)
  // {
  //   lEvtTree->GetEntry(j);
  //   for(std::vector<double>::const_iterator it = cenJetEt->begin(); it != cenJetEt->end(); ++it)
  //   {
  //     printf("%lf\n",*it);
  //   }
  //   printf("%d\n",l1_run);
  //   if(j > 100) break;
  // }
  lFile->Close();
  fFile->Close();
}

int main()
{
  matching_l12forest();
  return 0;
}
