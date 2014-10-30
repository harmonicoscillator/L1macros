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
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <map>

const int MAXJETS = 500;
const Int_t THRESHOLDS = 11;
const Double_t L1_THRESHOLD[THRESHOLDS] = {12,16,20,24,28,32,40,52,60,80,100};

void makeTurnPhotonOn_fromTree()
{
  TFile *inFile = TFile::Open("photon30_compTree.root");
  //TFile *inFile = TFile::Open("/export/d00/scratch/luck/jet55_data_compTree_combined.root");
  TTree *inTree = (TTree*)inFile->Get("l1_photon_tree");

  Int_t run, lumi, evt;

  Int_t nl1Egamma, nGen;
  Int_t l1Egamma_hwPt[MAXJETS], l1Egamma_hwEta[MAXJETS], l1Egamma_hwPhi[MAXJETS], l1Egamma_hwQual[MAXJETS];
  Int_t l1Egamma_hwIso[MAXJETS];
  Float_t l1Egamma_pt[MAXJETS], l1Egamma_eta[MAXJETS], l1Egamma_phi[MAXJETS];

  Int_t nPhoton;
  Float_t photon_pt[MAXJETS];
  Float_t photon_eta[MAXJETS];
  Float_t photon_phi[MAXJETS];
  Float_t photon_cc4[MAXJETS];
  Float_t photon_cr4[MAXJETS];
  Float_t photon_ct4PtCut20[MAXJETS];
  Float_t photon_trkSumPtHollowConeDR04[MAXJETS];
  Float_t photon_hcalTowerSumEtConeDR04[MAXJETS];
  Float_t photon_ecalRecHitSumEtConeDR04[MAXJETS];
  Float_t photon_hadronicOverEm[MAXJETS];
  Float_t photon_sigmaIetaIeta[MAXJETS];
  Float_t photon_isEle[MAXJETS];
  Float_t photon_sigmaIphiIphi[MAXJETS];
  Float_t photon_swissCrx[MAXJETS];
  Float_t photon_seedTime[MAXJETS];

  Bool_t goodEvent;
  Int_t hiBin;

  inTree->SetBranchAddress("run",&run);
  inTree->SetBranchAddress("lumi",&lumi);
  inTree->SetBranchAddress("evt",&evt);

  inTree->SetBranchAddress("goodEvent",&goodEvent);
  inTree->SetBranchAddress("hiBin",&hiBin);

  inTree->SetBranchAddress("nl1Egamma",&nl1Egamma);
  inTree->SetBranchAddress("l1Egamma_hwPt",l1Egamma_hwPt);
  inTree->SetBranchAddress("l1Egamma_hwEta",l1Egamma_hwEta);
  inTree->SetBranchAddress("l1Egamma_hwPhi",l1Egamma_hwPhi);
  inTree->SetBranchAddress("l1Egamma_hwQual",l1Egamma_hwQual);
  inTree->SetBranchAddress("l1Egamma_hwIso",l1Egamma_hwIso);
  inTree->SetBranchAddress("l1Egamma_pt",l1Egamma_pt);
  inTree->SetBranchAddress("l1Egamma_eta",l1Egamma_eta);
  inTree->SetBranchAddress("l1Egamma_phi",l1Egamma_phi);

  inTree->SetBranchAddress("nPhoton",&nPhoton);
  inTree->SetBranchAddress("photon_pt",photon_pt);
  inTree->SetBranchAddress("photon_eta",photon_eta);
  inTree->SetBranchAddress("photon_phi",photon_phi);
 
  inTree->SetBranchAddress("cc4",photon_cc4);
  inTree->SetBranchAddress("cr4",photon_cr4);
  inTree->SetBranchAddress("ct4PtCut20",photon_ct4PtCut20);
  inTree->SetBranchAddress("trkSumPtHollowConeDR04",photon_trkSumPtHollowConeDR04);
  inTree->SetBranchAddress("hcalTowerSumEtConeDR04",photon_hcalTowerSumEtConeDR04);
  inTree->SetBranchAddress("ecalRecHitSumEtConeDR04",photon_ecalRecHitSumEtConeDR04);
  inTree->SetBranchAddress("hadronicOverEm",photon_hadronicOverEm);
  inTree->SetBranchAddress("sigmaIetaIeta",photon_sigmaIetaIeta);
  inTree->SetBranchAddress("isEle",photon_isEle);
  inTree->SetBranchAddress("sigmaIphiIphi",photon_sigmaIphiIphi);
  inTree->SetBranchAddress("swissCrx",photon_swissCrx);
  inTree->SetBranchAddress("seedTime",photon_seedTime);

  Float_t gen_pt[MAXJETS], gen_eta[MAXJETS];//, gen_phi[MAXJETS];
  Float_t gen_iso[MAXJETS];
  Int_t gen_id[MAXJETS], gen_momId[MAXJETS];

  inTree->SetBranchAddress("nGen",&nGen);
  inTree->SetBranchAddress("gen_pt",gen_pt);
  inTree->SetBranchAddress("gen_eta",gen_eta);
  //inTree->SetBranchAddress("gen_phi",gen_phi);
  inTree->SetBranchAddress("gen_iso",gen_iso);
  inTree->SetBranchAddress("gen_id",gen_id);
  inTree->SetBranchAddress("gen_momId",gen_momId);


  TFile *outFile = new TFile(Form("hist_photon30_gen.root"),"RECREATE");
  outFile->cd();

  const int nBins = 25;
  const double maxPt = 100;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[2];
  fPt[0] = new TH1D("fPt_cen",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_periph");
  TH1D *accepted[THRESHOLDS][2];

  for(int i = 0; i < THRESHOLDS; ++i)
    for(int j = 0; j < 2; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",(int)L1_THRESHOLD[i],j),";offline p_{T}",nBins,0,maxPt);
    }

  TH2D *corr = new TH2D("corr",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);

  Long64_t entries = inTree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    inTree->GetEntry(j);

    double maxisopt = l1Egamma_pt[0];
    //double maxnonisopt = l1Egamma_pt[4];
    //double maxl1pt = std::max(maxisopt, maxnonisopt);
    double maxl1pt = maxisopt;

    double maxfpt = -1;
    for(int i = 0; i < nGen; ++i)
    {
      if(TMath::Abs(gen_momId[i]) <= 22)
      if(gen_id[i] == 22)
      if(gen_iso[i] < 5)
      if(gen_pt[i] > maxfpt)
	maxfpt = gen_pt[i];
    }
    //if(f_num > MAXJETS) std::cout << "TOO SMALL" << std::endl;
    l1Pt->Fill(maxl1pt);

    if(goodEvent)
    {
      if(hiBin < 60)
	fPt[0]->Fill(maxfpt);
      else if (hiBin >= 100)
	fPt[1]->Fill(maxfpt);

      corr->Fill(maxfpt,maxl1pt);

      for(int k = 0; k < THRESHOLDS; ++k)
      {
	if(maxl1pt>L1_THRESHOLD[k])
	{
	  if(hiBin < 60)
	    accepted[k][0]->Fill(maxfpt);
	  else if (hiBin >= 100)
	    accepted[k][1]->Fill(maxfpt);
	}
      }
    }
  }

  TGraphAsymmErrors *a[4][2];
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 2; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
      a[k][l]->SetName(Form("asymm_pt_%d_%d",(int)L1_THRESHOLD[k],l));
    }
  }

  l1Pt->Write();
  fPt[0]->Write();
  fPt[1]->Write();
  corr->Write();
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 2; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
  }

  inFile->Close();
  outFile->Close();
}

int main()
{
  makeTurnPhotonOn_fromTree();
  return 0;
}
