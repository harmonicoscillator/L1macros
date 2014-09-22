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

#include <vector>
#include <iostream>
#include <map>

const int MAXJETS = 500;
const Double_t L1_THRESHOLD[2] = {60, 100};
const Int_t THRESHOLDS = 2;

Long64_t makeKey(Int_t run, Int_t event){
  return (10000000000*(Long64_t)run + (Long64_t)event);
}

void matching_l12forest()
{
  const char *type = "akVs3CaloJets";
  //const TString l1_input = "/export/d00/scratch/luck/L1Tree_minbias_chunk1.root";
  const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbias_HI_and_PP_algos.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("HIdigis/L1UpgradeTree");
  //TTree *l1Tree = (TTree*)lFile->Get("PPdigis/L1UpgradeTree");

  Int_t l1_event, l1_run;
  Int_t l1_num;
  //Int_t l1_hwPt[MAXJETS], l1_hwEta[MAXJETS], l1_hwPhi[MAXJETS], l1_hwQual[MAXJETS];
  Double_t l1_pt[MAXJETS];//, l1_eta[MAXJETS], l1_phi[MAXJETS];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("nJet",&l1_num);
  // l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  // l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  // l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  // l1Tree->SetBranchAddress("jet_hwQual",l1_hwQual);
  l1Tree->SetBranchAddress("jet_pt",l1_pt);
  // l1Tree->SetBranchAddress("jet_eta",l1_eta);
  // l1Tree->SetBranchAddress("jet_phi",l1_phi);
  //l1Tree->SetBranchAddress("nEgamma",&l1_num);
  //l1Tree->SetBranchAddress("egamma_pt",l1_pt);

  //const TString forest_input = "/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged/0.root";
  const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim.root";
  TFile *fFile = TFile::Open(forest_input);
  //TTree *fTree = (TTree*)fFile->Get("akPu3CaloJetAnalyzer/t");
  TTree *fTree = (TTree*)fFile->Get("akVs3CaloJetAnalyzer/t");
  //TTree *fTree = (TTree*)fFile->Get("multiPhotonAnalyzer/photon");
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");

  Int_t f_evt, f_run, f_lumi;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t f_num;
  Float_t f_pt[MAXJETS];
  //Float_t f_eta[MAXJETS], f_phi[MAXJETS];

  //Float_t hadronicOverEm[MAXJETS], sigmaIetaIeta[MAXJETS], sigmaIphiIphi[MAXJETS];
  //Float_t cc4[MAXJETS], cr4[MAXJETS], ct4PtCut20[MAXJETS];
  //Int_t isEle[MAXJETS];
  //Float_t trkSumPtHollowConeDR04[MAXJETS], hcalTowerSumEtConeDR04[MAXJETS], ecalRecHitSumEtConeDR04[MAXJETS];
  //Float_t swissCrx[MAXJETS], seedTime[MAXJETS];
  fTree->SetBranchAddress("nref",&f_num);
  fTree->SetBranchAddress("jtpt",f_pt);
  // fTree->SetBranchAddress("jteta",f_eta);
  // fTree->SetBranchAddress("jtphi",f_phi);
  // fTree->SetBranchAddress("nPhotons",&f_num);
  // fTree->SetBranchAddress("pt",f_pt);
  // fTree->SetBranchAddress("eta",f_eta);
  // fTree->SetBranchAddress("phi",f_phi);
  // fTree->SetBranchAddress("cc4",cc4);
  // fTree->SetBranchAddress("cr4",cr4);
  // fTree->SetBranchAddress("ct4PtCut20",ct4PtCut20);
  // fTree->SetBranchAddress("trkSumPtHollowConeDR04",trkSumPtHollowConeDR04);
  // fTree->SetBranchAddress("hcalTowerSumEtConeDR04",hcalTowerSumEtConeDR04);
  // fTree->SetBranchAddress("ecalRecHitSumEtConeDR04",ecalRecHitSumEtConeDR04);
  // fTree->SetBranchAddress("hadronicOverEm",hadronicOverEm);
  // fTree->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
  // fTree->SetBranchAddress("isEle",isEle);
  // fTree->SetBranchAddress("sigmaIphiIphi",sigmaIphiIphi);
  // fTree->SetBranchAddress("swissCrx",swissCrx);
  // fTree->SetBranchAddress("seedTime",seedTime);

  TFile *outFile = new TFile(Form("hist_out_%s_cleaned.root",type),"RECREATE");

  std::map<Long64_t, Long64_t> kmap;

  // choose loop over l1 tree first (smaller)
  std::cout << "Begin making map." << std::endl;
  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);
    Long64_t key = makeKey(l1_run, l1_event);
    //std::cout << l1_run << "\t" << l1_event << "\t" << key << std::endl;

    std::pair<Long64_t,Long64_t> p(key,j);
    kmap.insert(p);
  }
  std::cout << "Finished making map." << std::endl;

  outFile->cd();

  const int nBins = 75;
  const double maxPt = 300;

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

  int count = 0;

  Long64_t entries = fTree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    fEvtTree->GetEntry(j);
    Long64_t key = makeKey(f_run, f_evt);
    //std::cout << f_run << "\t" << f_evt << "\t" << key << std::endl;

    std::map<Long64_t,Long64_t>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      l1Tree->GetEntry(got->second);
      kmap.erase(key);
      count++;

      double maxl1pt = -1;
      if(l1_num > 0)
	maxl1pt = l1_pt[0];
      //if(l1_num > MAXJETS) std::cout << "TOO SMALL" << std::endl;

      fTree->GetEntry(j);
      double maxfpt = -1;
      if(f_num > 0)
      {
	for(int i = 0; i < f_num; ++i)
	{
	  // if( (hadronicOverEm[i] < 0.1) && (!isEle[i]) &&
	  //     (swissCrx[i] < 0.9) && (seedTime[i] < 3) &&
	  //     (sigmaIphiIphi[i] > 0.002) && (sigmaIetaIeta[i] > 0.002))
	      maxfpt = f_pt[i]; break;
	}
      }
      //if(f_num > MAXJETS) std::cout << "TOO SMALL" << std::endl;
      l1Pt->Fill(maxl1pt);

      fSkimTree->GetEntry(j);
      if((pcollisionEventSelection == 1) && (pHBHENoiseFilter == 1))
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

  std::cout << "Matching entries: " << count << std::endl;

  lFile->Close();
  fFile->Close();
  outFile->Close();
}

int main()
{
  matching_l12forest();
  return 0;
}
