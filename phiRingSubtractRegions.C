#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <algorithm>

struct cand{
  int pt;
  int eta;
  int phi;
};

void SlidingWindowJetFinder(cand input[396], cand output[8]);

int deltaGctPhi(int phi1, int phi2)
{
  int diff = phi1 - phi2;
  if (std::abs(phi1 - phi2) == 17) { //18 regions in phi
    diff = -diff/std::abs(diff);
  }
  return diff;
}

void phiRingSubtractRegions()
{
  const TString l1_input = "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root";
  //const TString l1_input = "/export/d00/scratch/luck/Hydjet1p8_2760GeV_L1UpgradeAnalyzer_GT_run1_mc_HIon_L1UpgradeAnalyzer.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_event, l1_run, l1_lumi;
  Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];
  cand subRegions[396];
  cand outJets[8];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);

  TFile *outFile = new TFile(Form("Hydjet502_regionsub_moveForwardDef.root"),"RECREATE");
  TTree *outTree = new TTree("phi_ring_subtracted_tree","phi_ring_subtracted_tree");

  Int_t run, lumi, evt;

  Int_t region_hwPt_[396], region_hwPhi_[396], region_hwEta_[396];
  Int_t subregion_hwPt[396], subregion_hwPhi[396], subregion_hwEta[396];
  Int_t jet_hwPt[8], jet_hwEta[8], jet_hwPhi[8], jet_pt[8];

  int puLevelHI[22];
  float r_puLevelHI[22];

  outTree->Branch("run",&run,"run/I");
  outTree->Branch("lumi",&lumi,"lumi/I");
  outTree->Branch("evt",&evt,"evt/I");

  outTree->Branch("region_hwPt",region_hwPt_,"region_hwPt[396]/I");
  outTree->Branch("region_hwPhi",region_hwPhi_,"region_hwPhi[396]/I");
  outTree->Branch("region_hwEta",region_hwEta_,"region_hwEta[396]/I");
  outTree->Branch("subregion_hwPt",subregion_hwPt,"subregion_hwPt[396]/I");
  outTree->Branch("subregion_hwPhi",subregion_hwPhi,"subregion_hwPhi[396]/I");
  outTree->Branch("subregion_hwEta",subregion_hwEta,"subregion_hwEta[396]/I");

  outTree->Branch("jet_hwPt",jet_hwPt,"jet_hwPt[8]/I");
  outTree->Branch("jet_pt",jet_pt,"jet_pt[8]/I");
  outTree->Branch("jet_hwPhi",jet_hwPhi,"jet_hwPhi[8]/I");
  outTree->Branch("jet_hwEta",jet_hwEta,"jet_hwEta[8]/I");

  outTree->Branch("puLevelHI",puLevelHI,"puLevelHI[22]/I");
  outTree->Branch("r_puLevelHI",r_puLevelHI,"r_puLevelHI[22]/F");

  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);

    run = l1_run;
    lumi = l1_lumi;
    evt = l1_event;

    for(int i = 0; i < 396; ++i)
    {
      region_hwPt_[i] = region_hwPt[i];
      region_hwEta_[i] = region_hwEta[i];
      region_hwPhi_[i] = region_hwPhi[i];
    }

    /// --------------- For heavy ion -------------------------------------
    int etaCount[22];
    for(unsigned i = 0; i < 22; ++i)
    {
      puLevelHI[i] = 0;
      r_puLevelHI[i] = 0.0;
      etaCount[i] = 0;
    }

    for(int i = 0; i < 396; ++i){
      r_puLevelHI[region_hwEta[i]] += region_hwPt[i];
      etaCount[region_hwEta[i]]++;
    }

    for(unsigned i = 0; i < 22; ++i)
    {
      if(etaCount[i] != 18)
	std::cout << "ERROR: wrong number of regions in phi ring." << std::endl;
      puLevelHI[i] = floor(r_puLevelHI[i]/18 + 0.5);
    }

    for(int i = 0; i < 396; ++i){
      subregion_hwPt[i] = std::max(0, region_hwPt[i] - puLevelHI[region_hwEta[i]]);
      subregion_hwEta[i] = region_hwEta[i];
      subregion_hwPhi[i] = region_hwPhi[i];

      subRegions[i].pt = subregion_hwPt[i];
      subRegions[i].eta = subregion_hwEta[i] ;
      subRegions[i].phi = subregion_hwPhi[i];
    }
    // end copy of algorithm from Emulator

    SlidingWindowJetFinder(subRegions, outJets);

    for(int i = 0; i < 8; i++)
    {
      jet_hwPt[i] = outJets[i].pt;
      jet_pt[i] = outJets[i].pt * 4;
      jet_hwEta[i] = outJets[i].eta;
      jet_hwPhi[i] = outJets[i].phi;
    }

    outTree->Fill();
  }

  outTree->Write();

  lFile->Close();
  outFile->Close();
}

void SlidingWindowJetFinder(cand region[396], cand output[8])
{
  std::vector<cand> forjets;
  std::vector<cand> cenjets;

  for(int i = 0; i < 396; i++) {
    int regionET = region[i].pt;
    int regionEta = region[i].eta;
    int regionPhi = region[i].phi;
    int neighborN_et = 0;
    int neighborS_et = 0;
    int neighborE_et = 0;
    int neighborW_et = 0;
    int neighborNE_et = 0;
    int neighborSW_et = 0;
    int neighborNW_et = 0;
    int neighborSE_et = 0;
    unsigned int nNeighbors = 0;
    for(int j = 0; j < 396; j++) {
      int neighborET = region[j].pt;
      int neighborEta = region[j].eta;
      int neighborPhi = region[j].phi;
      if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	 (regionEta ) == neighborEta) {
	neighborN_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta    ) == neighborEta) {
	neighborS_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
	      (regionEta + 1) == neighborEta) {
	neighborE_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
	      (regionEta - 1) == neighborEta) {
	neighborW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	      (regionEta + 1) == neighborEta) {
	neighborNE_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta - 1) == neighborEta) {
	neighborSW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	      (regionEta - 1) == neighborEta) {
	neighborNW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta + 1) == neighborEta) {
	neighborSE_et = neighborET;
	nNeighbors++;
	continue;
      }
    }
    if(regionET > neighborN_et &&
       regionET > neighborNW_et &&
       regionET > neighborW_et &&
       regionET > neighborSW_et &&
       regionET >= neighborNE_et &&
       regionET >= neighborE_et &&
       regionET >= neighborSE_et &&
       regionET >= neighborS_et) {
      unsigned int jetET = regionET +
	neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;

      int jetPhi = regionPhi;
      int jetEta = regionEta;

      bool neighborCheck = (nNeighbors == 8);
      // On the eta edge we only expect 5 neighbors
      if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
	neighborCheck = true;

      if (!neighborCheck) {
	std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
      }

      cand theJet;
      theJet.pt = jetET / 8;
      theJet.eta = jetEta;
      theJet.phi = jetPhi;

      //first iteration, eta cut defines forward
      const bool forward = (jetEta <= 4 || jetEta >= 17);
      //const bool forward = (jetEta < 4 || jetEta > 17);

      if(forward)
	forjets.push_back(theJet);
      else
	cenjets.push_back(theJet);
    }
  }

  auto comp = [&](cand i, cand j)-> bool {
    return (i.pt > j.pt );
  };

  std::sort(forjets.begin(), forjets.end(), comp);
  std::sort(cenjets.begin(), cenjets.end(), comp);
  forjets.resize(4);
  cenjets.resize(4);

  output[0]=cenjets.at(0);
  output[1]=cenjets.at(1);
  output[2]=cenjets.at(2);
  output[3]=cenjets.at(3);
  output[4]=forjets.at(0);
  output[5]=forjets.at(1);
  output[6]=forjets.at(2);
  output[7]=forjets.at(3);
}


int main()
{
  phiRingSubtractRegions();
  return 0;
}
