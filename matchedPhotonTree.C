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

#include "CalculateIsolations.cc"

const int MAXJETS = 500;
//const Double_t L1_THRESHOLD[2] = {60, 100};
//const Int_t THRESHOLDS = 2;

Long64_t makeKey(Int_t run, Int_t event){
  return (10000000000*(Long64_t)run + (Long64_t)event);
}

void matchedPhotonTree(bool montecarlo)
{
  std::cout << "MC? " << montecarlo << std::endl;
  //const char *type = "akVs3CaloJets_eta2";
  //const TString l1_input = "/export/d00/scratch/luck/L1Tree_minbias_chunk1.root";
  //const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbias_HI_and_PP_algos.root";
  //const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbias_l1ntuple_HIPUM0.root";
  //const TString l1_input = "/export/d00/scratch/luck/dijet15_l1ntuple_20141022_v2.root";
  //const TString l1_input = "/export/d00/scratch/luck/photon30_l1ntuple_20141022_v3.root";
  //const TString l1_input = "/export/d00/scratch/luck/hydjet_l1ntuple_20141022_v2.root";
  //const TString l1_input = "/export/d00/scratch/luck/jet55_data_l1ntuple_20141022.root";
  //const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/jet55_data_l1ntuple_20141022.root";
  // const TString l1_input = "/export/d00/scratch/luck/minbias_HI_v2.root";
  //const TString l1_input = "/export/d00/scratch/luck/minbias_HI_rctCalibrations_v4_partial.root";
  //const TString l1_input = "/export/d00/scratch/luck/photon_data_l1ntuple_v2.root";
  const TString l1_input = "/export/d00/scratch/luck/photon_data_l1ntuple_noHCAL.root";
  //const TString l1_input = "/export/d00/scratch/luck/minbias_hydjet_l1ntuple.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");
  //TTree *l1Tree = (TTree*)lFile->Get("HIdigis/L1UpgradeTree");
  //TTree *l1Tree = (TTree*)lFile->Get("PPdigis/L1UpgradeTree");

  Int_t l1_event, l1_run;
  Int_t l1_num;
  Int_t l1_hwPt[MAXJETS], l1_hwEta[MAXJETS], l1_hwPhi[MAXJETS], l1_hwQual[MAXJETS];
  Double_t l1_pt[MAXJETS], l1_eta[MAXJETS], l1_phi[MAXJETS];
  Int_t l1_hwIso[MAXJETS];

  Int_t l1J_num;
  Int_t l1J_hwPt[MAXJETS], l1J_hwEta[MAXJETS], l1J_hwPhi[MAXJETS], l1J_hwQual[MAXJETS];
  Double_t l1J_pt[MAXJETS], l1J_eta[MAXJETS], l1J_phi[MAXJETS];

  Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];
  Int_t emcand_hwPt[144], emcand_hwPhi[144], emcand_hwEta[144], emcand_hwIso[144];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("nJet",&l1J_num);
  l1Tree->SetBranchAddress("jet_hwPt",l1J_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1J_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1J_hwPhi);
  l1Tree->SetBranchAddress("jet_hwQual",l1J_hwQual);
  l1Tree->SetBranchAddress("jet_pt",l1J_pt);
  l1Tree->SetBranchAddress("jet_eta",l1J_eta);
  l1Tree->SetBranchAddress("jet_phi",l1J_phi);
  l1Tree->SetBranchAddress("nEgamma",&l1_num);
  l1Tree->SetBranchAddress("egamma_hwPt",l1_hwPt);
  l1Tree->SetBranchAddress("egamma_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("egamma_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("egamma_hwQual",l1_hwQual);
  l1Tree->SetBranchAddress("egamma_pt",l1_pt);
  l1Tree->SetBranchAddress("egamma_eta",l1_eta);
  l1Tree->SetBranchAddress("egamma_phi",l1_phi);
  l1Tree->SetBranchAddress("egamma_hwIso",l1_hwIso);

  l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);
  l1Tree->SetBranchAddress("emcand_hwPt", emcand_hwPt);
  l1Tree->SetBranchAddress("emcand_hwEta", emcand_hwEta);
  l1Tree->SetBranchAddress("emcand_hwPhi", emcand_hwPhi);
  l1Tree->SetBranchAddress("emcand_hwIso", emcand_hwIso);

  //const TString forest_input = "/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged/0.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim_v3.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat15_Track9_Jet30_matchEqR_merged_forest_0.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/partial_PbPb_gammaJet_MC/HiForest_QCDPhoton30.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/ginnocen/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/HiMinBias_Forest_26June2014/d9ab4aca1923b3220eacf8ee0d550950/*.root";
  const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/HiForest_PbPb_photon2030.root";
  //const TString forest_input = "/mnt/hadoop/cms/store/user/dgulhan/HiForest_HydjetMB_730_53XBS/*.root";
  // TString forest_input[10];
  // TString base = "/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/";
  // forest_input[0] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi1_*.root";
  // forest_input[1] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi2_*.root";
  // forest_input[2] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi3a_*.root";
  // forest_input[3] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi3b_*.root";
  // forest_input[4] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi4_*.root";
  // forest_input[5] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi5_*.root";
  // forest_input[6] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi6_*.root";
  // forest_input[7] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi7_*.root";
  // forest_input[8] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi8_*.root";
  // forest_input[9] = base + "HiForest_jet55or65or80_JetRAA_v1_lumi9_*.root";

  // TFile *fFile = TFile::Open(forest_input);
  // // TTree *f1Tree = (TTree*)fFile->Get("akPu3CaloJetAnalyzer/t");
  // // TTree *f2Tree = (TTree*)fFile->Get("akVs3CaloJetAnalyzer/t");
  // TTree *fTree = (TTree*)fFile->Get("multiPhotonAnalyzer/photon");
  // TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  // TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");

  // TChain *f1Tree = new TChain("akPu3CaloJetAnalyzer/t","f1Tree");
  // TChain *f2Tree = new TChain("akVs3CaloJetAnalyzer/t","f2Tree");
  TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
  TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");
  TChain *fTree = new TChain("multiPhotonAnalyzer/photon");
  TChain *fTrigTree = new TChain("hltanalysis/HltTree");

  // f1Tree->Add(forest_input[sampleNum]);
  // f2Tree->Add(forest_input[sampleNum]);
  fEvtTree->Add(forest_input);
  fSkimTree->Add(forest_input);
  fTree->Add(forest_input);
  fTrigTree->Add(forest_input);

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("vz",&vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t f_num;//, f2_num;
  Float_t f_pt[MAXJETS];//, f2_pt[MAXJETS];
  Float_t f_eta[MAXJETS];//, f2_eta[MAXJETS];
  Float_t f_phi[MAXJETS];//, f2_phi[MAXJETS];
  //Float_t f1_rawpt[MAXJETS];//, f2_rawpt[MAXJETS];

  Int_t num_gen;
  Float_t genpt[MAXJETS], geneta[MAXJETS];//, genphi[MAXJETS];

  Float_t hadronicOverEm[MAXJETS], sigmaIetaIeta[MAXJETS], sigmaIphiIphi[MAXJETS];
  Float_t cc4[MAXJETS], cr4[MAXJETS], ct4PtCut20[MAXJETS];
  Int_t isEle[MAXJETS];
  Float_t trkSumPtHollowConeDR04[MAXJETS], hcalTowerSumEtConeDR04[MAXJETS], ecalRecHitSumEtConeDR04[MAXJETS];
  Float_t swissCrx[MAXJETS], seedTime[MAXJETS];

  Float_t genIso[MAXJETS];
  Int_t genId[MAXJETS], genMomId[MAXJETS];


  // f1Tree->SetBranchAddress("nref",&f1_num);
  // f1Tree->SetBranchAddress("jtpt",f1_pt);
  // f1Tree->SetBranchAddress("jteta",f1_eta);
  // f1Tree->SetBranchAddress("jtphi",f1_phi);
  // f1Tree->SetBranchAddress("rawpt",f1_rawpt);
  // f2Tree->SetBranchAddress("nref",&f2_num);
  // f2Tree->SetBranchAddress("jtpt",f2_pt);
  // f2Tree->SetBranchAddress("jteta",f2_eta);
  // f2Tree->SetBranchAddress("jtphi",f2_phi);
  // f2Tree->SetBranchAddress("rawpt",f2_rawpt);
  fTree->SetBranchAddress("nPhotons",&f_num);
  fTree->SetBranchAddress("pt",f_pt);
  fTree->SetBranchAddress("eta",f_eta);
  fTree->SetBranchAddress("phi",f_phi);
  fTree->SetBranchAddress("cc4",cc4);
  fTree->SetBranchAddress("cr4",cr4);
  fTree->SetBranchAddress("ct4PtCut20",ct4PtCut20);
  fTree->SetBranchAddress("trkSumPtHollowConeDR04",trkSumPtHollowConeDR04);
  fTree->SetBranchAddress("hcalTowerSumEtConeDR04",hcalTowerSumEtConeDR04);
  fTree->SetBranchAddress("ecalRecHitSumEtConeDR04",ecalRecHitSumEtConeDR04);
  fTree->SetBranchAddress("hadronicOverEm",hadronicOverEm);
  fTree->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
  fTree->SetBranchAddress("isEle",isEle);
  fTree->SetBranchAddress("sigmaIphiIphi",sigmaIphiIphi);
  fTree->SetBranchAddress("swissCrx",swissCrx);
  fTree->SetBranchAddress("seedTime",seedTime);

  if(montecarlo)
  {
    // f1Tree->SetBranchAddress("ngen",&num_gen);
    // f1Tree->SetBranchAddress("genpt",genpt);
    // f1Tree->SetBranchAddress("geneta",geneta);
    // f1Tree->SetBranchAddress("genphi",genphi);

    fTree->SetBranchAddress("nGp",&num_gen);
    fTree->SetBranchAddress("gpEt",genpt);
    fTree->SetBranchAddress("gpEta",geneta);
    //fTree->SetBranchAddress("gpPhi",genphi);
    fTree->SetBranchAddress("gpCalIsoDR04",genIso);
    fTree->SetBranchAddress("gpId",genId);
    fTree->SetBranchAddress("gpMomId",genMomId);
  }

  // // triggers of interest
  // Int_t L1_SingleEG2_BptxAND, L1_SingleEG5_BptxAND, L1_SingleEG8_BptxAND, L1_SingleEG12;
  // Int_t HLT_HISinglePhoton20_v2, HLT_HISinglePhoton30_v2, HLT_HISinglePhoton40_v2;
  // Int_t HLT_HIJet55_v1;

  // fTrigTree->SetBranchAddress("L1_SingleEG2_BptxAND",&L1_SingleEG2_BptxAND);
  // fTrigTree->SetBranchAddress("L1_SingleEG5_BptxAND",&L1_SingleEG5_BptxAND);
  // fTrigTree->SetBranchAddress("L1_SingleEG8_BptxAND",&L1_SingleEG8_BptxAND);
  // fTrigTree->SetBranchAddress("L1_SingleEG12",&L1_SingleEG12);
  // fTrigTree->SetBranchAddress("HLT_HISinglePhoton20_v2",&HLT_HISinglePhoton20_v2);
  // fTrigTree->SetBranchAddress("HLT_HISinglePhoton30_v2",&HLT_HISinglePhoton30_v2);
  // fTrigTree->SetBranchAddress("HLT_HISinglePhoton40_v2",&HLT_HISinglePhoton40_v2);
  // fTrigTree->SetBranchAddress("HLT_HIJet55_v1",&HLT_HIJet55_v1);


  TFile *outFile = new TFile(Form("photon_data_noHCAL_compTree.root"),"RECREATE");
  TTree *outTree = new TTree("l1_photon_tree","l1_photon_tree");

  Int_t run, lumi, evt;

  Int_t nl1Jet;
  Int_t l1Jet_hwPt[MAXJETS], l1Jet_hwEta[MAXJETS], l1Jet_hwPhi[MAXJETS], l1Jet_hwQual[MAXJETS];
  Float_t l1Jet_pt[MAXJETS], l1Jet_eta[MAXJETS], l1Jet_phi[MAXJETS];
  // Float_t PuJet_pt[MAXJETS], PuJet_eta[MAXJETS], PuJet_phi[MAXJETS], PuJet_rawpt[MAXJETS];
  // Float_t VsJet_pt[MAXJETS], VsJet_eta[MAXJETS], VsJet_phi[MAXJETS], VsJet_rawpt[MAXJETS];
  // Float_t genJet_pt[MAXJETS], genJet_eta[MAXJETS], genJet_phi[MAXJETS];
  Bool_t goodEvent;
  Int_t hiBinOut;

  Int_t nl1Egamma, nGen;
  Int_t l1Egamma_hwPt[MAXJETS], l1Egamma_hwEta[MAXJETS], l1Egamma_hwPhi[MAXJETS], l1Egamma_hwQual[MAXJETS];
  Int_t l1Egamma_hwIso[MAXJETS];
  Float_t l1Egamma_pt[MAXJETS], l1Egamma_eta[MAXJETS], l1Egamma_phi[MAXJETS];

  Int_t emcand_hwPt_[144], emcand_hwPhi_[144], emcand_hwEta_[144], emcand_hwIso_[144];
  Int_t iso3x3[144], isoCross[144], isoFourPoint[144], isoSingle[144], isoHolePunch[144];

  Int_t region_hwPt_[396], region_hwPhi_[396], region_hwEta_[396];

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

  Float_t gen_pt[MAXJETS], gen_eta[MAXJETS];//, gen_phi[MAXJETS];
  Float_t gen_iso[MAXJETS];
  Int_t gen_id[MAXJETS];
  Int_t gen_momId[MAXJETS];

  // Int_t L1_SingleEG2_BptxAND_, L1_SingleEG5_BptxAND_, L1_SingleEG8_BptxAND_, L1_SingleEG12_;
  // Int_t HLT_HISinglePhoton20_v2_, HLT_HISinglePhoton30_v2_, HLT_HISinglePhoton40_v2_;
  // Int_t HLT_HIJet55_v1_;


  outTree->Branch("run",&run,"run/I");
  outTree->Branch("lumi",&lumi,"lumi/I");
  outTree->Branch("evt",&evt,"evt/I");

  outTree->Branch("goodEvent",&goodEvent,"goodEvent/O");
  outTree->Branch("hiBin",&hiBinOut,"hiBin/I");
  outTree->Branch("nl1Jet",&nl1Jet,"nl1Jet/I");
  outTree->Branch("l1Jet_hwPt",l1Jet_hwPt,"l1Jet_hwPt[nl1Jet]/I");
  outTree->Branch("l1Jet_hwEta",l1Jet_hwEta,"l1Jet_hwEta[nl1Jet]/I");
  outTree->Branch("l1Jet_hwPhi",l1Jet_hwPhi,"l1Jet_hwPhi[nl1Jet]/I");
  outTree->Branch("l1Jet_hwQual",l1Jet_hwQual,"l1Jet_hwQual[nl1Jet]/I");
  outTree->Branch("l1Jet_pt",l1Jet_pt,"l1Jet_pt[nl1Jet]/F");
  outTree->Branch("l1Jet_eta",l1Jet_eta,"l1Jet_eta[nl1Jet]/F");
  outTree->Branch("l1Jet_phi",l1Jet_phi,"l1Jet_phi[nl1Jet]/F");
  outTree->Branch("nl1Egamma",&nl1Egamma,"nl1Egamma/I");
  outTree->Branch("l1Egamma_hwPt",l1Egamma_hwPt,"l1Egamma_hwPt[nl1Egamma]/I");
  outTree->Branch("l1Egamma_hwEta",l1Egamma_hwEta,"l1Egamma_hwEta[nl1Egamma]/I");
  outTree->Branch("l1Egamma_hwPhi",l1Egamma_hwPhi,"l1Egamma_hwPhi[nl1Egamma]/I");
  outTree->Branch("l1Egamma_hwQual",l1Egamma_hwQual,"l1Egamma_hwQual[nl1Egamma]/I");
  outTree->Branch("l1Egamma_hwIso",l1Egamma_hwIso,"l1Egamma_hwIso[nl1Egamma]/I");
  outTree->Branch("l1Egamma_pt",l1Egamma_pt,"l1Egamma_pt[nl1Egamma]/F");
  outTree->Branch("l1Egamma_eta",l1Egamma_eta,"l1Egamma_eta[nl1Egamma]/F");
  outTree->Branch("l1Egamma_phi",l1Egamma_phi,"l1Egamma_phi[nl1Egamma]/F");

  outTree->Branch("emcand_hwPt",emcand_hwPt_,"emcand_hwPt[144]/I");
  outTree->Branch("emcand_hwPhi",emcand_hwPhi_,"emcand_hwPhi[144]/I");
  outTree->Branch("emcand_hwEta",emcand_hwEta_,"emcand_hwEta[144]/I");
  outTree->Branch("emcand_hwIso",emcand_hwIso_,"emcand_hwIso[144]/I");
  outTree->Branch("iso3x3",iso3x3,"iso3x3[144]/I");
  outTree->Branch("isoCross",isoCross,"isoCross[144]/I");
  outTree->Branch("isoFourPoint",isoFourPoint,"isoFourPoint[144]/I");
  outTree->Branch("isoSingle",isoSingle,"isoSingle[144]/I");
  outTree->Branch("isoHolePunch",isoHolePunch,"isoHolePunch[144]/I");

  outTree->Branch("region_hwPt",region_hwPt_,"region_hwPt[396]/I");
  outTree->Branch("region_hwPhi",region_hwPhi_,"region_hwPhi[396]/I");
  outTree->Branch("region_hwEta",region_hwEta_,"region_hwEta[396]/I");

  outTree->Branch("nPhoton",&nPhoton,"nPhoton/I");
  outTree->Branch("photon_pt",photon_pt,"photon_pt[nPhoton]/F");
  outTree->Branch("photon_eta",photon_eta,"photon_eta[nPhoton]/F");
  outTree->Branch("photon_phi",photon_phi,"photon_phi[nPhoton]/F");

  outTree->Branch("cc4",photon_cc4,"cc4[nPhoton]/F");
  outTree->Branch("cr4",photon_cr4,"cr4[nPhoton]/F");
  outTree->Branch("ct4PtCut20",photon_ct4PtCut20,"ct4PtCut20[nPhoton]/F");
  outTree->Branch("trkSumPtHollowConeDR04",photon_trkSumPtHollowConeDR04,"trkSumPtHollowConeDR04[nPhoton]/F");
  outTree->Branch("hcalTowerSumEtConeDR04",photon_hcalTowerSumEtConeDR04,"hcalTowerSumEtConeDR04[nPhoton]/F");
  outTree->Branch("ecalRecHitSumEtConeDR04",photon_ecalRecHitSumEtConeDR04,"ecalRecHitSumEtConeDR04[nPhoton]/F");
  outTree->Branch("hadronicOverEm",photon_hadronicOverEm,"hadronicOverEm[nPhoton]/F");
  outTree->Branch("sigmaIetaIeta",photon_sigmaIetaIeta,"sigmaIetaIeta[nPhoton]/F");
  outTree->Branch("isEle",photon_isEle,"isEle[nPhoton]/F");
  outTree->Branch("sigmaIphiIphi",photon_sigmaIphiIphi,"sigmaIphiIphi[nPhoton]/F");
  outTree->Branch("swissCrx",photon_swissCrx,"swissCrx[nPhoton]/F");
  outTree->Branch("seedTime",photon_seedTime,"seedTime[nPhoton]/F");


  if(montecarlo)
  {
    // outTree->Branch("nGenJet",&nGenJet,"nGenJet/I");
    // outTree->Branch("genJet_pt",genJet_pt,"genJet_pt[nGenJet]/F");
    // outTree->Branch("genJet_eta",genJet_eta,"genJet_eta[nGenJet]/F");
    // outTree->Branch("genJet_phi",genJet_phi,"genJet_phi[nGenJet]/F");

    outTree->Branch("nGen",&nGen,"nGen/I");
    outTree->Branch("gen_pt",gen_pt,"gen_pt[nGen]/F");
    outTree->Branch("gen_eta",gen_eta,"gen_eta[nGen]/F");
    //outTree->Branch("gen_phi",gen_phi,"gen_phi[nGen]/F");
    outTree->Branch("gen_iso",gen_iso,"gen_iso[nGen]/F");
    outTree->Branch("gen_id",gen_id,"gen_id[nGen]/I");
    outTree->Branch("gen_momId",gen_momId,"gen_momId[nGen]/I");

  }

  // outTree->Branch("L1_SingleEG2_BptxAND",&L1_SingleEG2_BptxAND_,"L1_SingleEG2_BptxAND/I");
  // outTree->Branch("L1_SingleEG5_BptxAND",&L1_SingleEG5_BptxAND_,"L1_SingleEG5_BptxAND/I");
  // outTree->Branch("L1_SingleEG8_BptxAND",&L1_SingleEG8_BptxAND_,"L1_SingleEG8_BptxAND/I");
  // outTree->Branch("L1_SingleEG12",&L1_SingleEG12_,"L1_SingleEG12/I");
  // outTree->Branch("HLT_HISinglePhoton20_v2",&HLT_HISinglePhoton20_v2_,"HLT_HISinglePhoton20_v2/I");
  // outTree->Branch("HLT_HISinglePhoton30_v2",&HLT_HISinglePhoton30_v2_,"HLT_HISinglePhoton30_v2/I");
  // outTree->Branch("HLT_HISinglePhoton40_v2",&HLT_HISinglePhoton40_v2_,"HLT_HISinglePhoton40_v2/I");
  // outTree->Branch("HLT_HIJet55_v1",&HLT_HIJet55_v1_,"HLT_HIJet55_v1/I");

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

      run = f_run;
      lumi = 0; //until I fix the analyzer
      evt = f_evt;
      hiBinOut = hiBin;

      nl1Jet = l1J_num;
      for(int i = 0; i < l1J_num; i++)
      {
	l1Jet_hwPt[i] = l1J_hwPt[i];
	l1Jet_hwEta[i] = l1J_hwEta[i];
	l1Jet_hwPhi[i] = l1J_hwPhi[i];
	l1Jet_hwQual[i] = l1J_hwQual[i];

	l1Jet_pt[i] = l1J_pt[i];
	l1Jet_eta[i] = l1J_eta[i];
	l1Jet_phi[i] = l1J_phi[i];
      }


      nl1Egamma = l1_num;
      for(int i = 0; i < l1_num; i++)
      {
	l1Egamma_hwPt[i] = l1_hwPt[i];
	l1Egamma_hwEta[i] = l1_hwEta[i];
	l1Egamma_hwPhi[i] = l1_hwPhi[i];
	l1Egamma_hwQual[i] = l1_hwQual[i];

	l1Egamma_pt[i] = l1_pt[i];
	l1Egamma_eta[i] = l1_eta[i];
	l1Egamma_phi[i] = l1_phi[i];

	l1Egamma_hwIso[i] = l1_hwIso[i];
      }

      // emcand info and iso
      Isolation isocal;
      isocal.FillRegions(region_hwPt, region_hwEta, region_hwPhi);
      for(int i = 0; i < 144; ++i)
      {
	emcand_hwPt_[i] = emcand_hwPt[i];
	emcand_hwEta_[i] = emcand_hwEta[i];
	emcand_hwPhi_[i] = emcand_hwPhi[i];
	emcand_hwIso_[i] = emcand_hwIso_[i];

	iso3x3[i] = isocal.iso3x3(emcand_hwEta_[i], emcand_hwPhi_[i]);
	isoCross[i] = isocal.isoCross(emcand_hwEta_[i], emcand_hwPhi_[i]);
	isoFourPoint[i] = isocal.isoFourPoint(emcand_hwEta_[i], emcand_hwPhi_[i]);
	isoSingle[i] = isocal.isoSingle(emcand_hwEta_[i], emcand_hwPhi_[i]);
	isoHolePunch[i] = isocal.isoHolePunch(emcand_hwEta_[i], emcand_hwPhi_[i]);
      }

      for(int i = 0; i < 396; ++i)
      {
	region_hwPt_[i] = region_hwPt[i];
	region_hwEta_[i] = region_hwEta[i];
	region_hwPhi_[i] = region_hwPhi[i];
      }

      fTrigTree->GetEntry(j);
      // L1_SingleEG2_BptxAND_ = L1_SingleEG2_BptxAND;
      // L1_SingleEG5_BptxAND_ = L1_SingleEG5_BptxAND;
      // L1_SingleEG8_BptxAND_ = L1_SingleEG8_BptxAND;
      // L1_SingleEG12_ = L1_SingleEG12;
      // HLT_HISinglePhoton20_v2_ = HLT_HISinglePhoton20_v2;
      // HLT_HISinglePhoton30_v2_ = HLT_HISinglePhoton30_v2;
      // HLT_HISinglePhoton40_v2_ = HLT_HISinglePhoton40_v2;
      // HLT_HIJet55_v1_ = HLT_HIJet55_v1;

      // f1Tree->GetEntry(j);
      // f2Tree->GetEntry(j);
      fTree->GetEntry(j);

      // nPuJet = f1_num;
      // for(int i = 0; i < f1_num; i++)
      // {
      // 	PuJet_pt[i] = f1_pt[i];
      // 	PuJet_rawpt[i] = f1_rawpt[i];
      // 	PuJet_eta[i] = f1_eta[i];
      // 	PuJet_phi[i] = f1_phi[i];
      // }

      // nVsJet = f2_num;
      // for(int i = 0; i < f2_num; i++)
      // {
      // 	VsJet_pt[i] = f2_pt[i];
      // 	VsJet_rawpt[i] = f2_rawpt[i];
      // 	VsJet_eta[i] = f2_eta[i];
      // 	VsJet_phi[i] = f2_phi[i];
      // }
      nPhoton = f_num;
      for(int i = 0; i < f_num; i++)
      {
	photon_pt[i] = f_pt[i];
	photon_eta[i] = f_eta[i];
	photon_phi[i] = f_phi[i];
	photon_cc4[i] = cc4[i];
	photon_cr4[i] = cr4[i];
	photon_ct4PtCut20[i] = ct4PtCut20[i];
	photon_trkSumPtHollowConeDR04[i] = trkSumPtHollowConeDR04[i];
	photon_hcalTowerSumEtConeDR04[i] = hcalTowerSumEtConeDR04[i];
	photon_ecalRecHitSumEtConeDR04[i] = ecalRecHitSumEtConeDR04[i];
	photon_hadronicOverEm[i] = hadronicOverEm[i];
	photon_sigmaIetaIeta[i] = sigmaIetaIeta[i];
	photon_isEle[i] = isEle[i];
	photon_sigmaIphiIphi[i] = sigmaIphiIphi[i];
	photon_swissCrx[i] = swissCrx[i];
	photon_seedTime[i] = seedTime[i];
      }

      if(montecarlo)
      {
	// nGenJet = num_gen;
	// for(int i = 0; i < num_gen; i++)
	// {
	//   genJet_pt[i] = genpt[i];
	//   genJet_eta[i] = geneta[i];
	//   genJet_phi[i] = genphi[i];
	// }

	nGen = num_gen;
	for(int i = 0; i < num_gen; i++)
	{
	  gen_pt[i] = genpt[i];
	  gen_eta[i] = geneta[i];
	  //gen_phi[i] = genphi[i];
	  gen_iso[i] = genIso[i];
	  gen_id[i] = genId[i];
	  gen_momId[i] = genMomId[i];

	}
      }

      fSkimTree->GetEntry(j);
      goodEvent = false;
      if((pcollisionEventSelection == 1) && (montecarlo || (pHBHENoiseFilter == 1)) && (TMath::Abs(vz) < 15))
      {
	goodEvent = true;
      }

      outTree->Fill();
    }
  }

  outTree->Write();

  std::cout << "Matching entries: " << count << std::endl;

  lFile->Close();
  //fFile->Close();
  outFile->Close();
}

int main(int argc, char **argv)
{
  // bool montecarlo = false;
  // int sampleNum = 0;
  // if(argc == 3)
  // {
  //   montecarlo = atoi(argv[2]);
  //   sampleNum = atoi(argv[1]);
  // }
  bool montecarlo = 0;
  if(argc == 2)
    montecarlo = atoi(argv[1]);

  matchedPhotonTree(montecarlo);
  return 0;
}
