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
const Int_t THRESHOLDS = 8;
const Double_t L1_THRESHOLD[THRESHOLDS] = {2, 4, 5, 6, 7, 8, 10, 12};

double getPhysicalEta(int gtEta)
{
  int etaIndex = gtEta;
  const double rgnEtaValues[11] = {
    0.174, // HB and inner HE bins are 0.348 wide
    0.522,
    0.870,
    1.218,
    1.566,
    1.956, // Last two HE bins are 0.432 and 0.828 wide
    2.586,
    3.250, // HF bins are 0.5 wide
    3.750,
    4.250,
    4.750
  };
  if(etaIndex < 11) {
    return -rgnEtaValues[-(etaIndex - 10)]; // 0-10 are negative eta values
  }
  else if (etaIndex < 22) {
    return rgnEtaValues[etaIndex - 11]; // 11-21 are positive eta values
  }
  return -9;
}
double getPhysicalPhi(int phiIndex)
{
  if (phiIndex < 10)
    return 2. * M_PI * phiIndex / 18.;
  if (phiIndex < 18)
    return -M_PI + 2. * M_PI * (phiIndex - 9) / 18.;
  return -9;
}

void makeTurnPhotonOn_fromTree()
{
  TFile *inFile = TFile::Open("photon_data_compTree.root");
  //TFile *inFile = TFile::Open("/export/d00/scratch/luck/jet55_data_compTree_combined.root");
  TTree *inTree = (TTree*)inFile->Get("l1_photon_tree");

  Int_t run, lumi, evt;

  Int_t nl1Egamma;
  Int_t l1Egamma_hwPt[MAXJETS], l1Egamma_hwEta[MAXJETS], l1Egamma_hwPhi[MAXJETS], l1Egamma_hwQual[MAXJETS];
  Int_t l1Egamma_hwIso[MAXJETS];
  Float_t l1Egamma_pt[MAXJETS], l1Egamma_eta[MAXJETS], l1Egamma_phi[MAXJETS];

  Int_t nl1Jet;
  Int_t l1Jet_hwPt[MAXJETS], l1Jet_hwEta[MAXJETS], l1Jet_hwPhi[MAXJETS], l1Jet_hwQual[MAXJETS];
  Float_t l1Jet_pt[MAXJETS], l1Jet_eta[MAXJETS], l1Jet_phi[MAXJETS];

  Int_t emcand_hwPt[144], emcand_hwPhi[144], emcand_hwEta[144], emcand_hwIso[144];
  Int_t iso3x3[144], isoCross[144], isoFourPoint[144], isoSingle[144], isoHolePunch[144];

  Int_t nPhoton;
  Float_t photon_pt[MAXJETS];
  Float_t photon_eta[MAXJETS];
  Float_t photon_phi[MAXJETS];
  Float_t cc4[MAXJETS];
  Float_t cr4[MAXJETS];
  Float_t ct4PtCut20[MAXJETS];
  Float_t trkSumPtHollowConeDR04[MAXJETS];
  Float_t hcalTowerSumEtConeDR04[MAXJETS];
  Float_t ecalRecHitSumEtConeDR04[MAXJETS];
  Float_t hadronicOverEm[MAXJETS];
  Float_t sigmaIetaIeta[MAXJETS];
  Float_t isEle[MAXJETS];
  Float_t sigmaIphiIphi[MAXJETS];
  Float_t swissCrx[MAXJETS];
  Float_t seedTime[MAXJETS];

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

  inTree->SetBranchAddress("nl1Jet",&nl1Jet);
  inTree->SetBranchAddress("l1Jet_hwPt",l1Jet_hwPt);
  inTree->SetBranchAddress("l1Jet_hwEta",l1Jet_hwEta);
  inTree->SetBranchAddress("l1Jet_hwPhi",l1Jet_hwPhi);
  inTree->SetBranchAddress("l1Jet_hwQual",l1Jet_hwQual);
  inTree->SetBranchAddress("l1Jet_pt",l1Jet_pt);
  inTree->SetBranchAddress("l1Jet_eta",l1Jet_eta);
  inTree->SetBranchAddress("l1Jet_phi",l1Jet_phi);

  inTree->SetBranchAddress("emcand_hwPt",emcand_hwPt);
  inTree->SetBranchAddress("emcand_hwPhi",emcand_hwPhi);
  inTree->SetBranchAddress("emcand_hwEta",emcand_hwEta);
  inTree->SetBranchAddress("emcand_hwIso",emcand_hwIso);
  inTree->SetBranchAddress("iso3x3",iso3x3);
  inTree->SetBranchAddress("isoCross",isoCross);
  inTree->SetBranchAddress("isoFourPoint",isoFourPoint);
  inTree->SetBranchAddress("isoSingle",isoSingle);
  inTree->SetBranchAddress("isoHolePunch",isoHolePunch);


  inTree->SetBranchAddress("nPhoton",&nPhoton);
  inTree->SetBranchAddress("photon_pt",photon_pt);
  inTree->SetBranchAddress("photon_eta",photon_eta);
  inTree->SetBranchAddress("photon_phi",photon_phi);

  inTree->SetBranchAddress("cc4",cc4);
  inTree->SetBranchAddress("cr4",cr4);
  inTree->SetBranchAddress("ct4PtCut20",ct4PtCut20);
  inTree->SetBranchAddress("trkSumPtHollowConeDR04",trkSumPtHollowConeDR04);
  inTree->SetBranchAddress("hcalTowerSumEtConeDR04",hcalTowerSumEtConeDR04);
  inTree->SetBranchAddress("ecalRecHitSumEtConeDR04",ecalRecHitSumEtConeDR04);
  inTree->SetBranchAddress("hadronicOverEm",hadronicOverEm);
  inTree->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
  inTree->SetBranchAddress("isEle",isEle);
  inTree->SetBranchAddress("sigmaIphiIphi",sigmaIphiIphi);
  inTree->SetBranchAddress("swissCrx",swissCrx);
  inTree->SetBranchAddress("seedTime",seedTime);

  // Int_t nGen;
  // Float_t gen_pt[MAXJETS], gen_eta[MAXJETS];//, gen_phi[MAXJETS];
  // Float_t gen_iso[MAXJETS];
  // Int_t gen_id[MAXJETS], gen_momId[MAXJETS];

  // inTree->SetBranchAddress("nGen",&nGen);
  // inTree->SetBranchAddress("gen_pt",gen_pt);
  // inTree->SetBranchAddress("gen_eta",gen_eta);
  // //inTree->SetBranchAddress("gen_phi",gen_phi);
  // inTree->SetBranchAddress("gen_iso",gen_iso);
  // inTree->SetBranchAddress("gen_id",gen_id);
  // inTree->SetBranchAddress("gen_momId",gen_momId);


  TFile *outFile = new TFile(Form("hist_photon_iso_data.root"),"RECREATE");
  outFile->cd();

  const int nBins = 50;
  const double maxPt = 100;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[2];
  fPt[0] = new TH1D("fPt_cen",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_periph");
  TH1D *accepted[THRESHOLDS][2];
  TH1D *isoDistribution = new TH1D("isoDistribution",";isolation energy (GeV)",100, -5, 45);
  TH1D *badIsoDistribution = (TH1D*)isoDistribution->Clone("badIsoDistribution");

  TH1D *jetSpectra = new TH1D("jetSpectra","l1 jet (GeV)",64, 0, 256);
  TH1D *badJetSpectra = (TH1D*)jetSpectra->Clone("badJetSpectra");
  TH1D *goodJetSpectra = (TH1D*)jetSpectra->Clone("goodJetSpectra");
  

  for(int i = 0; i < THRESHOLDS; ++i)
    for(int j = 0; j < 2; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",(int)L1_THRESHOLD[i],j),";offline p_{T}",nBins,0,maxPt);
    }

  TH2D *corr = new TH2D("corr",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);
  TH2D *matching = new TH2D("matching",";#Delta #eta;#Delta #phi",100,-5,5,100,0, TMath::Pi() );
  TH2D *matched_bad = (TH2D*)matching->Clone("matched_bad");

  Long64_t entries = inTree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    inTree->GetEntry(j);

    double maxl1eta = -999;
    double maxl1phi = -999;
    double maxl1pt = -1;
    double iso = -999;

    for(int i = 0; i < 144; ++i)
    {
      if( emcand_hwPt[i] > maxl1pt)
	if( ((0.5*iso3x3[i])-emcand_hwPt[i]) < 1.0 ) // 0.5 factor to convert to GeV, hacky
	{
	  maxl1pt = emcand_hwPt[i];
	  maxl1eta = getPhysicalEta(emcand_hwEta[i]);
	  maxl1phi = getPhysicalPhi(emcand_hwPhi[i]);
	  iso = (0.5*iso3x3[i])-emcand_hwPt[i];
	}
    }
    // double maxisopt = l1Egamma_pt[0];
    // double maxnonisopt = l1Egamma_pt[4];
    // double maxl1pt = std::max(maxisopt, maxnonisopt);
    //double maxl1pt = maxisopt;

    double maxfpt = -1;
    double maxfeta = -999;
    double maxfphi = -999;
    for(int i = 0; i < nPhoton; ++i)
    {
      // if(TMath::Abs(gen_momId[i]) <= 22)
      // if(gen_id[i] == 22)
      // if(gen_iso[i] < 5)
      // if(gen_pt[i] > maxfpt)
      if(photon_pt[i] > maxfpt)
      if((cc4[i] + cr4[i] + ct4PtCut20[i]) < 1)
      if(TMath::Abs(photon_eta[i]) < 1.4791)
      if(!isEle[i])
      if(TMath::Abs(seedTime[i])<3)
      if(swissCrx[i] < 0.9)
      if(sigmaIetaIeta[i] > 0.002)
      if(sigmaIphiIphi[i] > 0.002)
      if(hadronicOverEm[i] < 0.1)
      {
	maxfpt = photon_pt[i];
	maxfeta = photon_eta[i];
	maxfphi = photon_phi[i];
      }
    }
    //if(f_num > MAXJETS) std::cout << "TOO SMALL" << std::endl;
    l1Pt->Fill(maxl1pt);

    // if(goodEvent && maxl1pt < 8 && maxfpt > 30)
    // {
    //   std::cout << "***********" << std::endl;
    //   std::cout << "run lumi evt" << std::endl;
    //   std::cout << run << " " << lumi << " " << evt << std::endl;
    //   std::cout << "l1pt fpt" << std::endl;
    //   std::cout << maxl1pt << " " << maxfpt << std::endl;
    //   std::cout << "***********" << std::endl;
    // }

    if(goodEvent)
    {
      if(hiBin < 60)
	fPt[0]->Fill(maxfpt);
      else if (hiBin >= 100)
	fPt[1]->Fill(maxfpt);

      corr->Fill(maxfpt,maxl1pt);
      double diffphi = TMath::Abs(maxfphi-maxl1phi);
      if(diffphi > TMath::Pi())
	diffphi = (2*TMath::Pi()) - diffphi;
      double diffeta = maxfeta-maxl1eta;
      matching->Fill(diffeta, diffphi);
      isoDistribution->Fill(iso);
      jetSpectra->Fill(l1Jet_pt[0]);

      if(maxl1pt <= 5 && maxfpt > 40)
      {
	matched_bad->Fill(maxfeta-maxl1eta, diffphi);
	goodJetSpectra->Fill(l1Jet_pt[0]);
      }
      //double deltaR = TMath::Sqrt(diffphi*diffphi + diffeta*diffeta);
      //if(deltaR > 0.5)
      if(maxfpt < 15)
      {
	badIsoDistribution->Fill(iso);
	badJetSpectra->Fill(l1Jet_pt[0]);
      }

      for(int k = 0; k < THRESHOLDS; ++k)
      {
	if(maxl1pt>L1_THRESHOLD[k] || l1Jet_pt[0] > 36)
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
  matching->Write();
  matched_bad->Write();
  isoDistribution->Write();
  badIsoDistribution->Write();
  jetSpectra->Write();
  badJetSpectra->Write();
  goodJetSpectra->Write();
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
