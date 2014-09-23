#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>

#include <vector>
#include <iostream>
#include <map>

const int NREG = 396;

Long64_t makeKey(Int_t run, Int_t event){
  return (10000000000*(Long64_t)run + (Long64_t)event);
}

void makePUM0Table()
{
  const TString l1_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasHIanalyzer_withregions.root";
  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_event, l1_run;
  Int_t l1_num;
  Int_t region_hwPt[NREG], region_hwEta[NREG], region_hwPhi[NREG], region_tauVeto[NREG];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("nRegions",&l1_num);
  l1Tree->SetBranchAddress("region_hwPt",region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta",region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi",region_hwPhi);
  l1Tree->SetBranchAddress("region_tauVeto",region_tauVeto);

  const TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim.root";
  TFile *fFile = TFile::Open(forest_input);
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

  TFile *outFile = new TFile("HI_PUM0_evtsel_out.root","RECREATE");

  std::map<Long64_t, Long64_t> kmap;

  // choose loop over l1 tree first (smaller)
  //std::cout << "Begin making map." << std::endl;
  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);
    Long64_t key = makeKey(l1_run, l1_event);

    std::pair<Long64_t,Long64_t> p(key,j);
    kmap.insert(p);
  }
  //std::cout << "Finished making map." << std::endl;

  outFile->cd();

  TH1I *hists[22][18]; // [eta][pu bin], arbitrary value of 18 for # bins in pu
  for(int i = 0; i < 22; ++i)
    for(int j = 0; j < 18; ++j)
    {
      hists[i][j] = new TH1I(Form("hist_%d_%d",i,j),"", 1024,0,1024);
    }

  TH2I *centPUM = new TH2I("cenPUM","",200,0,200,396,0,396);
  int count = 0;

  Long64_t entries = fEvtTree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    //if(j % 10000 == 0)
    //  printf("%lld / %lld\n",j,entries);

    fEvtTree->GetEntry(j);
    Long64_t key = makeKey(f_run, f_evt);

    std::map<Long64_t,Long64_t>::const_iterator got = kmap.find(key);
    if(got == kmap.end() ) {
      continue;
    } else {
      l1Tree->GetEntry(got->second);
      kmap.erase(key);
      count++;

      fSkimTree->GetEntry(j);
      if((pcollisionEventSelection == 1) && (pHBHENoiseFilter == 1))
      {

	//int pubin = (int) ( (double)hiBin * (18.0/200.0));
	int PUM0 = 0;
	for(int i = 0; i < NREG; ++i)
	{
	  if(region_hwPt[i] > 0)
	    ++PUM0;
	}
	int pubin = PUM0/22;
	if(pubin == 18) pubin = 17; //special case for every region firing
	for(int i = 0; i < NREG; ++i)
	{
	  hists[region_hwEta[i]][pubin]->Fill(region_hwPt[i]);
	}
	centPUM->Fill(hiBin,PUM0);

      }
    }
  }

  std::cout << "cms.vdouble(";
  TH1D *hists_eta[22];
  for(int i = 0; i < 22; ++i)
  {
    hists_eta[i] = new TH1D(Form("hists_eta_%d",i),"",18,0,17);
    for(int j = 0; j < 18; ++j)
    {
      double Mean = hists[i][j]->GetMean();
      double MeanError = hists[i][j]->GetMeanError();
      hists_eta[i]->SetBinContent(j,Mean);
      hists_eta[i]->SetBinError(j,MeanError);
      std::cout << Mean*0.5;
      if(j != 17) std::cout << ", ";
    }
  }
  std::cout << ")" << std::endl;

  for(int i = 0; i < 22; ++i)
  {
    hists_eta[i]->Write();
    for(int j = 0; j < 18; ++j)
    {
      hists[i][j]->Write();
    }
  }
  centPUM->Write();


  //std::cout << "Matching entries: " << count << std::endl;

  lFile->Close();
  fFile->Close();
  outFile->Close();
}

int main()
{
  makePUM0Table();
  return 0;
}
