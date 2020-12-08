#include <vector>
#include <iostream>
#include <map>
#include <iterator>
#include <numeric>
#include <string>
#include <functional>

#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TF1.h"

using std::cout;
using std::endl;
using std::vector;

map <string, TH1*> histograms;

//============================ standalone functions===============================

double weight_DeepCSV(TH1F* eff, double pt){
	Double_t weight = 0;
	for(Int_t i_bin=1;i_bin<eff->GetNbinsX()+1;i_bin++){
		if(pt>eff->GetBinLowEdge(i_bin) && pt < eff->GetBinLowEdge(i_bin+1)){
			weight = eff->GetBinContent(i_bin);
        	if(TMath::IsNaN(weight)) weight = 1.0;
			continue;
		}
	}
	return weight;
}
//--------------------------------------------------------------------
double weight(TEfficiency* eff, double x){
  TH1 *total  = eff->GetCopyTotalHisto();
  TH1 *passed  = eff->GetCopyPassedHisto();

  Double_t xmin = 0;
  xmin = total->GetBinLowEdge(1);
  double binWidth = 0;
  binWidth = total->GetBinLowEdge(2) - xmin;  

  double xRem = 0;
  xRem = x - xmin;

  int rBin = 1;
  rBin = int(xRem/binWidth)+1;

  double weight = 0;
  weight = eff->GetEfficiency(rBin);

  delete total;
  delete passed;
  return weight;
}
//--------------------------------------------------------------------
double weight_error_down(TEfficiency* eff, double x){
	TH1 *total  = eff->GetCopyTotalHisto();

	Double_t xmin = 0;
  	xmin = total->GetBinLowEdge(1);
  	double binWidth = 0;
  	binWidth = total->GetBinLowEdge(2) - xmin;  

  	double xRem = 0;
  	xRem = x - xmin;

  	int rBin = 1;
  	rBin = int(xRem/binWidth)+1;

  	double weight_error_down = 0;
  	weight_error_down = eff->GetEfficiencyErrorLow(rBin);

 	delete total;
	return weight_error_down;
}
//--------------------------------------------------------------------
double weight_error_up(TEfficiency* eff, double x){
	TH1 *total  = eff->GetCopyTotalHisto();

	Double_t xmin = 0;
  	xmin = total->GetBinLowEdge(1);
  	double binWidth = 0;
  	binWidth = total->GetBinLowEdge(2) - xmin;  

  	double xRem = 0;
  	xRem = x - xmin;

  	int rBin = 1;
  	rBin = int(xRem/binWidth)+1;

  	double weight_error_up = 0;
  	weight_error_up = eff->GetEfficiencyErrorUp(rBin);

 	delete total;
	return weight_error_up;
}
//=========================Main program========================================================

void Uncertainties (TString name, TString MCSF, Double_t xsec){
	TString fname = "/fdata/hepx/store/user/delgado_andrea/Zprime94XNtuples/"+name+".root";
	TFile *f = TFile::Open(fname, "READ");
	TTree *t = (TTree*) f->Get("tree");

	//-------------------------------Open weights/SFs
	TFile *mc_eff = TFile::Open(MCSF);	
	TFile *w = TFile::Open("/fdata/hepx/store/user/delgado_andrea/94X/updatedNtuplesJune19/weightHistosAll.root");
	TFile *w45 = TFile::Open("/fdata/hepx/store/user/delgado_andrea/94X/updatedNtuplesJune19/weightHistosAllQuad45.root");
	TFile *trig = TFile::Open("ConfInterval.root", "READ");

  	TEfficiency *l1 = (TEfficiency*) w->Get("L1Quad30_caloT");
	TEfficiency *t1 = (TEfficiency*) w->Get("Quad30Calo_caloT");
  	TEfficiency *t2 = (TEfficiency*) w->Get("Di90Calo_caloT");
  	TEfficiency *t3 = (TEfficiency*) w->Get("Di90Di30LeadCSV3_min_PF");
  	TEfficiency *t4 = (TEfficiency*) w->Get("Quad30PF_PfT");  
  	TEfficiency *t5 = (TEfficiency*) w->Get("Di90PF_PfT");

  	TEfficiency *l145 = (TEfficiency*) w45->Get("L1Quad45_caloT");
  	TEfficiency *t145 = (TEfficiency*) w45->Get("Quad45Calo_caloT");
  	TEfficiency *t245 = (TEfficiency*) w45->Get("Quad45LeadCSV3_min_PF");
  	TEfficiency *t345 = (TEfficiency*) w45->Get("Quad45PF_PfT");  

  	TF1 *fit1      = (TF1*) trig->Get("central_L1Quad30_caloT");
  	TF1 *fit1_up   = (TF1*) trig->Get("upper_L1Quad30_caloT");
  	TF1 *fit1_down = (TF1*) trig->Get("lower_L1Quad30_caloT");
  	TF1 *fit1_unc   = (TF1*) trig->Get("uncertainty_L1Quad30_caloT");
  	TF1 *fit2      = (TF1*) trig->Get("central_Quad30Calo_caloT");
  	TF1 *fit2_up   = (TF1*) trig->Get("upper_Quad30Calo_caloT");
  	TF1 *fit2_down = (TF1*) trig->Get("lower_Quad30Calo_caloT");
  	TF1 *fit2_unc  = (TF1*) trig->Get("uncertainty_Quad30Calo_caloT");
  	TF1 *fit3      = (TF1*) trig->Get("central_Di90Calo_caloT");
  	TF1 *fit3_up   = (TF1*) trig->Get("upper_Di90Calo_caloT");
  	TF1 *fit3_down = (TF1*) trig->Get("lower_Di90Calo_caloT");
  	TF1 *fit3_unc  = (TF1*) trig->Get("uncertainty_Di90Calo_caloT");
  	TF1 *fit4      = (TF1*) trig->Get("central_Quad30PF_PfT");
  	TF1 *fit4_up   = (TF1*) trig->Get("upper_Quad30PF_PfT");
  	TF1 *fit4_down = (TF1*) trig->Get("lower_Quad30PF_PfT");
  	TF1 *fit4_unc  = (TF1*) trig->Get("uncertainty_Quad30PF_PfT");
  	TF1 *fit5      = (TF1*) trig->Get("central_Di90PF_PfT");
  	TF1 *fit5_up   = (TF1*) trig->Get("upper_Di90PF_PfT");
  	TF1 *fit5_down = (TF1*) trig->Get("lower_Di90PF_PfT");
  	TF1 *fit5_unc  = (TF1*) trig->Get("uncertainty_Di90PF_PfT");
  	TF1 *fitcsv    = (TF1*) trig->Get("central_Di90Di30LeadCSV3_min_PF");
  	TF1 *fitcsv_up = (TF1*) trig->Get("upper_Di90Di30LeadCSV3_min_PF");
  	TF1 *fitcsv_down = (TF1*) trig->Get("lower_Di90Di30LeadCSV3_min_PF");
  	TF1 *fitcsv_unc  = (TF1*) trig->Get("uncertainty_Di90Di30LeadCSV3_min_PF");


   	if (t == 0 ){
		std::cout << "No tree found!" <<std::endl;
		return;
   	}

   	Double_t lumi = 35.9;

	Float_t GenJet_phi[50];
	Float_t GenJet_eta[50];
	Float_t GenJet_pt[50];
	Float_t GenJet_mass[50];	
	Float_t Jet_rawPt[50];
	Float_t Jet_pt[50];		
	Float_t Jet_eta[50];		
	Float_t Jet_phi[50];		
	Float_t Jet_mass[50];
	Int_t Jet_hadronFlavour[50];	
	Float_t Jet_btagCMVA[50];				
	Float_t Jet_btagCSV[50];		
	Float_t Jet_btagDeepCSVb[50];			
	Float_t Jet_btagDeepCSVbb[50];		
 	Float_t Jet_chHEF[50];
  	Float_t Jet_neHEF[50];
  	Float_t Jet_chEmEF[50];
  	Float_t Jet_neEmEF[50];
  	Float_t Jet_corr[50];
  	Float_t Jet_corr_JECUp[50];
  	Float_t Jet_corr_JECDown[50];
  	Float_t Jet_corr_TotalUp[50];
  	Float_t Jet_corr_TotalDown[50];
  	Float_t Jet_corr_JER[50];
  	Float_t Jet_corr_JERUp[50];
  	Float_t Jet_corr_JERDown[50];  	  
  	Int_t Jet_numberOfDaughters[50];
  	Float_t Jet_muEF[50];
  	Int_t Jet_chMult[50];
  	Int_t Jet_nhMult[50];		
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v;
	Float_t Jet_btagDeepCSVL_SF_up[50];
	Float_t Jet_btagDeepCSVL_SF_down[50];
	Float_t Jet_btagDeepCSVL_SF[50];
	Float_t Jet_btagDeepCSVM_SF[50];
	Float_t Jet_btagDeepCSVM_SF_up[50];
	Float_t Jet_btagDeepCSVM_SF_down[50];
	Float_t Jet_btagDeepCSVT_SF[50];
	Float_t Jet_btagDeepCSVT_SF_up[50];
	Float_t Jet_btagDeepCSVT_SF_down[50];
	Float_t trgObjects_caloJets_pt[50];
	Float_t trgObjects_pfJets_pt[50];
	Float_t LHE_weights_pdf_wgt[150];

	Int_t nJet;
	Int_t nGenJet;	
	Int_t ntrgObjects_caloJets;
	Int_t ntrgObjects_pfJets;
	Float_t puWeight;
	Float_t puWeightUp;
	Float_t puWeightDown;

	t->SetBranchAddress("LHE_weights_pdf_wgt", LHE_weights_pdf_wgt);
	t->SetBranchAddress("puWeight", &puWeight);
	t->SetBranchAddress("puWeightUp", &puWeightUp);
	t->SetBranchAddress("puWeightDown", &puWeightDown);
	t->SetBranchAddress("ntrgObjects_caloJets", &ntrgObjects_caloJets);
	t->SetBranchAddress("ntrgObjects_pfJets", &ntrgObjects_pfJets);
	t->SetBranchAddress("trgObjects_caloJets_pt", trgObjects_caloJets_pt);
	t->SetBranchAddress("trgObjects_pfJets_pt", trgObjects_pfJets_pt);
	t->SetBranchAddress("GenJet_eta",GenJet_eta);			
	t->SetBranchAddress("GenJet_phi",GenJet_phi);			
	t->SetBranchAddress("GenJet_pt",GenJet_pt);			
	t->SetBranchAddress("GenJet_mass",GenJet_mass);	
	t->SetBranchAddress("Jet_rawPt", Jet_rawPt);		
	t->SetBranchAddress("Jet_pt", Jet_pt);	
	t->SetBranchAddress("Jet_eta", Jet_eta);	
	t->SetBranchAddress("Jet_phi", Jet_phi);	
	t->SetBranchAddress("Jet_mass", Jet_mass);	
	t->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
	t->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);			
	t->SetBranchAddress("Jet_btagDeepCSVb", Jet_btagDeepCSVb);				
	t->SetBranchAddress("Jet_btagDeepCSVbb", Jet_btagDeepCSVbb);					
	t->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA);
	t->SetBranchAddress("Jet_chHEF", Jet_chHEF);
  	t->SetBranchAddress("Jet_neHEF", Jet_neHEF);
  	t->SetBranchAddress("Jet_chEmEF", Jet_chEmEF);
  	t->SetBranchAddress("Jet_neEmEF", Jet_neEmEF);
  	t->SetBranchAddress("Jet_numberOfDaughters", &Jet_numberOfDaughters);
  	t->SetBranchAddress("Jet_muEF", Jet_muEF);
  	t->SetBranchAddress("Jet_chMult",&Jet_chMult);
  	t->SetBranchAddress("Jet_nhMult", &Jet_nhMult);
	t->SetBranchAddress("nJet",&nJet);			
	t->SetBranchAddress("nGenJet",&nGenJet);			
	t->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v);
	t->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v);	
	t->SetBranchAddress("Jet_btagDeepCSVT_SF", Jet_btagDeepCSVT_SF);
	t->SetBranchAddress("Jet_btagDeepCSVT_SF_up", Jet_btagDeepCSVT_SF_up);
	t->SetBranchAddress("Jet_btagDeepCSVT_SF_down", Jet_btagDeepCSVT_SF_down);
	t->SetBranchAddress("Jet_btagDeepCSVM_SF", Jet_btagDeepCSVM_SF);
	t->SetBranchAddress("Jet_btagDeepCSVM_SF_up", Jet_btagDeepCSVM_SF_up);
	t->SetBranchAddress("Jet_btagDeepCSVM_SF_down", Jet_btagDeepCSVM_SF_down);
	t->SetBranchAddress("Jet_btagDeepCSVL_SF", Jet_btagDeepCSVL_SF);	
	t->SetBranchAddress("Jet_btagDeepCSVL_SF_up", Jet_btagDeepCSVL_SF_up);		
	t->SetBranchAddress("Jet_btagDeepCSVL_SF_down", Jet_btagDeepCSVL_SF_down);
	t->SetBranchAddress("Jet_corr", Jet_corr);
	t->SetBranchAddress("Jet_corr_JECUp",Jet_corr_JECUp);
	t->SetBranchAddress("Jet_corr_JECDown",Jet_corr_JECDown);
	t->SetBranchAddress("Jet_corr_TotalDown", Jet_corr_TotalDown);
	t->SetBranchAddress("Jet_corr_TotalUp", Jet_corr_TotalUp);	
	t->SetBranchAddress("Jet_corr_JER", Jet_corr_JER);
	t->SetBranchAddress("Jet_corr_JERUp", Jet_corr_JERUp);
	t->SetBranchAddress("Jet_corr_JERDown", Jet_corr_JERDown);	

	histograms["mjjsr"] = new TH1F("mjjsr", "mjj from two leading jets in event for signal region ", 120,0, 1200);	
	histograms["mjjcr1"] = new TH1F("mjjcr1", "mjj from two leading jets in event for control region 1", 120,0, 1200);	
	histograms["mjjcr2"] = new TH1F("mjjcr2", "mjj from two leading jets in event for control region 2", 120,0, 1200);		
	histograms["mjjcr3"] = new TH1F("mjjcr3", "mjj from two leading jets in event for control region 3", 120,0, 1200);

	//-----------------------histograms for DeepCSV SF variation 
	histograms["mjjsr_DeepCSVUp"] = new TH1F("mjjsr_DeepCSVUp", "mjj from two leading jets in event for signal region after DeepCSVReWeight + SF up", 120,0, 1200);	
	histograms["mjjcr1_DeepCSVUp"] = new TH1F("mjjcr1_DeepCSVUp", "mjj from two leading jets in event for control region 1 after DeepCSVReWeight + SF up", 120,0, 1200);	
	histograms["mjjcr2_DeepCSVUp"] = new TH1F("mjjcr2_DeepCSVUp", "mjj from two leading jets in event for control region 2 after DeepCSVReWeight + SF up", 120,0, 1200);		
	histograms["mjjcr3_DeepCSVUp"] = new TH1F("mjjcr3_DeepCSVUp", "mjj from two leading jets in event for control region 3 after DeepCSVReWeight + SF up", 120,0, 1200);			
	
	histograms["mjjsr_DeepCSVDown"] = new TH1F("mjjsr_DeepCSVDown", "mjj from two leading jets in event for signal region after DeepCSVReWeight + SF down", 120,0, 1200);	
	histograms["mjjcr1_DeepCSVDown"] = new TH1F("mjjcr1_DeepCSVDown", "mjj from two leading jets in event for control region 1 after DeepCSVReWeight + SF down", 120,0, 1200);	
	histograms["mjjcr2_DeepCSVDown"] = new TH1F("mjjcr2_DeepCSVDown", "mjj from two leading jets in event for control region 2 after DeepCSVReWeight + SF down", 120,0, 1200);		
	histograms["mjjcr3_DeepCSVDown"] = new TH1F("mjjcr3_DeepCSVDown", "mjj from two leading jets in event for control region 3 after DeepCSVReWeight + SF down", 120,0, 1200);			
	//---------------------------histograms for trigger weight variations
	histograms["mjjsr_TriggerUp"] = new TH1F("mjjsr_TriggerUp", "mjj from two leading jets in event for signal region using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr1_TriggerUp"] = new TH1F("mjjcr1_TriggerUp", "mjj from two leading jets in event for control region 1 using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr2_TriggerUp"] = new TH1F("mjjcr2_TriggerUp", "mjj from two leading jets in event for control region 2 using trigger weights + shift up", 120,0, 1200);		
	histograms["mjjcr3_TriggerUp"] = new TH1F("mjjcr3_TriggerUp", "mjj from two leading jets in event for control region 3 using trigger weights + shift up", 120,0, 1200);			

	histograms["mjjsr_TriggerDown"] = new TH1F("mjjsr_TriggerDown", "mjj from two leading jets in event for signal region using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr1_TriggerDown"] = new TH1F("mjjcr1_TriggerDown", "mjj from two leading jets in event for control region 1 using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr2_TriggerDown"] = new TH1F("mjjcr2_TriggerDown", "mjj from two leading jets in event for control region 2 using trigger weights + shift down", 120,0, 1200);		
	histograms["mjjcr3_TriggerDown"] = new TH1F("mjjcr3_TriggerDown", "mjj from two leading jets in event for control region 3 using trigger weights + shift down", 120,0, 1200);			
	//--------------------------histograms for JES variations
	histograms["mjjsr_jesUp"] = new TH1F("mjjsr_jesUp", "mjj from two leading jets in event for signal region +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1_jesUp"] = new TH1F("mjjcr1_jesUp", "mjj from two leading jets in event for control region 1 +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2_jesUp"] = new TH1F("mjjcr2_jesUp", "mjj from two leading jets in event for control region 2 +1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3_jesUp"] = new TH1F("mjjcr3_jesUp", "mjj from two leading jets in event for control region 3 +1#sigma JES", 120,0, 1200);			
	
	histograms["mjjsr_jesDown"] = new TH1F("mjjsr_jesDown", "mjj from two leading jets in event for signal region -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1_jesDown"] = new TH1F("mjjcr1_jesDown", "mjj from two leading jets in event for control region 1 -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2_jesDown"] = new TH1F("mjjcr2_jesDown", "mjj from two leading jets in event for control region 2 -1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3_jesDown"] = new TH1F("mjjcr3_jesDown", "mjj from two leading jets in event for control region 3 -1#sigma JES", 120,0, 1200);			
	//--------------------------histograms for JER variations
	histograms["mjjsr_jerUp"] = new TH1F("mjjsr_jerUp", "mjj from two leading jets in event for signal region +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1_jerUp"] = new TH1F("mjjcr1_jerUp", "mjj from two leading jets in event for control region 1 +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2_jerUp"] = new TH1F("mjjcr2_jerUp", "mjj from two leading jets in event for control region 2 +1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3_jerUp"] = new TH1F("mjjcr3_jerUp", "mjj from two leading jets in event for control region 3 +1#sigma JER", 120,0, 1200);			

	histograms["mjjsr_jerDown"] = new TH1F("mjjsr_jerDown", "mjj from two leading jets in event for signal region -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1_jerDown"] = new TH1F("mjjcr1_jerDown", "mjj from two leading jets in event for control region 1 -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2_jerDown"] = new TH1F("mjjcr2_jerDown", "mjj from two leading jets in event for control region 2 -1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3_jerDown"] = new TH1F("mjjcr3_jerDown", "mjj from two leading jets in event for control region 3 -1#sigma JER", 120,0, 1200);			
	//--------------------------histograms for PU weight variations
	histograms["mjjsr_PUUp"] = new TH1F("mjjsr_PUUp", "mjj from two leading jets in event for signal region +4.6p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1_PUUp"] = new TH1F("mjjcr1_PUUp", "mjj from two leading jets in event for control region 1 +4.6p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2_PUUp"] = new TH1F("mjjcr2_PUUp", "mjj from two leading jets in event for control region 2 +4.6p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3_PUUp"] = new TH1F("mjjcr3_PUUp", "mjj from two leading jets in event for control region 3 +4.6p minbias xsec", 120,0, 1200);			

	histograms["mjjsr_PUDown"] = new TH1F("mjjsr_PUDown", "mjj from two leading jets in event for signal region -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1_PUDown"] = new TH1F("mjjcr1_PUDown", "mjj from two leading jets in event for control region 1 -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2_PUDown"] = new TH1F("mjjcr2_PUDown", "mjj from two leading jets in event for control region 2 -7p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3_PUDown"] = new TH1F("mjjcr3_PUDown", "mjj from two leading jets in event for control region 3 -7p minbias xsec", 120,0, 1200);			
	//--------------------------histograms for PDF variation
	histograms["mjjsr_PDFUp"] = new TH1F("mjjsr_PDFUp", "mjj from two leading jets in event for signal region + PDF up", 120,0, 1200);	
	histograms["mjjsr_PDFDown"] = new TH1F("mjjsr_PDFDown", "mjj from two leading jets in event for signal region + PDF down", 120,0, 1200);	
	//--------------------------histograms for alpha variation
	histograms["mjjsr_alphaUp"] = new TH1F("mjjsr_alphaUp", "mjj from two leading jets in event for signal region + alpha up", 120,0, 1200);	
	histograms["mjjsr_alphaDown"] = new TH1F("mjjsr_alphaDown", "mjj from two leading jets in event for signal region + alpha down", 120,0, 1200);

//================================================================================================================
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << t->GetEntries() << endl;

	TH1F *Fl0_WPT = (TH1F*) mc_eff->Get("Fl0_WPT"); 
	TH1F *Fl0_WPM = (TH1F*) mc_eff->Get("Fl0_WPM");
	TH1F *Fl0_WPL = (TH1F*) mc_eff->Get("Fl0_WPL");
	TH1F *Fl4_WPT = (TH1F*) mc_eff->Get("Fl4_WPT"); 
	TH1F *Fl4_WPM = (TH1F*) mc_eff->Get("Fl4_WPM");
	TH1F *Fl4_WPL = (TH1F*) mc_eff->Get("Fl4_WPL");
	TH1F *Fl5_WPT = (TH1F*) mc_eff->Get("Fl5_WPT"); 
	TH1F *Fl5_WPM = (TH1F*) mc_eff->Get("Fl5_WPM");
	TH1F *Fl5_WPL = (TH1F*) mc_eff->Get("Fl5_WPL");

	Double_t nevts = t->GetEntries();
//=======================================Loop over entries===========================================================
	for (Int_t indx=0;indx<nevts;indx++){
	//for (Int_t indx=0;indx<10;indx++){
		t->GetEntry(indx);

		//------------------------------------------PDF Weigt calculation
		Double_t pdf_Weight = 0.0, temp_pdf = 0.0, pdf_Weight_up = 0.0, pdf_Weight_down = 0.0;
		for(unsigned int i = 1; i<101; i++){
			temp_pdf += (LHE_weights_pdf_wgt[i]-LHE_weights_pdf_wgt[0])*(LHE_weights_pdf_wgt[i]-LHE_weights_pdf_wgt[0]);
		}	
		pdf_Weight = TMath::Sqrt(temp_pdf);	
		pdf_Weight_up = 1.0 + pdf_Weight;
		pdf_Weight_down = 1.0 - pdf_Weight;
		//-----------------------------------------alpha_s variation

		Double_t alpha = 0.0, alpha_up = 0.0, alpha_down = 0.0;
		alpha = (LHE_weights_pdf_wgt[102] - LHE_weights_pdf_wgt[101])/2;
		alpha_up = 1.0 + alpha;
		alpha_down = 1.0 - alpha;


		bool passT4 = false, passT3 = false, passT2 = false, FiT = false, FiM = false, FiL = false, SeT = false,
		SeM = false, SeL = false, ThT = false, ThM = false, ThL = false, FoT = false, FoM = false, FoL = false;

		Int_t nJet30 = 0;

		float MaxBtag=0., MinBtag=1., PreMinBtag=1.;
    	vector<TLorentzVector> PFjets;
		vector<float> CSV;

		Double_t w_DeepCSV_MC = 0, w_DeepCSV_Data = 0, w_DeepCSV_Data_up = 0, w_DeepCSV_Data_down = 0;

		std::vector<double> jet_DeepCSVweight_MC;	
		std::vector<double> jet_DeepCSVweight_Data;	
		std::vector<double> jet_DeepCSVweight_Data_up;	
		std::vector<double> jet_DeepCSVweight_Data_down;

		//Loop over jets
		for (int j=0; j<nJet; j++){

			if(Jet_pt[j]<30 || fabs(Jet_eta[j])>2.4) continue;

         	TLorentzVector jet;
      		const float DeepCSV=Jet_btagDeepCSVb[j]+Jet_btagDeepCSVbb[j];
      		jet.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],Jet_mass[j]);
         	Double_t temp_weight = 0, temp_weight_data = 0, temp_weight_data_up = 0, temp_weight_data_down = 0;//, DeepCSV = 0;
         	//DeepCSV = Jet_btagDeepCSVb[j]+Jet_btagDeepCSVbb[j];

         	if(DeepCSV>0.8958){
         		if(Jet_hadronFlavour[j]==0){
         			temp_weight = weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==4){
         			temp_weight = weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==5){
         			temp_weight = weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);         			         			
         		}
         		jet_DeepCSVweight_MC.push_back(temp_weight);
         		jet_DeepCSVweight_Data.push_back(temp_weight_data);
         		jet_DeepCSVweight_Data_up.push_back(temp_weight_data_up);
         		jet_DeepCSVweight_Data_down.push_back(temp_weight_data_down);         		         		
         	}
         	else if(DeepCSV>0.6324 && DeepCSV < 0.8958){
         		if(Jet_hadronFlavour[j]==0){
         			temp_weight = weight_DeepCSV(Fl0_WPM, Jet_pt[j]) - weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl0_WPT, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==4){
         			temp_weight = weight_DeepCSV(Fl4_WPM, Jet_pt[j]) - weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl4_WPT, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==5){
         			temp_weight = weight_DeepCSV(Fl5_WPM, Jet_pt[j]) - weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         		    temp_weight_data = Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         		    temp_weight_data_up = Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_up[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);
         		    temp_weight_data_down = Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]) - Jet_btagDeepCSVT_SF_down[j]*weight_DeepCSV(Fl5_WPT, Jet_pt[j]);         		             		    
         		}
         		jet_DeepCSVweight_MC.push_back(temp_weight);
         		jet_DeepCSVweight_Data.push_back(temp_weight_data);
         		jet_DeepCSVweight_Data_up.push_back(temp_weight_data_up);
         		jet_DeepCSVweight_Data_down.push_back(temp_weight_data_down);         		         		
         	}
         	else if(DeepCSV>0.2219 && DeepCSV < 0.6324){
         		if(Jet_hadronFlavour[j]==0){
         			temp_weight = weight_DeepCSV(Fl0_WPL, Jet_pt[j]) - weight_DeepCSV(Fl0_WPM, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]);
         			temp_weight_data_up = Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]);
         			temp_weight_data_down = Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl0_WPM, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==4){
         			temp_weight = weight_DeepCSV(Fl4_WPL, Jet_pt[j]) - weight_DeepCSV(Fl4_WPM, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]);         			
         			temp_weight_data_up = Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]);         			
         			temp_weight_data_down = Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl4_WPM, Jet_pt[j]);         			         			         			
         		}
         		if(Jet_hadronFlavour[j]==5){
         			temp_weight = weight_DeepCSV(Fl5_WPL, Jet_pt[j]) - weight_DeepCSV(Fl5_WPM, Jet_pt[j]);
         			temp_weight_data = Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]);         		
         			temp_weight_data_up = Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_up[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]);         		
         			temp_weight_data_down = Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]) - Jet_btagDeepCSVM_SF_down[j]*weight_DeepCSV(Fl5_WPM, Jet_pt[j]);         		         			         			
         		}
         		jet_DeepCSVweight_MC.push_back(temp_weight);
         		jet_DeepCSVweight_Data.push_back(temp_weight_data);
         		jet_DeepCSVweight_Data_up.push_back(temp_weight_data_up);
         		jet_DeepCSVweight_Data_down.push_back(temp_weight_data_down);         		         		
         	}
         	else if(DeepCSV<0.2219){
         		if(Jet_hadronFlavour[j]==0){
         			temp_weight_data = 1.0 - Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]);
         			temp_weight_data_up = 1.0 - Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]);
         			temp_weight_data_down = 1.0 - Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl0_WPL, Jet_pt[j]);         			         			
         			temp_weight = 1.0 - weight_DeepCSV(Fl0_WPL, Jet_pt[j]);
         		}
         		if(Jet_hadronFlavour[j]==4){
         			temp_weight = 1.0 - weight_DeepCSV(Fl4_WPL, Jet_pt[j]);
         			temp_weight_data = 1.0 - Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]);
         			temp_weight_data_up = 1.0 - Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]);
         			temp_weight_data_down = 1.0 - Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl4_WPL, Jet_pt[j]);         			         			
         		}
         		if(Jet_hadronFlavour[j]==5){
         			temp_weight = 1.0 - weight_DeepCSV(Fl5_WPL, Jet_pt[j]);
         			temp_weight_data = 1.0 - Jet_btagDeepCSVL_SF[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]);         			
         			temp_weight_data_up = 1.0 - Jet_btagDeepCSVL_SF_up[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]);         			
         			temp_weight_data_down = 1.0 - Jet_btagDeepCSVL_SF_down[j]*weight_DeepCSV(Fl5_WPL, Jet_pt[j]);         			         			         			
         		}
         		jet_DeepCSVweight_MC.push_back(temp_weight);
         		jet_DeepCSVweight_Data.push_back(temp_weight_data);
         		jet_DeepCSVweight_Data_up.push_back(temp_weight_data_up);
         		jet_DeepCSVweight_Data_down.push_back(temp_weight_data_down);         		         		
         	}

         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9){
	               			nJet30++;
	               		    if(PFjets.size()<4 && DeepCSV>MaxBtag) MaxBtag=DeepCSV;                     
              				PFjets.push_back(jet);
              				CSV.push_back(DeepCSV);
	               		}
               		}
               		else{
               			nJet30++;
            			PFjets.push_back(jet);CSV.push_back(DeepCSV);
            		}
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30++;
    	 
		} //Loop over jets

		Double_t multi = 1.0, multi2 = 1.0, multi2_up = 1.0, multi2_down = 1.0;

		for (unsigned int i=0; i<jet_DeepCSVweight_MC.size(); i++){
        	if(TMath::IsNaN(jet_DeepCSVweight_MC[i])) jet_DeepCSVweight_MC[i] = 1.0;
			multi *= jet_DeepCSVweight_MC[i];
		}
		for (unsigned int i=0; i<jet_DeepCSVweight_Data.size(); i++){
        	if(TMath::IsNaN(jet_DeepCSVweight_Data[i])) jet_DeepCSVweight_Data[i] = 1.0;
			multi2 *= jet_DeepCSVweight_Data[i];
		}
		for (unsigned int i=0; i<jet_DeepCSVweight_Data_up.size(); i++){
        	if(TMath::IsNaN(jet_DeepCSVweight_Data_up[i])) jet_DeepCSVweight_Data_up[i] = 1.0;
			multi2_up *= jet_DeepCSVweight_Data_up[i];
		}
		for (unsigned int i=0; i<jet_DeepCSVweight_Data_down.size(); i++){
        	if(TMath::IsNaN(jet_DeepCSVweight_Data_down[i])) jet_DeepCSVweight_Data_down[i] = 1.0;
			multi2_down *= jet_DeepCSVweight_Data_down[i];
		}

		Double_t DeepCSVEventWeight = 0, DeepCSVEventWeight_up = 0, DeepCSVEventWeight_down = 0;				
		DeepCSVEventWeight = multi2/multi;
		DeepCSVEventWeight_up = multi2_up/multi;
		DeepCSVEventWeight_down = multi2_down/multi;
		if(TMath::IsNaN(DeepCSVEventWeight)) DeepCSVEventWeight = 1.0;
		if(TMath::IsNaN(DeepCSVEventWeight_up)) DeepCSVEventWeight_up = 1.0;
		if(TMath::IsNaN(DeepCSVEventWeight_down)) DeepCSVEventWeight_down = 1.0;

		//====================b-tagging selection=====================================
		if(nJet30>3){
        	if(Jet_pt[0]>30 && fabs(Jet_eta[0])<2.4){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL=true;
         	}
        	if(Jet_pt[1]>30 && fabs(Jet_eta[1])<2.4){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL=true;
         	}

			if(Jet_pt[2]>30 && fabs(Jet_eta[2])<2.4){
            	if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT=true;
            	else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM=true;
            	else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL=true;
			}
			if(Jet_pt[3]>30 && fabs(Jet_eta[3])<2.4){
            	if(Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.8958) FoT=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.6324) FoM=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.2219) FoL=true;
      		}	
      	}
		if(FiT && SeT){
			passT2=true;
				if(ThT){
					passT3=true;
					if(FoT) passT4=true;
					else if(FoM) passT4=true;					
					else if(FoL) passT4=true;
					else if(!FoL && !FoM && !FoT) passT4=true;											
				} 
				else if(ThM){
					passT3=true;
					if(FoT) passT4=true;
					else if(FoM) passT4=true;
				} 
				else if(ThL){
					passT3=true;
					if(FoT) passT4=true;				
				}
		}					
		else if(FiT && SeM){
				passT2=true;
				if(ThT){
					passT3=true;
					if(FoT) passT4=true;
					else if(FoM) passT4=true;					
				}
				else if(ThM){
					passT3=true;
					if(FoT) passT4=true;
				}				
			}					
		else if(FiM && SeT){
				passT2=true;
				if(ThT){
					passT3=true;
					if(FoT) passT4=true;
					else if(FoM) passT4=true;					
				}
				else if(ThM){
					passT3=true;				
					if(FoT) passT4=true;
				}
		} 
		else if(FiM && SeM){
				passT2=true;
				if(ThT){
					passT3=true;	
					if(FoT) passT4=true;
				}
		}
//========================JES/JER variations=============================================
	//-------------------------JES Up 

		bool passT4_Up = false, passT3_Up = false, passT2_Up = false, FiT_Up = false, FiM_Up = false,
		FiL_Up = false, SeT_Up = false, SeM_Up = false, SeL_Up = false, ThT_Up = false, ThM_Up = false,
		ThL_Up = false, FoT_Up = false, FoM_Up = false, FoL_Up = false;
		
		Int_t nJet30_Up = 0;
		
		Double_t Jet_pt0_Up = Jet_corr_TotalUp[0] * Jet_rawPt[0] *Jet_corr_JER[0];
		Double_t Jet_pt1_Up = Jet_corr_TotalUp[1] * Jet_rawPt[1] *Jet_corr_JER[1];
		Double_t Jet_pt2_Up = Jet_corr_TotalUp[2] * Jet_rawPt[2] *Jet_corr_JER[2];
		Double_t Jet_pt3_Up = Jet_corr_TotalUp[3] * Jet_rawPt[3] *Jet_corr_JER[3];


		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_up = 0;
			Jet_pt_up = Jet_corr_TotalUp[j]* Jet_rawPt[j]*Jet_corr_JER[j];
			if(Jet_pt_up<30 || fabs(Jet_eta[j])>2.4) continue;
         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9) nJet30_Up++;
               		}
               		else nJet30_Up++;
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30_Up++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30_Up++;
		}
		if(nJet30_Up>3){
        	if(Jet_pt0_Up>30 && fabs(Jet_eta[0])<2.4){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_Up=true;
         	}
        	if(Jet_pt1_Up>30 && fabs(Jet_eta[1])<2.4){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_Up=true;
         	}
			if(Jet_pt2_Up>30 && fabs(Jet_eta[2])<2.4){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_Up=true;
			}
			if(Jet_pt3_Up>30 && fabs(Jet_eta[3])<2.4){
            	if(Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.8958) FoT_Up=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.6324) FoM_Up=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.2219) FoL_Up=true;
      		}
      	}	
		if(FiT_Up && SeT_Up){
			passT2_Up=true;
			if(ThT_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;
				else if(FoM_Up) passT4_Up=true;					
				else if(FoL_Up) passT4_Up=true;
				else if(!FoL_Up && !FoM_Up && !FoT_Up) passT4_Up=true;											
			} 
			else if(ThM_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;
				else if(FoM_Up) passT4_Up=true;
			} 
			else if(ThL_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;				
			}
		}					
		else if(FiT_Up && SeM_Up){
			passT2_Up=true;
			if(ThT_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;
				else if(FoM_Up) passT4_Up=true;					
			}
			else if(ThM_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;
			}				
		}					
		else if(FiM_Up && SeT_Up){
			passT2_Up=true;
			if(ThT_Up){
				passT3_Up=true;
				if(FoT_Up) passT4_Up=true;
				else if(FoM_Up) passT4_Up=true;					
			}
			else if(ThM_Up){
				passT3_Up=true;				
				if(FoT_Up) passT4_Up=true;
			}
		} 
		else if(FiM_Up && SeM_Up){
			passT2_Up=true;
			if(ThT_Up){
				passT3_Up=true;	
				if(FoT_Up) passT4_Up=true;
			}
		}

	//---------------------------JEC Down 

		bool passT4_Down = false, passT3_Down = false, passT2_Down = false, FiT_Down = false, FiM_Down = false,
		FiL_Down = false, SeT_Down = false, SeM_Down = false, SeL_Down = false, ThT_Down = false,
		ThM_Down = false, ThL_Down = false, FoT_Down = false, FoM_Down = false, FoL_Down = false;
		
		Int_t nJet30_Down = 0;
		
		Double_t Jet_pt0_Down = Jet_corr_TotalDown[0] * Jet_rawPt[0]*Jet_corr_JER[0];
		Double_t Jet_pt1_Down = Jet_corr_TotalDown[1] * Jet_rawPt[1]*Jet_corr_JER[1] ;
		Double_t Jet_pt2_Down = Jet_corr_TotalDown[2] * Jet_rawPt[2]*Jet_corr_JER[2];
		Double_t Jet_pt3_Down = Jet_corr_TotalDown[3] * Jet_rawPt[3]*Jet_corr_JER[3];

		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_Down = 0;
			Jet_pt_Down = Jet_corr_TotalDown[j]* Jet_rawPt[j]*Jet_corr_JER[j];
			if(Jet_pt_Down<30 || fabs(Jet_eta[j])>2.4) continue;
         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9) nJet30_Down++;
               		}
               		else nJet30_Down++;
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30_Down++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30_Down++;
		}
		if(nJet30_Down>3){
        	if(Jet_pt0_Down>30 && fabs(Jet_eta[0])<2.4){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_Down=true;
         	}
        	if(Jet_pt1_Down>30 && fabs(Jet_eta[1])<2.4){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_Down=true;
         	}
			if(Jet_pt2_Down>30 && fabs(Jet_eta[2])<2.4){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_Down=true;
			}
			if(Jet_pt3_Down>30 && fabs(Jet_eta[3])<2.4){
            	if(Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.8958) FoT_Down=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.6324) FoM_Down=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.2219) FoL_Down=true;
      		}
      	}	
		if(FiT_Down && SeT_Down){
			passT2_Down=true;
			if(ThT_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;
				else if(FoM_Down) passT4_Down=true;					
				else if(FoL_Down) passT4_Down=true;
				else if(!FoL_Down && !FoM_Down && !FoT_Down) passT4_Down=true;											
			} 
			else if(ThM_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;
				else if(FoM_Down) passT4_Down=true;
			} 
			else if(ThL_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;				
			}
		}					
		else if(FiT_Down && SeM_Down){
			passT2_Down=true;
			if(ThT_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;
				else if(FoM_Down) passT4_Down=true;					
			}
			else if(ThM_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;
			}				
		}					
		else if(FiM_Down && SeT_Down){
			passT2_Down=true;
			if(ThT_Down){
				passT3_Down=true;
				if(FoT_Down) passT4_Down=true;
				else if(FoM_Down) passT4_Down=true;					
			}
			else if(ThM_Down){
				passT3_Down=true;				
				if(FoT_Down) passT4_Down=true;
			}
		} 
		else if(FiM_Down && SeM_Down){
			passT2_Down=true;
			if(ThT_Down){
				passT3_Down=true;	
				if(FoT_Down) passT4_Down=true;
			}
		}
	//--------------------------------- JER Up 

		bool passT4_JER_Up = false, passT3_JER_Up = false, passT2_JER_Up = false, FiT_JER_Up = false,
		FiM_JER_Up = false, FiL_JER_Up = false, SeT_JER_Up = false, SeM_JER_Up = false, SeL_JER_Up = false,
		ThT_JER_Up = false, ThM_JER_Up = false, ThL_JER_Up = false, FoT_JER_Up = false, FoM_JER_Up = false, FoL_JER_Up = false;
		
		Int_t nJet30_JER_Up = 0;
		
		Double_t Jet_pt0_JER_Up = Jet_corr_JERUp[0] * Jet_rawPt[0] * Jet_corr_JER[0];
		Double_t Jet_pt1_JER_Up = Jet_corr_JERUp[1] * Jet_rawPt[1]* Jet_corr_JER[1];
		Double_t Jet_pt2_JER_Up = Jet_corr_JERUp[2] * Jet_rawPt[2]* Jet_corr_JER[2];
		Double_t Jet_pt3_JER_Up = Jet_corr_JERUp[3] * Jet_rawPt[3]* Jet_corr_JER[3];


		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_up = 0;
			Jet_pt_up = Jet_corr_JERUp[j]* Jet_rawPt[j] *Jet_corr_JER[j];
			if(Jet_pt_up<30 || fabs(Jet_eta[j])>2.4) continue;
         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9) nJet30_JER_Up++;
               		}
               		else nJet30_JER_Up++;
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30_JER_Up++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30_JER_Up++;
		}
		if(nJet30_JER_Up>3){
        	if(Jet_pt0_JER_Up>30 && fabs(Jet_eta[0])<2.4){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_JER_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_JER_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_JER_Up=true;
         	}
        	if(Jet_pt1_JER_Up>30 && fabs(Jet_eta[1])<2.4){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_JER_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_JER_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_JER_Up=true;
         	}
			if(Jet_pt2_JER_Up>30 && fabs(Jet_eta[2])<2.4){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_JER_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_JER_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_JER_Up=true;
			}
			if(Jet_pt3_JER_Up>30 && fabs(Jet_eta[3])<2.4){
            	if(Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.8958) FoT_JER_Up=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.6324) FoM_JER_Up=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.2219) FoL_JER_Up=true;
      		}
      	}	
		if(FiT_JER_Up && SeT_JER_Up){
			passT2_JER_Up=true;
			if(ThT_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;
				else if(FoM_JER_Up) passT4_JER_Up=true;					
				else if(FoL_JER_Up) passT4_JER_Up=true;
				else if(!FoL_JER_Up && !FoM_JER_Up && !FoT_JER_Up) passT4_JER_Up=true;											
			} 
			else if(ThM_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;
				else if(FoM_JER_Up) passT4_JER_Up=true;
			} 
			else if(ThL_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;				
			}
		}					
		else if(FiT_JER_Up && SeM_JER_Up){
			passT2_JER_Up=true;
			if(ThT_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;
				else if(FoM_JER_Up) passT4_JER_Up=true;					
			}
			else if(ThM_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;
			}				
		}					
		else if(FiM_JER_Up && SeT_JER_Up){
			passT2_JER_Up=true;
			if(ThT_JER_Up){
				passT3_JER_Up=true;
				if(FoT_JER_Up) passT4_JER_Up=true;
				else if(FoM_JER_Up) passT4_JER_Up=true;					
			}
			else if(ThM_JER_Up){
				passT3_JER_Up=true;				
				if(FoT_JER_Up) passT4_JER_Up=true;
			}
		} 
		else if(FiM_JER_Up && SeM_JER_Up){
			passT2_JER_Up=true;
			if(ThT_JER_Up){
				passT3_JER_Up=true;	
				if(FoT_JER_Up) passT4_JER_Up=true;
			}
		}
	//------------------------------JER Down

		bool passT4_JER_Down = false, passT3_JER_Down = false, passT2_JER_Down = false, FiT_JER_Down = false,
		FiM_JER_Down = false, FiL_JER_Down = false, SeT_JER_Down = false, SeM_JER_Down = false, SeL_JER_Down = false,
		ThT_JER_Down = false, ThM_JER_Down = false, ThL_JER_Down = false, FoT_JER_Down = false, FoM_JER_Down = false, FoL_JER_Down = false;
		
		Int_t nJet30_JER_Down = 0;
		
		Double_t Jet_pt0_JER_Down = Jet_corr_JERDown[0] * Jet_rawPt[0] * Jet_corr[0];
		Double_t Jet_pt1_JER_Down = Jet_corr_JERDown[1] * Jet_rawPt[1] * Jet_corr[1];
		Double_t Jet_pt2_JER_Down = Jet_corr_JERDown[2] * Jet_rawPt[2] * Jet_corr[2];
		Double_t Jet_pt3_JER_Down = Jet_corr_JERDown[3] * Jet_rawPt[3] * Jet_corr[3];


		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_up = 0;
			Jet_pt_up = Jet_corr_JERDown[j]* Jet_rawPt[j] * Jet_corr[j];
			if(Jet_pt_up<30 || fabs(Jet_eta[j])>2.4) continue;
         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9) nJet30_JER_Down++;
               		}
               		else nJet30_JER_Down++;
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30_JER_Down++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30_JER_Down++;
		}
		if(nJet30_JER_Down>3){
        	if(Jet_pt0_JER_Down>30 && fabs(Jet_eta[0])<2.4){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_JER_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_JER_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_JER_Down=true;
         	}
        	if(Jet_pt1_JER_Down>30 && fabs(Jet_eta[1])<2.4){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_JER_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_JER_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_JER_Down=true;
         	}
			if(Jet_pt2_JER_Down>30 && fabs(Jet_eta[2])<2.4){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_JER_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_JER_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_JER_Down=true;
			}
			if(Jet_pt3_JER_Down>30 && fabs(Jet_eta[3])<2.4){
            	if(Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.8958) FoT_JER_Down=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.6324) FoM_JER_Down=true;
            	else if (Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3]>0.2219) FoL_JER_Down=true;
      		}
      	}	
		if(FiT_JER_Down && SeT_JER_Down){
			passT2_JER_Down=true;
			if(ThT_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;
				else if(FoM_JER_Down) passT4_JER_Down=true;					
				else if(FoL_JER_Down) passT4_JER_Down=true;
				else if(!FoL_JER_Down && !FoM_JER_Down && !FoT_JER_Down) passT4_JER_Down=true;											
			} 
			else if(ThM_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;
				else if(FoM_JER_Down) passT4_JER_Down=true;
			} 
			else if(ThL_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;				
			}
		}					
		else if(FiT_JER_Down && SeM_JER_Down){
			passT2_JER_Down=true;
			if(ThT_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;
				else if(FoM_JER_Down) passT4_JER_Down=true;					
			}
			else if(ThM_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;
			}				
		}					
		else if(FiM_JER_Down && SeT_JER_Down){
			passT2_JER_Down=true;
			if(ThT_JER_Down){
				passT3_JER_Down=true;
				if(FoT_JER_Down) passT4_JER_Down=true;
				else if(FoM_JER_Down) passT4_JER_Down=true;					
			}
			else if(ThM_JER_Down){
				passT3_JER_Down=true;				
				if(FoT_JER_Down) passT4_JER_Down=true;
			}
		} 
		else if(FiM_JER_Down && SeM_JER_Down){
			passT2_JER_Down=true;
			if(ThT_JER_Down){
				passT3_JER_Down=true;	
				if(FoT_JER_Down) passT4_JER_Down=true;
			}
		}
	//=====================Trigger weight calculation========================================
		vector<float> CSV2=CSV;
    	std::sort(CSV2.begin(),CSV2.end(),[](float a, float b){return a>b;});
    	if(CSV2.size()){
      		MinBtag=CSV2[min((int)CSV2.size()-1,2)];
      		PreMinBtag=CSV2[min((int)CSV2.size()-1,1)];
    	}
		float SumCalo=0.;
    	for(unsigned int c=0; c<min(ntrgObjects_caloJets,4); ++c) SumCalo+=trgObjects_caloJets_pt[c];

	//----------------------Weights for Di90Di30------------------------------
		Double_t totalWeight = 0.0, totalWeight_up = 0.0, totalWeight_down = 0.0;
    	vector<double> weights;
    	vector<double> weights_up;
    	vector<double> weights_down;

    	weights.push_back(weight(l1, SumCalo));
    	weights.push_back(weight(t1, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	weights.push_back(weight(t2, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
    	weights.push_back(weight(t3, MinBtag));
    	weights.push_back(weight(t4, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    	weights.push_back(weight(t5, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));

    	weights_up.push_back(weight_error_up(l1, SumCalo));
    	weights_up.push_back(weight_error_up(t1, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	weights_up.push_back(weight_error_up(t2, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
    	weights_up.push_back(weight_error_up(t3, MinBtag));
    	weights_up.push_back(weight_error_up(t4, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    	weights_up.push_back(weight_error_up(t5, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));

    	weights_down.push_back(weight_error_down(l1, SumCalo));
    	weights_down.push_back(weight_error_down(t1, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	weights_down.push_back(weight_error_down(t2, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
    	weights_down.push_back(weight_error_down(t3, MinBtag));
    	weights_down.push_back(weight_error_down(t4, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    	weights_down.push_back(weight_error_down(t5, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));
    
    	for (unsigned int i = 0; i < weights.size(); ++i){
    		if(TMath::IsNaN(weights[i])) weights[i] = 1.0;
       		if(TMath::IsNaN(weights_up[i])) weights_up[i] = 1.0;
       		if(TMath::IsNaN(weights_down[i])) weights_down[i] = 1.0;
    	}
    	
    	std::vector<double> var_up;
    	std::vector<double> var_down;
    	totalWeight = fit1->Eval(SumCalo)*fit2->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)])*fit3->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)])*fitcsv->Eval(MinBtag)*fit4->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)])*fit5->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]);
    	
    	if(totalWeight< 0.0) totalWeight =0.0;
    	var_up.push_back(fit1->Eval(SumCalo)+fit1_unc->Eval(SumCalo));
    	var_up.push_back(fit2->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)])+fit2_unc->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	var_up.push_back(fit3->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)])+fit3_unc->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
    	var_up.push_back(fitcsv->Eval(MinBtag)+fitcsv_unc->Eval(MinBtag));
    	var_up.push_back(fit4->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)])+fit4_unc->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    	var_up.push_back(fit5->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)])+fit5_unc->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));

    	var_down.push_back(fit1->Eval(SumCalo)-fit1_unc->Eval(SumCalo));
    	var_down.push_back(fit2->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)])-fit2_unc->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	var_down.push_back(fit3->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)])-fit3_unc->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
    	var_down.push_back(fitcsv->Eval(MinBtag)-fitcsv_unc->Eval(MinBtag));
    	var_down.push_back(fit4->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)])-fit4_unc->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    	var_down.push_back(fit5->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)])-fit5_unc->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));

    	for(unsigned int i = 0; i < var_down.size(); i++){
    		if(var_up[i] >= 1.0) var_up[i] = 1.0;
    		if(var_down[i] < 0.0 ) var_down[i] = 0.0; 
    	}

    	totalWeight_up = var_up[0]*var_up[1]*var_up[2]*var_up[3]*var_up[4]*var_up[5];
    	totalWeight_down = var_down[0]*var_down[1]*var_down[2]*var_down[3]*var_down[4]*var_down[5];

    	if(totalWeight==0 || totalWeight==-0){
    		totalWeight_up = 0.0;
    		totalWeight_down = 0.0;
    	}

    	//totalWeight = weights[0] *weights[1] *weights[2] *weights[3] *weights[4] *weights[5];
    	//totalWeight_up   = fit1_up->Eval(SumCalo)  *fit2_up->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)])  *fit3_up->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)])  *fitcsv_up->Eval(MinBtag)*fit4_up->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)])  *fit5_up->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]);
    	//totalWeight_down = fit1_down->Eval(SumCalo)*fit2_down->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)])*fit3_down->Eval(trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)])*fitcsv_down->Eval(MinBtag)*fit4_down->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)])*fit5_down->Eval(trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]);

    	if(totalWeight_up < totalWeight){
    		std::cout << "Weight: " << totalWeight << endl; 
    		std::cout << "Weight up: " << totalWeight_up << endl; 
    		std::cout << "Weight down: " << totalWeight_down << endl; 
    		std::cout << SumCalo << endl;
    		std::cout << trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)] << endl;
    		std::cout << trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)] << endl;
    		std::cout << MinBtag << endl;
    		std::cout << trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)] << endl;
    		std::cout << trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)] << endl;

    	}
    	/*totalWeight = weights[0] *weights[1] *weights[2] *weights[3] *weights[4] *weights[5];
    	totalWeight_up = (weights[0]+weights_up[0])*(weights[1]+weights_up[1])*(weights[2]+weights_up[2])*(weights[3]+weights_up[3])*(weights[4]+weights_up[4])*(weights[5]+weights_up[5]);
       	totalWeight_down = (weights[0]-weights_down[0])*(weights[1]-weights_down[1])*(weights[2]-weights_down[2])*(weights[3]-weights_down[3])*(weights[4]-weights_down[4])*(weights[5]-weights_down[5]);

       	WeightSum += totalWeight;
       	WeightSum_up += totalWeight_up;
       	WeightSum_down += totalWeight_down;*/
//--------------------------Weights for Quad45-------------------------------------------------
		Double_t totalWeight45 = 0.0, totalWeight45_up = 0.0, totalWeight45_down = 0.0;
    	vector<double> weights45;
    	vector<double> weights45_up;
    	vector<double> weights45_down;

      	weights45.push_back(weight(l145, SumCalo));
      	weights45.push_back(weight(t145, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
      	weights45.push_back(weight(t245, MinBtag));
      	weights45.push_back(weight(t345, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    
    	weights45_up.push_back(weight_error_up(l145, SumCalo));
    	weights45_up.push_back(weight_error_up(t145, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	weights45_up.push_back(weight_error_up(t245, MinBtag));
    	weights45_up.push_back(weight_error_up(t345, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));

    	weights45_down.push_back(weight_error_down(l145, SumCalo));
    	weights45_down.push_back(weight_error_down(t145, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
    	weights45_down.push_back(weight_error_down(t245, MinBtag));
    	weights45_down.push_back(weight_error_down(t345, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
    
    	for (unsigned int i = 0; i < weights45.size(); ++i){
    		if(TMath::IsNaN(weights45[i])) weights45[i] = 1.0;
       		if(TMath::IsNaN(weights45_up[i])) weights45_up[i] = 1.0;
       		if(TMath::IsNaN(weights45_down[i])) weights45_down[i] = 1.0;
    	}
    	totalWeight45 = weights45[0] *weights45[1] *weights45[2] *weights45[3];
    	totalWeight45_up = (weights45[0]+weights45_up[0])*(weights45[1]+weights45_up[1])*(weights45[2]+weights45_up[2])*(weights45[3]+weights45_up[3]);
       	totalWeight45_down = (weights45[0]-weights45_down[0])*(weights45[1]-weights45_down[1])*(weights45[2]-weights45_down[2])*(weights45[3]-weights45_down[3]);

       	/*WeightSum45 += totalWeight45;
       	WeightSum45_up += totalWeight45_up;
       	WeightSum45_down += totalWeight45_down;*/

	//=================Fill histograms====================================
		TLorentzVector j1;
		TLorentzVector j2;
		TLorentzVector j1_JESUp, j2_JESUp, j1_JESDown, j2_JESDown;
		TLorentzVector j1_JERUp, j2_JERUp, j1_JERDown, j2_JERDown;
		Double_t xsec_weight = 0.0;
		xsec_weight = xsec*lumi*1000./nevts;

		j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
		j1_JESUp.SetPtEtaPhiM(Jet_pt0_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2_JESUp.SetPtEtaPhiM(Jet_pt1_Up, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
		j1_JESDown.SetPtEtaPhiM(Jet_pt0_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2_JESDown.SetPtEtaPhiM(Jet_pt1_Down, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
		j1_JERUp.SetPtEtaPhiM(Jet_pt0_JER_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2_JERUp.SetPtEtaPhiM(Jet_pt1_JER_Up, Jet_eta[1], Jet_phi[1], Jet_mass[1]);			
		j1_JERDown.SetPtEtaPhiM(Jet_pt0_JER_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2_JERDown.SetPtEtaPhiM(Jet_pt1_JER_Down, Jet_eta[1], Jet_phi[1], Jet_mass[1]);		
		TLorentzVector s = j1+j2;
		TLorentzVector s_Up = j1_JESUp+j2_JESUp;
		TLorentzVector s_Down = j1_JESDown+j2_JESDown;
		TLorentzVector s_JERUp = j1_JERUp+j2_JERUp;
		TLorentzVector s_JERDown = j1_JERDown+j2_JERDown;	

		if(passT4){
			histograms["mjjsr"]->Fill(s.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_PUUp"]->Fill(s.M(), xsec_weight*puWeightUp*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_PUDown"]->Fill(s.M(), xsec_weight*puWeightDown*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_jesUp"]->Fill(s_Up.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_jesDown"]->Fill(s_Down.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_jerUp"]->Fill(s_JERUp.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_jerDown"]->Fill(s_JERDown.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_DeepCSVUp"]->Fill(s.M(),  xsec_weight*puWeight*DeepCSVEventWeight_up*totalWeight);
			histograms["mjjsr_DeepCSVDown"]->Fill(s.M(),  xsec_weight*puWeight*DeepCSVEventWeight_down*totalWeight);
			histograms["mjjsr_TriggerUp"]->Fill(s.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight_up);
			histograms["mjjsr_TriggerDown"]->Fill(s.M(), xsec_weight*puWeight*DeepCSVEventWeight*totalWeight_down);
			histograms["mjjsr_PDFUp"]->Fill(s.M(), xsec_weight*pdf_Weight_up*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_PDFDown"]->Fill(s.M(), xsec_weight*pdf_Weight_down*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_alphaUp"]->Fill(s.M(), xsec_weight*alpha_up*puWeight*DeepCSVEventWeight*totalWeight);
			histograms["mjjsr_alphaDown"]->Fill(s.M(), xsec_weight*alpha_down*puWeight*DeepCSVEventWeight*totalWeight);

			histograms["mjjcr2"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45);
			histograms["mjjcr2_DeepCSVUp"]->Fill(s.M(),  puWeight*DeepCSVEventWeight_up*totalWeight45);
			histograms["mjjcr2_DeepCSVDown"]->Fill(s.M(), puWeight*DeepCSVEventWeight_down*totalWeight45);
			histograms["mjjcr2_TriggerUp"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45_up);
			histograms["mjjcr2_TriggerDown"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45_down);
		}
		else{
			if(passT2){
				histograms["mjjcr1"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight);
				histograms["mjjcr1_DeepCSVUp"]->Fill(s.M(),  puWeight*DeepCSVEventWeight_up*totalWeight);
				histograms["mjjcr1_DeepCSVDown"]->Fill(s.M(), puWeight*DeepCSVEventWeight_down*totalWeight);
				histograms["mjjcr1_TriggerUp"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight_up);
				histograms["mjjcr1_TriggerDown"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight_down);

				histograms["mjjcr3"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45);
				histograms["mjjcr3_DeepCSVUp"]->Fill(s.M(),  puWeight*DeepCSVEventWeight_up*totalWeight45);
				histograms["mjjcr3_DeepCSVDown"]->Fill(s.M(),puWeight*DeepCSVEventWeight_down*totalWeight45);
				histograms["mjjcr3_TriggerUp"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45_up);
				histograms["mjjcr3_TriggerDown"]->Fill(s.M(), puWeight*DeepCSVEventWeight*totalWeight45_down);
			}
		}	

	}//loop over entries

	TFile* outfile = new TFile("/fdata/hepx/store/user/delgado_andrea/mjjShapes/mjj_syst_"+name+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;
}
