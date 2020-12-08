#include <vector>
#include <iostream>
#include <map>
#include <iterator>

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

using std::cout;
using std::endl;
using std::vector;

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
//--------------------------------------------------------------------
map <string, TH1*> histograms;

void TrigUnc (TString fname, TString name){
	TFile *f = TFile::Open(fname, "READ");
	TTree *t = (TTree*) f->Get("tree");

	TFile *w = TFile::Open("/fdata/hepx/store/user/delgado_andrea/94X/updatedNtuplesJune19/weightHistosAll.root");
	TFile *w45 = TFile::Open("/fdata/hepx/store/user/delgado_andrea/94X/updatedNtuplesJune19/weightHistosAllQuad45.root");

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

   	if (t == 0 ){
		std::cout << "No tree found!" <<std::endl;
		return;
   	} 
	Float_t GenJet_phi[50];
	Float_t GenJet_eta[50];
	Float_t GenJet_pt[50];
	Float_t GenJet_mass[50];	
	Float_t Jet_rawPt[50];
	Float_t Jet_pt[50];		
	Float_t Jet_eta[50];		
	Float_t Jet_phi[50];		
	Float_t Jet_mass[50];	
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
  	Float_t Jet_rawPtAfterSmearing[50];	
  	Int_t Jet_numberOfDaughters[50];
  	Float_t Jet_muEF[50];
  	Int_t Jet_chMult[50];
  	Int_t Jet_nhMult[50];		
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v;
	Float_t trgObjects_caloJets_pt[50];
	Float_t trgObjects_pfJets_pt[50];

	Int_t nJet;
	Int_t nGenJet;	
	Int_t ntrgObjects_caloJets;
	Int_t ntrgObjects_pfJets;

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
	t->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);			
	t->SetBranchAddress("Jet_btagDeepCSVb", Jet_btagDeepCSVb);				
	t->SetBranchAddress("Jet_btagDeepCSVbb", Jet_btagDeepCSVbb);					
	t->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA);
	t->SetBranchAddress("Jet_chHEF", Jet_chHEF);
  	t->SetBranchAddress("Jet_neHEF", Jet_neHEF);
  	t->SetBranchAddress("Jet_chEmEF", Jet_chEmEF);
  	t->SetBranchAddress("Jet_neEmEF", Jet_neHEF);
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
	histograms["mjjcr0p"] = new TH1F("mjjcr0p", "mjj frommjjcr3 two leading jets in event for control region 0'", 120,0, 1200);	
	histograms["mjjcr1p"] = new TH1F("mjjcr1p", "mjj from two leading jets in event for control region 1'", 120,0, 1200);	
	histograms["mjjcr2p"] = new TH1F("mjjcr2p", "mjj from two leading jets in event for control region 2'", 120,0, 1200);		
	histograms["mjjcr3p"] = new TH1F("mjjcr3p", "mjj from two leading jets in event for control region 3'", 120,0, 1200);

	histograms["mjjsrtw"] = new TH1F("mjjsrtw", "mjj from two leading jets in event for signal region using trigger weights", 120,0, 1200);	
	histograms["mjjcr1tw"] = new TH1F("mjjcr1tw", "mjj from two leading jets in event for control region 1 using trigger weights", 120,0, 1200);	
	histograms["mjjcr2tw"] = new TH1F("mjjcr2tw", "mjj from two leading jets in event for control region 2 using trigger weights", 120,0, 1200);		
	histograms["mjjcr3tw"] = new TH1F("mjjcr3tw", "mjj from two leading jets in event for control region 3 using trigger weights", 120,0, 1200);			
	histograms["mjjcr0ptw"] = new TH1F("mjjcr0ptw", "mjj frommjjcr3 two leading jets in event for control region 0' using trigger weights", 120,0, 1200);	
	histograms["mjjcr1ptw"] = new TH1F("mjjcr1ptw", "mjj from two leading jets in event for control region 1' using trigger weightsS", 120,0, 1200);	
	histograms["mjjcr2ptw"] = new TH1F("mjjcr2ptw", "mjj from two leading jets in event for control region 2' using trigger weights", 120,0, 1200);		
	histograms["mjjcr3ptw"] = new TH1F("mjjcr3ptw", "mjj from two leading jets in event for control region 3' using trigger weightsS", 120,0, 1200);

	histograms["mjjsrtw_up"] = new TH1F("mjjsrtw_up", "mjj from two leading jets in event for signal region using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr1tw_up"] = new TH1F("mjjcr1tw_up", "mjj from two leading jets in event for control region 1 using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr2tw_up"] = new TH1F("mjjcr2tw_up", "mjj from two leading jets in event for control region 2 using trigger weights + shift up", 120,0, 1200);		
	histograms["mjjcr3tw_up"] = new TH1F("mjjcr3tw_up", "mjj from two leading jets in event for control region 3 using trigger weights + shift up", 120,0, 1200);			
	histograms["mjjcr0ptw_up"] = new TH1F("mjjcr0ptw_up", "mjj frommjjcr3 two leading jets in event for control region 0' using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr1ptw_up"] = new TH1F("mjjcr1ptw_up", "mjj from two leading jets in event for control region 1' using trigger weights + shift up", 120,0, 1200);	
	histograms["mjjcr2ptw_up"] = new TH1F("mjjcr2ptw_up", "mjj from two leading jets in event for control region 2' using trigger weights + shift up", 120,0, 1200);		
	histograms["mjjcr3ptw_up"] = new TH1F("mjjcr3ptw_up", "mjj from two leading jets in event for control region 3' using trigger weights + shift up", 120,0, 1200);

	histograms["mjjsrtw_down"] = new TH1F("mjjsrtw_down", "mjj from two leading jets in event for signal region using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr1tw_down"] = new TH1F("mjjcr1tw_down", "mjj from two leading jets in event for control region 1 using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr2tw_down"] = new TH1F("mjjcr2tw_down", "mjj from two leading jets in event for control region 2 using trigger weights + shift down", 120,0, 1200);		
	histograms["mjjcr3tw_down"] = new TH1F("mjjcr3tw_down", "mjj from two leading jets in event for control region 3 using trigger weights + shift down", 120,0, 1200);			
	histograms["mjjcr0ptw_down"] = new TH1F("mjjcr0ptw_down", "mjj frommjjcr3 two leading jets in event for control region 0' using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr1ptw_down"] = new TH1F("mjjcr1ptw_down", "mjj from two leading jets in event for control region 1' using trigger weights + shift down", 120,0, 1200);	
	histograms["mjjcr2ptw_down"] = new TH1F("mjjcr2ptw_down", "mjj from two leading jets in event for control region 2' using trigger weights + shift down", 120,0, 1200);		
	histograms["mjjcr3ptw_down"] = new TH1F("mjjcr3ptw_down", "mjj from two leading jets in event for control region 3' using trigger weights + shift down", 120,0, 1200);

	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << t->GetEntries() << endl;
	Int_t nPass = 0;
	Double_t nPass_Trig = 0.0;
	Double_t nPass_Trig_up = 0.0;
	Double_t nPass_Trig_down = 0.0;	
	Double_t WeightSum = 0.0;
	Double_t WeightSum45 = 0.0;
	Double_t WeightSum_up = 0.0;
	Double_t WeightSum45_up = 0.0;
	Double_t WeightSum_down = 0.0;
	Double_t WeightSum45_down = 0.0;
	
	for (Int_t indx=0;indx<t->GetEntries();indx++){ 			
	//for (Int_t indx=0;indx<10;indx++){ 
		t->GetEntry(indx);

		//============================================================

		bool passT4 = false;
		bool passT3 = false;
		bool passT2 = false;						

		bool FiT = false;
		bool FiM = false;
		bool FiL = false;

		bool SeT = false;
		bool SeM = false;
		bool SeL = false;

		bool ThT = false;
		bool ThM = false;
		bool ThL = false;

		bool FoT = false;
		bool FoM = false;
		bool FoL = false;
		Int_t nJet30 = 0;

    	float MaxBtag=0.;
    	float MinBtag=1.;
    	float PreMinBtag=1.;
    	vector<TLorentzVector> PFjets;
		vector<float> CSV;

		for (int j=0; j<nJet; j++){						

			if(Jet_pt[j]<30 || fabs(Jet_eta[j])>2.4) continue;
      		TLorentzVector jet;
      		const float DeepCSV=Jet_btagDeepCSVb[j]+Jet_btagDeepCSVbb[j];
      		jet.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],Jet_mass[j]);
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
		}
		
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
//==================Input variables for trigger weights=================
    	vector<float> CSV2=CSV;
    	std::sort(CSV2.begin(),CSV2.end(),[](float a, float b){return a>b;});
    	if(CSV2.size()){
      		MinBtag=CSV2[min((int)CSV2.size()-1,2)];
      		PreMinBtag=CSV2[min((int)CSV2.size()-1,1)];
    	}

//=============================Trigger Weights=================
		float SumCalo=0.;
    	for(unsigned c=0; c<min(ntrgObjects_caloJets,4); ++c) SumCalo+=trgObjects_caloJets_pt[c];

//=============================Weights for Di90Di30======================================================
		Double_t totalWeight = 0.0;
		Double_t totalWeight_up = 0.0;
		Double_t totalWeight_down = 0.0;
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
    
    	for (int i = 0; i < weights.size(); ++i){
    		if(TMath::IsNaN(weights[i])) weights[i] = 1.0;
       		if(TMath::IsNaN(weights_up[i])) weights_up[i] = 1.0;
       		if(TMath::IsNaN(weights_down[i])) weights_down[i] = 1.0;
    	}
    	totalWeight = weights[0] *weights[1] *weights[2] *weights[3] *weights[4] *weights[5];
    	totalWeight_up = (weights[0]+weights_up[0])*(weights[1]+weights_up[1])*(weights[2]+weights_up[2])*(weights[3]+weights_up[3])*(weights[4]+weights_up[4])*(weights[5]+weights_up[5]);
       	totalWeight_down = (weights[0]-weights_down[0])*(weights[1]-weights_down[1])*(weights[2]-weights_down[2])*(weights[3]-weights_down[3])*(weights[4]-weights_down[4])*(weights[5]-weights_down[5]);

       	/*WeightSum += totalWeight;
       	WeightSum_up += totalWeight_up;
       	WeightSum_down += totalWeight_down;*/

       	if(totalWeight > 0.0) WeightSum += totalWeight;
       	else WeightSum += 1.0;

       	if(totalWeight_up > 0.0) WeightSum_up += totalWeight_up;
 		else WeightSum_up += 1.0;

       	if(totalWeight_down > 0.0) WeightSum_down += totalWeight_down;
       	else WeightSum_down += 1.0;
//=============================Weights for Quad45=========================================================
		Double_t totalWeight45 = 0.0;
		Double_t totalWeight45_up = 0.0;
		Double_t totalWeight45_down = 0.0;
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
    
    	for (int i = 0; i < weights45.size(); ++i){
    		if(TMath::IsNaN(weights45[i])) weights45[i] = 1.0;
       		if(TMath::IsNaN(weights45_up[i])) weights45_up[i] = 1.0;
       		if(TMath::IsNaN(weights45_down[i])) weights45_down[i] = 1.0;
    	}
    	totalWeight45 = weights45[0] *weights45[1] *weights45[2] *weights45[3];
    	totalWeight45_up = (weights45[0]+weights45_up[0])*(weights45[1]+weights45_up[1])*(weights45[2]+weights45_up[2])*(weights45[3]+weights45_up[3]);
       	totalWeight45_down = (weights45[0]-weights45_down[0])*(weights45[1]-weights45_down[1])*(weights45[2]-weights45_down[2])*(weights45[3]-weights45_down[3]);

       	WeightSum45 += totalWeight45;
       	WeightSum45_up += totalWeight45_up;
       	WeightSum45_down += totalWeight45_down;
//============================= mjj Calculation ===============================
		TLorentzVector j1;
		TLorentzVector j2;														

		j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
		j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);

		TLorentzVector s = j1+j2;
		if(passT4){
			nPass_Trig += totalWeight;
			nPass_Trig_up += totalWeight_up;
			nPass_Trig_down += totalWeight_down;
			histograms["mjjsrtw"]->Fill(s.M(), totalWeight);
			histograms["mjjsrtw_up"]->Fill(s.M(), totalWeight_up);
			histograms["mjjsrtw_down"]->Fill(s.M(), totalWeight_down);

			histograms["mjjcr2tw"]->Fill(s.M(), totalWeight45);
			histograms["mjjcr2tw_up"]->Fill(s.M(), totalWeight45_up);
			histograms["mjjcr2tw_down"]->Fill(s.M(), totalWeight45_down);
		}
		else{
			if(passT2){
				histograms["mjjcr1tw"]->Fill(s.M(), totalWeight);
				histograms["mjjcr1tw_up"]->Fill(s.M(), totalWeight_up);
				histograms["mjjcr1tw_down"]->Fill(s.M(), totalWeight_down);

				histograms["mjjcr3tw"]->Fill(s.M(), totalWeight45);
				histograms["mjjcr3tw_up"]->Fill(s.M(), totalWeight45_up);
				histograms["mjjcr3tw_down"]->Fill(s.M(), totalWeight45_down);
			}
			if(passT3){
				histograms["mjjcr0ptw"]->Fill(s.M(),totalWeight);
				histograms["mjjcr0ptw_up"]->Fill(s.M(),totalWeight_up);
				histograms["mjjcr0ptw_down"]->Fill(s.M(),totalWeight_down);

				histograms["mjjcr2ptw"]->Fill(s.M(), totalWeight45);
				histograms["mjjcr2ptw_up"]->Fill(s.M(), totalWeight45_up);
				histograms["mjjcr2ptw_down"]->Fill(s.M(), totalWeight45_down);
			}
			else if(passT2){
				histograms["mjjcr1ptw"]->Fill(s.M(),totalWeight);
				histograms["mjjcr1ptw_up"]->Fill(s.M(),totalWeight_up);
				histograms["mjjcr1ptw_down"]->Fill(s.M(),totalWeight_down);

				histograms["mjjcr3ptw"]->Fill(s.M(), totalWeight45);
				histograms["mjjcr3ptw_up"]->Fill(s.M(), totalWeight45_up);
				histograms["mjjcr3ptw_down"]->Fill(s.M(), totalWeight_down);
			}
		}

		if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1){						
			//=================  Plots ====================
			if(passT4){
				histograms["mjjsr"]->Fill(s.M());
				nPass++;
			} 
			else{
				if(passT2){
					histograms["mjjcr1"]->Fill(s.M());
					//histograms["test"]->Fill(indx,1.);						
				}
				if(passT3){
					histograms["mjjcr0p"]->Fill(s.M());
				}
				else if(passT2){
					histograms["mjjcr1p"]->Fill(s.M());
				}
			} 							
		} 
		else if(HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v==1){
			if(passT4){
				histograms["mjjcr2"]->Fill(s.M());
			} 
			else{
				if(passT2){
					histograms["mjjcr3"]->Fill(s.M());
				}
				if(passT3){
					histograms["mjjcr2p"]->Fill(s.M());
				}
				else if(passT2){
					histograms["mjjcr3p"]->Fill(s.M());
				}
			}
			
		}// HLT Quad45

	}//loop over entries

	Int_t nevts = 0;
	nevts = t->GetEntries();

	cout << "Events passing selection: " << nPass << endl;
	cout << "Acceptance: " << nPass*100./nevts << " +- " << 100.*TMath::Sqrt(nPass*(1.0 - (nPass/nevts)))/nevts << endl;
	/*cout << "Total weighted events: " << WeightSum << endl;
	cout << "Total weighted events up: " << WeightSum_up << endl;
	cout << "Total weighted events down: " << WeightSum_down << endl;
	cout << "Total weighted events 45: " << WeightSum45 << endl;
	cout << "Total weighted events up 45: " << WeightSum45_up << endl;
	cout << "Total weighted events down 45: " << WeightSum45_down << endl;*/
	cout << "After weights: " << nPass_Trig*100./WeightSum << " +- " << 100.*TMath::Sqrt(nPass_Trig*(1.0 - (nPass_Trig/WeightSum)))/WeightSum << endl;
	cout << "Shift up: " << nPass_Trig_up*100./WeightSum_up << " +- " << 100.*TMath::Sqrt(nPass_Trig_up*(1.0 - (nPass_Trig_up/WeightSum_up)))/WeightSum_up << endl;
	cout << "Shift down: " << nPass_Trig_down*100./WeightSum_down << " +- " << 100.*TMath::Sqrt(nPass_Trig_down*(1.0 - (nPass_Trig_down/WeightSum_down)))/WeightSum_down << endl;
	
	cout << "JER shift up(%): " << 100.*((nPass_Trig*100./WeightSum)-(nPass_Trig_up*100./WeightSum_up))/(nPass_Trig*100./WeightSum) << endl;
	cout << "JER shift down(%): " << 100.*((nPass_Trig*100./WeightSum)-(nPass_Trig_up*100./WeightSum_down))/(nPass_Trig*100./WeightSum) << endl;

	/*cout << "After weights: " << nPass_Trig*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_Trig*(1.0 - (nPass_Trig/nevts)))/nevts << endl;
	cout << "Shift up: " << nPass_Trig_up*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_Trig_up*(1.0 - (nPass_Trig_up/nevts)))/nevts << endl;
	cout << "Shift down: " << nPass_Trig_down*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_Trig_down*(1.0 - (nPass_Trig_down/nevts)))/nevts << endl;
	
	cout << "JER shift up(%): " << 100.*((nPass_Trig*100./nevts)-(nPass_Trig_up*100./nevts))/(nPass_Trig*100./nevts) << endl;
	cout << "JER shift down(%): " << 100.*((nPass_Trig*100./nevts)-(nPass_Trig_up*100./nevts))/(nPass_Trig*100./nevts) << endl;*/

	//Save invariant mass histogram to a file
	TFile* outfile = new TFile("/fdata/hepx/store/user/delgado_andrea/mjjShapes/Trigger/mjj_DenisSel_trigger_"+name+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;
}
