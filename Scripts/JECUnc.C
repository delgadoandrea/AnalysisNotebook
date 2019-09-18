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

using std::cout;
using std::endl;
using std::vector;

class dRComb{
public:

	dRComb(int a, int b, double c): i(a), j(b), dR(c){}
	int i;
	int j;
	double dR;

};

bool sortBySmallestDR (dRComb i, dRComb j){
	return (i.dR < j.dR);
}

bool sortByLargestM (dRComb i, dRComb j){
	return (i.dR > j.dR);
}

map <string, TH1*> histograms;
map <string, TH1*> histograms2D;

void JECUnc (TChain *chain, TString name){
   	if (chain == 0 ){
		std::cout << "No tree found!" <<std::endl;
		return;
   	} 
	Float_t GenJet_phi[50];
	Float_t GenJet_eta[50];
	Float_t GenJet_pt[50];
	Float_t GenJet_mass[50];	
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
  	Float_t Jet_corr_JECUp[50];
  	Float_t Jet_corr_JECDown[50];
  	Int_t Jet_numberOfDaughters[50];
  	Float_t Jet_muEF[50];
  	Int_t Jet_chMult[50];
  	Int_t Jet_nhMult[50];		
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v;
	Int_t HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v;

	Int_t nJet;
	Int_t nGenJet;	

	chain->SetBranchAddress("GenJet_eta",GenJet_eta);			
	chain->SetBranchAddress("GenJet_phi",GenJet_phi);			
	chain->SetBranchAddress("GenJet_pt",GenJet_pt);			
	chain->SetBranchAddress("GenJet_mass",GenJet_mass);		
	chain->SetBranchAddress("Jet_pt", Jet_pt);	
	chain->SetBranchAddress("Jet_eta", Jet_eta);	
	chain->SetBranchAddress("Jet_phi", Jet_phi);	
	chain->SetBranchAddress("Jet_mass", Jet_mass);	
	chain->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);			
	chain->SetBranchAddress("Jet_btagDeepCSVb", Jet_btagDeepCSVb);				
	chain->SetBranchAddress("Jet_btagDeepCSVbb", Jet_btagDeepCSVbb);					
	chain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA);
	chain->SetBranchAddress("Jet_chHEF", Jet_chHEF);
  	chain->SetBranchAddress("Jet_neHEF", Jet_neHEF);
  	chain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF);
  	chain->SetBranchAddress("Jet_neEmEF", Jet_neHEF);
  	chain->SetBranchAddress("Jet_numberOfDaughters", &Jet_numberOfDaughters);
  	chain->SetBranchAddress("Jet_muEF", Jet_muEF);
  	chain->SetBranchAddress("Jet_chMult",&Jet_chMult);
  	chain->SetBranchAddress("Jet_nhMult", &Jet_nhMult);
	chain->SetBranchAddress("nJet",&nJet);			
	chain->SetBranchAddress("nGenJet",&nGenJet);			
	chain->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v);
	chain->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v);
	chain->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v);
	chain->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v", &HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v);	
	chain->SetBranchAddress("Jet_corr_JECUp",Jet_corr_JECUp);
	chain->SetBranchAddress("Jet_corr_JECDown",Jet_corr_JECDown);

	histograms["mjjsr"] = new TH1F("mjjsr", "mjj from two leading jets in event for signal region ", 120,0, 1200);	
	histograms["mjjcr1"] = new TH1F("mjjcr1", "mjj from two leading jets in event for control region 1", 120,0, 1200);	
	histograms["mjjcr2"] = new TH1F("mjjcr2", "mjj from two leading jets in event for control region 2", 120,0, 1200);		
	histograms["mjjcr3"] = new TH1F("mjjcr3", "mjj from two leading jets in event for control region 3", 120,0, 1200);			
	histograms["mjjcr0p"] = new TH1F("mjjcr0p", "mjj frommjjcr3 two leading jets in event for control region 0'", 120,0, 1200);	
	histograms["mjjcr1p"] = new TH1F("mjjcr1p", "mjj from two leading jets in event for control region 1'", 120,0, 1200);	
	histograms["mjjcr2p"] = new TH1F("mjjcr2p", "mjj from two leading jets in event for control region 2'", 120,0, 1200);		
	histograms["mjjcr3p"] = new TH1F("mjjcr3p", "mjj from two leading jets in event for control region 3'", 120,0, 1200);

	histograms["mjjsr_jup"] = new TH1F("mjjsr_jup", "mjj from two leading jets in event for signal region +1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr1_jup"] = new TH1F("mjjcr1_jup", "mjj from two leading jets in event for control region 1 +1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr2_jup"] = new TH1F("mjjcr2_jup", "mjj from two leading jets in event for control region 2 +1#sigma JEC", 120,0, 1200);		
	histograms["mjjcr3_jup"] = new TH1F("mjjcr3_jup", "mjj from two leading jets in event for control region 3 +1#sigma JEC", 120,0, 1200);			
	histograms["mjjcr0p_jup"] = new TH1F("mjjcr0p_jup", "mjj frommjjcr3 two leading jets in event for control region 0' +1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr1p_jup"] = new TH1F("mjjcr1p_jup", "mjj from two leading jets in event for control region 1' +1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr2p_jup"] = new TH1F("mjjcr2p_jup", "mjj from two leading jets in event for control region 2' +1#sigma JEC", 120,0, 1200);		
	histograms["mjjcr3p_jup"] = new TH1F("mjjcr3p_jup", "mjj from two leading jets in event for control region 3' +1#sigma JEC", 120,0, 1200);

	histograms["mjjsr_jdown"] = new TH1F("mjjsr_jdown", "mjj from two leading jets in event for signal region -1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr1_jdown"] = new TH1F("mjjcr1_jdown", "mjj from two leading jets in event for control region 1 -1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr2_jdown"] = new TH1F("mjjcr2_jdown", "mjj from two leading jets in event for control region 2 -1#sigma JEC", 120,0, 1200);		
	histograms["mjjcr3_jdown"] = new TH1F("mjjcr3_jdown", "mjj from two leading jets in event for control region 3 -1#sigma JEC", 120,0, 1200);			
	histograms["mjjcr0p_jdown"] = new TH1F("mjjcr0p_jdown", "mjj frommjjcr3 two leading jets in event for control region 0' -1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr1p_jdown"] = new TH1F("mjjcr1p_jdown", "mjj from two leading jets in event for control region 1' -1#sigma JEC", 120,0, 1200);	
	histograms["mjjcr2p_jdown"] = new TH1F("mjjcr2p_jdown", "mjj from two leading jets in event for control region 2' -1#sigma JEC", 120,0, 1200);		
	histograms["mjjcr3p_jdown"] = new TH1F("mjjcr3p_jdown", "mjj from two leading jets in event for control region 3' -1#sigma JEC", 120,0, 1200);

	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << chain->GetEntries() << endl;
	Int_t nPass = 0;	
	
	for (Int_t indx=0;indx<chain->GetEntries();indx++){ 			
	//for (Int_t indx=0;indx<10;indx++){ 
		chain->GetEntry(indx);

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
		
		for (int j=0; j<nJet; j++){
			if(Jet_pt[j]<30 || fabs(Jet_eta[j])>2.6) continue;
         	if(fabs(Jet_eta[j])<2.7){
	        	if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
	            	if(fabs(Jet_eta[j])<2.4){
	               		if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9) nJet30++;
               		}
               		else nJet30++;
            	}
         	}
        	else if(fabs(Jet_eta[j])<3.0){
            	if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) nJet30++;
         	}
         	else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) nJet30++;
		}
		
		if(nJet30>3){
        	if(Jet_pt[0]>30 && fabs(Jet_eta[0])<2.6){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL=true;
         	}
        	if(Jet_pt[1]>30 && fabs(Jet_eta[1])<2.6){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL=true;
         	}
			if(Jet_pt[2]>30 && fabs(Jet_eta[2])<2.6){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL=true;
			}
			if(Jet_pt[3]>30 && fabs(Jet_eta[3])<2.6){
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

//===================== JEC Up =======================================

		bool passT4_Up = false;
		bool passT3_Up = false;
		bool passT2_Up = false;						

		bool FiT_Up = false;
		bool FiM_Up = false;
		bool FiL_Up = false;

		bool SeT_Up = false;
		bool SeM_Up = false;
		bool SeL_Up = false;

		bool ThT_Up = false;
		bool ThM_Up = false;
		bool ThL_Up = false;

		bool FoT_Up = false;
		bool FoM_Up = false;
		bool FoL_Up = false;
		Int_t nJet30_Up = 0;
		Float_t Jet_pt0_Up = Jet_corr_JECUp[0] * Jet_pt[0];
		Double_t Jet_pt1_Up = Jet_corr_JECUp[1] * Jet_pt[1];
		Double_t Jet_pt2_Up = Jet_corr_JECUp[2] * Jet_pt[2];
		Double_t Jet_pt3_Up = Jet_corr_JECUp[3] * Jet_pt[3];


		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_up = 0;
			Jet_pt_up = Jet_corr_JECUp[j]* Jet_pt[j];
			if(Jet_pt_up<30 || fabs(Jet_eta[j])>2.6) continue;
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
        	if(Jet_pt0_Up>30 && fabs(Jet_eta[0])<2.6){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_Up=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_Up=true;
         	}
        	if(Jet_pt1_Up>30 && fabs(Jet_eta[1])<2.6){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_Up=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_Up=true;
         	}
			if(Jet_pt2_Up>30 && fabs(Jet_eta[2])<2.6){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_Up=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_Up=true;
			}
			if(Jet_pt3_Up>30 && fabs(Jet_eta[3])<2.6){
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

//===================== JEC Down =======================================

		bool passT4_Down = false;
		bool passT3_Down = false;
		bool passT2_Down = false;						

		bool FiT_Down = false;
		bool FiM_Down = false;
		bool FiL_Down = false;

		bool SeT_Down = false;
		bool SeM_Down = false;
		bool SeL_Down = false;

		bool ThT_Down = false;
		bool ThM_Down = false;
		bool ThL_Down = false;

		bool FoT_Down = false;
		bool FoM_Down = false;
		bool FoL_Down = false;
		Int_t nJet30_Down = 0;
		Float_t Jet_pt0_Down = Jet_corr_JECDown[0]*Jet_pt[0];
		Double_t Jet_pt1_Down = Jet_corr_JECDown[1] * Jet_pt[1];
		Double_t Jet_pt2_Down = Jet_corr_JECDown[2] * Jet_pt[2];
		Double_t Jet_pt3_Down = Jet_corr_JECDown[3] * Jet_pt[3];

		for (int j=0; j<nJet; j++){
			Double_t Jet_pt_Down = 0;
			Jet_pt_Down = Jet_corr_JECDown[j]* Jet_pt[j];
			if(Jet_pt_Down<30 || fabs(Jet_eta[j])>2.6) continue;
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
        	if(Jet_pt0_Down>30 && fabs(Jet_eta[0])<2.6){
            	if(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.8958) FiT_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.6324) FiM_Down=true;
            	else if (Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]>0.2219) FiL_Down=true;
         	}
        	if(Jet_pt1_Down>30 && fabs(Jet_eta[1])<2.6){
            	if(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.8958) SeT_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.6324) SeM_Down=true;
            	else if (Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]>0.2219) SeL_Down=true;
         	}
			if(Jet_pt2_Down>30 && fabs(Jet_eta[2])<2.6){
            if(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.8958) ThT_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.6324) ThM_Down=true;
            else if (Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]>0.2219) ThL_Down=true;
			}
			if(Jet_pt3_Down>30 && fabs(Jet_eta[3])<2.6){
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

//============================= mjj Calculation ===============================
		
		if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1){						
			TLorentzVector j1;
			TLorentzVector j2;
			TLorentzVector j1_Up;
			TLorentzVector j2_Up;
			TLorentzVector j1_Down;
			TLorentzVector j2_Down;

			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_Up.SetPtEtaPhiM(Jet_pt0_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_Up.SetPtEtaPhiM(Jet_pt1_Up, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_Down.SetPtEtaPhiM(Jet_pt0_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_Down.SetPtEtaPhiM(Jet_pt1_Down, Jet_eta[1], Jet_phi[1], Jet_mass[1]);

			TLorentzVector s = j1+j2;
			TLorentzVector s_Up = j1_Up+j2_Up;
			TLorentzVector s_Down = j1_Down+j2_Down;

			/*if (passT2 && !passT2_Down){
				cout << "JEC's affecting selection: " << indx << endl;
				if(nJet30>3 && nJet30_Down<4) cout << "	Not enough jets passing selction" <<endl;
				else{
					cout << "JEt pt 0: " << Jet_pt0_Down << endl;
					cout << "JEt pt 1: " << Jet_pt1_Down << endl;
				}
			}*/
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
			//=================  Plots Up ====================
			if(passT4_Up){
				histograms["mjjsr_jup"]->Fill(s_Up.M());
			} 
			else{
				if(passT2_Up){
					histograms["mjjcr1_jup"]->Fill(s_Up.M());						
				}
				if(passT3_Up){
					histograms["mjjcr0p_jup"]->Fill(s_Up.M());
				}
				else if(passT2_Up){
					histograms["mjjcr1p_jup"]->Fill(s_Up.M());
				}
			} 
			//=================  Plots Down ====================
			if(passT4_Down){
				histograms["mjjsr_jdown"]->Fill(s_Down.M());
			} 
			else{
				if(passT2_Down){
					histograms["mjjcr1_jdown"]->Fill(s_Down.M());						
				}
				if(passT3_Down){
					histograms["mjjcr0p_jdown"]->Fill(s_Down.M());
				}
				else if(passT2_Down){
					histograms["mjjcr1p_jdown"]->Fill(s_Down.M());
				}
			} 
		} 
		else if(HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v==1){
			TLorentzVector j1;
			TLorentzVector j2;
			TLorentzVector j1_Up;
			TLorentzVector j2_Up;
			TLorentzVector j1_Down;
			TLorentzVector j2_Down;

			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_Up.SetPtEtaPhiM(Jet_pt0_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_Up.SetPtEtaPhiM(Jet_pt1_Up, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_Down.SetPtEtaPhiM(Jet_pt0_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_Down.SetPtEtaPhiM(Jet_pt1_Down, Jet_eta[1], Jet_phi[1], Jet_mass[1]);

			TLorentzVector s = j1+j2;
			TLorentzVector s_Up = j1_Up+j2_Up;
			TLorentzVector s_Down = j1_Down+j2_Down;
			
			if(passT4){
				histograms["mjjcr2"]->Fill(s.M());
				//histograms["test"]->Fill(indx,2.);					
			} 
			else{
				if(passT2){
					histograms["mjjcr3"]->Fill(s.M());
					//histograms["test"]->Fill(indx,3.);						
				}
				if(passT3){
					histograms["mjjcr2p"]->Fill(s.M());
				}
				else if(passT2){
					histograms["mjjcr3p"]->Fill(s.M());
				}
			}
			//=================  Plots Up ====================
			if(passT4_Up){
				histograms["mjjcr2_jup"]->Fill(s_Up.M());					
			} 
			else{
				if(passT2_Up){
					histograms["mjjcr3_jup"]->Fill(s_Up.M());						
				}
				if(passT3_Up){
					histograms["mjjcr2p_jup"]->Fill(s_Up.M());
				}
				else if(passT2_Up){
					histograms["mjjcr3p_jup"]->Fill(s_Up.M());
				}
			}

			//=================  Plots Down ====================
			if(passT4_Down){
				histograms["mjjcr2_jdown"]->Fill(s_Down.M());					
			} 
			else{
				if(passT2_Down){
					histograms["mjjcr3_jdown"]->Fill(s_Down.M());						
				}
				if(passT3_Down){
					histograms["mjjcr2p_jdown"]->Fill(s_Down.M());
				}
				else if(passT2_Down){
					histograms["mjjcr3p_jdown"]->Fill(s_Down.M());
				}
			}

		}// HLT Quad45

	}//loop over entries
	cout << "Events passing selection: " << nPass << endl;
	cout << "Aceptance: " << nPass*100./chain->GetEntries();

	//Save invariant mass histogram to a file
	TFile* outfile = new TFile("/fdata/hepx/store/user/delgado_andrea/freshStartAug2019/BTagCSV/Plots/mjj_DenisSel_"+name+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;
}
