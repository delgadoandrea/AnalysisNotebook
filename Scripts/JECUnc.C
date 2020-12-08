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

void JECUnc (TString fname, TString name){
	TFile *f = TFile::Open(fname, "READ");
	TTree *t = (TTree*) f->Get("tree");

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

	Int_t nJet;
	Int_t nGenJet;	

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

	histograms["mjjsr_jesup"] = new TH1F("mjjsr_jesup", "mjj from two leading jets in event for signal region +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1_jesup"] = new TH1F("mjjcr1_jesup", "mjj from two leading jets in event for control region 1 +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2_jesup"] = new TH1F("mjjcr2_jesup", "mjj from two leading jets in event for control region 2 +1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3_jesup"] = new TH1F("mjjcr3_jesup", "mjj from two leading jets in event for control region 3 +1#sigma JES", 120,0, 1200);			
	histograms["mjjcr0p_jesup"] = new TH1F("mjjcr0p_jesup", "mjj frommjjcr3 two leading jets in event for control region 0' +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1p_jesup"] = new TH1F("mjjcr1p_jesup", "mjj from two leading jets in event for control region 1' +1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2p_jesup"] = new TH1F("mjjcr2p_jesup", "mjj from two leading jets in event for control region 2' +1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3p_jesup"] = new TH1F("mjjcr3p_jesup", "mjj from two leading jets in event for control region 3' +1#sigma JES", 120,0, 1200);

	histograms["mjjsr_jesdown"] = new TH1F("mjjsr_jesdown", "mjj from two leading jets in event for signal region -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1_jesdown"] = new TH1F("mjjcr1_jesdown", "mjj from two leading jets in event for control region 1 -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2_jesdown"] = new TH1F("mjjcr2_jesdown", "mjj from two leading jets in event for control region 2 -1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3_jesdown"] = new TH1F("mjjcr3_jesdown", "mjj from two leading jets in event for control region 3 -1#sigma JES", 120,0, 1200);			
	histograms["mjjcr0p_jesdown"] = new TH1F("mjjcr0p_jesdown", "mjj frommjjcr3 two leading jets in event for control region 0' -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr1p_jesdown"] = new TH1F("mjjcr1p_jesdown", "mjj from two leading jets in event for control region 1' -1#sigma JES", 120,0, 1200);	
	histograms["mjjcr2p_jesdown"] = new TH1F("mjjcr2p_jesdown", "mjj from two leading jets in event for control region 2' -1#sigma JES", 120,0, 1200);		
	histograms["mjjcr3p_jesdown"] = new TH1F("mjjcr3p_jesdown", "mjj from two leading jets in event for control region 3' -1#sigma JES", 120,0, 1200);

	histograms["mjjsr_jerdown"] = new TH1F("mjjsr_jerdown", "mjj from two leading jets in event for signal region -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1_jerdown"] = new TH1F("mjjcr1_jerdown", "mjj from two leading jets in event for control region 1 -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2_jerdown"] = new TH1F("mjjcr2_jerdown", "mjj from two leading jets in event for control region 2 -1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3_jerdown"] = new TH1F("mjjcr3_jerdown", "mjj from two leading jets in event for control region 3 -1#sigma JER", 120,0, 1200);			
	histograms["mjjcr0p_jerdown"] = new TH1F("mjjcr0p_jerdown", "mjj frommjjcr3 two leading jets in event for control region 0' -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1p_jerdown"] = new TH1F("mjjcr1p_jerdown", "mjj from two leading jets in event for control region 1' -1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2p_jerdown"] = new TH1F("mjjcr2p_jerdown", "mjj from two leading jets in event for control region 2' -1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3p_jerdown"] = new TH1F("mjjcr3p_jerdown", "mjj from two leading jets in event for control region 3' -1#sigma JER", 120,0, 1200);

	histograms["mjjsr_jerup"] = new TH1F("mjjsr_jerup", "mjj from two leading jets in event for signal region +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1_jerup"] = new TH1F("mjjcr1_jerup", "mjj from two leading jets in event for control region 1 +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2_jerup"] = new TH1F("mjjcr2_jerup", "mjj from two leading jets in event for control region 2 +1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3_jerup"] = new TH1F("mjjcr3_jerup", "mjj from two leading jets in event for control region 3 +1#sigma JER", 120,0, 1200);			
	histograms["mjjcr0p_jerup"] = new TH1F("mjjcr0p_jerup", "mjj frommjjcr3 two leading jets in event for control region 0' +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr1p_jerup"] = new TH1F("mjjcr1p_jerup", "mjj from two leading jets in event for control region 1' +1#sigma JER", 120,0, 1200);	
	histograms["mjjcr2p_jerup"] = new TH1F("mjjcr2p_jerup", "mjj from two leading jets in event for control region 2' +1#sigma JER", 120,0, 1200);		
	histograms["mjjcr3p_jerup"] = new TH1F("mjjcr3p_jerup", "mjj from two leading jets in event for control region 3' +1#sigma JER", 120,0, 1200);


	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << t->GetEntries() << endl;
	Int_t nPass = 0;	
	Int_t nPass_JESUp = 0;
	Int_t nPass_JESDown = 0;
	Int_t nPass_JER = 0;
	Int_t nPass_JERUp = 0;			
	Int_t nPass_JERDown = 0;				
	
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
		
		for (int j=0; j<nJet; j++){						

			if(Jet_pt[j]<30 || fabs(Jet_eta[j])>2.4) continue;
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

//===================== JES Up =======================================

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

//===================== JER Up =======================================

		bool passT4_JER_Up = false;
		bool passT3_JER_Up = false;
		bool passT2_JER_Up = false;						

		bool FiT_JER_Up = false;
		bool FiM_JER_Up = false;
		bool FiL_JER_Up = false;

		bool SeT_JER_Up = false;
		bool SeM_JER_Up = false;
		bool SeL_JER_Up = false;

		bool ThT_JER_Up = false;
		bool ThM_JER_Up = false;
		bool ThL_JER_Up = false;

		bool FoT_JER_Up = false;
		bool FoM_JER_Up = false;
		bool FoL_JER_Up = false;
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
//===================== JER Down =======================================

		bool passT4_JER_Down = false;
		bool passT3_JER_Down = false;
		bool passT2_JER_Down = false;						

		bool FiT_JER_Down = false;
		bool FiM_JER_Down = false;
		bool FiL_JER_Down = false;

		bool SeT_JER_Down = false;
		bool SeM_JER_Down = false;
		bool SeL_JER_Down = false;

		bool ThT_JER_Down = false;
		bool ThM_JER_Down = false;
		bool ThL_JER_Down = false;

		bool FoT_JER_Down = false;
		bool FoM_JER_Down = false;
		bool FoL_JER_Down = false;
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
//============================= mjj Calculation ===============================
		
		if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1){						
			TLorentzVector j1;
			TLorentzVector j2;

			TLorentzVector j1_JESUp;
			TLorentzVector j2_JESUp;			
			TLorentzVector j1_JESDown;
			TLorentzVector j2_JESDown;								
			TLorentzVector j1_JERUp;						
			TLorentzVector j2_JERUp;												
			TLorentzVector j1_JERDown;						
			TLorentzVector j2_JERDown;															

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
				histograms["mjjsr_jesup"]->Fill(s_Up.M());
				nPass_JESUp++;
			} 
			else{
				if(passT2_Up){
					histograms["mjjcr1_jesup"]->Fill(s_Up.M());						
				}
				if(passT3_Up){
					histograms["mjjcr0p_jesup"]->Fill(s_Up.M());
				}
				else if(passT2_Up){
					histograms["mjjcr1p_jesup"]->Fill(s_Up.M());
				}
			} 
			//=================  Plots Down ====================
			if(passT4_Down){
				histograms["mjjsr_jesdown"]->Fill(s_Down.M());
				nPass_JESDown++;
			} 
			else{
				if(passT2_Down){
					histograms["mjjcr1_jesdown"]->Fill(s_Down.M());						
				}
				if(passT3_Down){
					histograms["mjjcr0p_jesdown"]->Fill(s_Down.M());
				}
				else if(passT2_Down){
					histograms["mjjcr1p_jesdown"]->Fill(s_Down.M());
				}
			} 	
			//=================  Plots JER Up ====================
			if(passT4_JER_Up){
				histograms["mjjsr_jerup"]->Fill(s_JERUp.M());
				nPass_JERUp++;
			} 
			else{
				if(passT2_JER_Up){
					histograms["mjjcr1_jerup"]->Fill(s_JERUp.M());						
				}
				if(passT3_JER_Up){
					histograms["mjjcr0p_jerup"]->Fill(s_JERUp.M());
				}
				else if(passT2_JER_Up){
					histograms["mjjcr1p_jerup"]->Fill(s_JERUp.M());
				}
			} 
			//=================  Plots JER Down ====================
			if(passT4_JER_Down){
				histograms["mjjsr_jerdown"]->Fill(s_JERDown.M());
				nPass_JERDown++;
			} 
			else{
				if(passT2_JER_Down){
					histograms["mjjcr1_jerdown"]->Fill(s_JERDown.M());						
				}
				if(passT3_JER_Down){
					histograms["mjjcr0p_jerdown"]->Fill(s_JERDown.M());
				}
				else if(passT2_JER_Down){
					histograms["mjjcr1p_jerdown"]->Fill(s_JERDown.M());
				}
			}								
		} 
		else if(HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v==1){
			TLorentzVector j1;
			TLorentzVector j2;

			TLorentzVector j1_JESUp;
			TLorentzVector j2_JESUp;			
			TLorentzVector j1_JESDown;
			TLorentzVector j2_JESDown;
			TLorentzVector j1_JERUp;						
			TLorentzVector j2_JERUp;												
			TLorentzVector j1_JERDown;						
			TLorentzVector j2_JERDown;															

			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_JESUp.SetPtEtaPhiM(Jet_pt0_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_JESUp.SetPtEtaPhiM(Jet_pt1_Up, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_JESDown.SetPtEtaPhiM(Jet_pt0_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_JESDown.SetPtEtaPhiM(Jet_pt1_Down, Jet_eta[1], Jet_phi[1], Jet_mass[1]);
			j1_JERUp.SetPtEtaPhiM(Jet_pt0_JER_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_JERUp.SetPtEtaPhiM(Jet_pt1_JER_Up, Jet_eta[0], Jet_phi[0], Jet_mass[0]);			
			j1_JERDown.SetPtEtaPhiM(Jet_pt0_JER_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2_JERDown.SetPtEtaPhiM(Jet_pt1_JER_Down, Jet_eta[0], Jet_phi[0], Jet_mass[0]);							

			TLorentzVector s = j1+j2;
			TLorentzVector s_Up = j1_JESUp+j2_JESUp;
			TLorentzVector s_Down = j1_JESDown+j2_JESDown;
			TLorentzVector s_JERUp = j1_JERUp+j2_JERUp;
			TLorentzVector s_JERDown = j1_JERDown+j2_JERDown;			

			
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
			//=================  Plots Up ====================
			if(passT4_Up){
				histograms["mjjcr2_jesup"]->Fill(s_Up.M());					
			} 
			else{
				if(passT2_Up){
					histograms["mjjcr3_jesup"]->Fill(s_Up.M());						
				}
				if(passT3_Up){
					histograms["mjjcr2p_jesup"]->Fill(s_Up.M());
				}
				else if(passT2_Up){
					histograms["mjjcr3p_jesup"]->Fill(s_Up.M());
				}
			}

			//=================  Plots Down ====================
			if(passT4_Down){
				histograms["mjjcr2_jesdown"]->Fill(s_Down.M());					
			} 
			else{
				if(passT2_Down){
					histograms["mjjcr3_jesdown"]->Fill(s_Down.M());						
				}
				if(passT3_Down){
					histograms["mjjcr2p_jesdown"]->Fill(s_Down.M());
				}
				else if(passT2_Down){
					histograms["mjjcr3p_jesdown"]->Fill(s_Down.M());
				}
			}
			//=================  Plots JER Up ====================
			if(passT4_JER_Up){
				histograms["mjjcr2_jerup"]->Fill(s_JERUp.M());					
			} 
			else{
				if(passT2_JER_Up){
					histograms["mjjcr3_jerup"]->Fill(s_JERUp.M());						
				}
				if(passT3_JER_Up){
					histograms["mjjcr2p_jerup"]->Fill(s_JERUp.M());
				}
				else if(passT2_JER_Up){
					histograms["mjjcr3p_jerup"]->Fill(s_JERUp.M());
				}
			}			
			//=================  Plots JER Down ====================
			if(passT4_JER_Down){
				histograms["mjjcr2_jerdown"]->Fill(s_JERDown.M());					
			} 
			else{
				if(passT2_JER_Down){
					histograms["mjjcr3_jerdown"]->Fill(s_JERDown.M());						
				}
				if(passT3_JER_Down){
					histograms["mjjcr2p_jerdown"]->Fill(s_JERDown.M());
				}
				else if(passT2_JER_Down){
					histograms["mjjcr3p_jerdown"]->Fill(s_JERDown.M());
				}
			}

		}// HLT Quad45

	}//loop over entries

	Int_t nevts = 0;
	nevts = t->GetEntries();

	cout << "Events passing selection: " << nPass << endl;
	cout << "Aceptance: " << nPass*100./nevts << " +- " << 100.*TMath::Sqrt(nPass*(1.0 - (nPass/nevts)))/nevts << endl;
	cout << "JES Up: " << nPass_JESUp*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_JESUp*(1.0 - (nPass_JESUp/nevts)))/nevts << endl;
	cout << "JES Down: " << nPass_JESDown*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_JESDown*(1.0 - (nPass_JESDown/nevts)))/nevts << endl;
	cout << "JES shift up(%): " << 100.*((nPass*100./nevts)-(nPass_JESUp*100./nevts))/(nPass*100./nevts) << endl;
	cout << "JES shift down(%): " << 100.*((nPass*100./nevts)-(nPass_JESDown*100./nevts))/(nPass*100./nevts) << endl;	
	cout << "JER Up: " << nPass_JERUp*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_JERUp*(1.0 - (nPass_JERUp/nevts)))/nevts << endl;
	cout << "JER Down: " << nPass_JERDown*100./nevts << " +- " << 100.*TMath::Sqrt(nPass_JERDown*(1.0 - (nPass_JERDown/nevts)))/nevts << endl;	
	cout << "JER shift up(%): " << 100.*((nPass*100./nevts)-(nPass_JERUp*100./nevts))/(nPass*100./nevts) << endl;
	cout << "JER shift down(%): " << 100.*((nPass*100./nevts)-(nPass_JERDown*100./nevts))/(nPass*100./nevts) << endl;	

	//Save invariant mass histogram to a file
	TFile* outfile = new TFile("/fdata/hepx/store/user/delgado_andrea/mjjShapes/JESJER/mjj_DenisSel_JESJER_"+name+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;
}
