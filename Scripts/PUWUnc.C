//Code that calculates the variations in accceptance when shifting the PU weights up/down
//Takes as input the ntuples that had the pu Weights calculated from a second pass od input PU histos
//Also creates histograms to study differences in shape after weight variations
//------------------------------------------------------

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

map <string, TH1*> histograms;

void PUWUnc (TChain *chain, TString name){

	if (chain == 0 ){
		std::cout << "No tree found!" <<std::endl;
		return;
   	} 

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
	Float_t puWeight;
	Float_t puWeightUp;
	Float_t puWeightDown;

	chain->SetBranchAddress("puWeight", &puWeight);
	chain->SetBranchAddress("puWeightUp", &puWeightUp);
	chain->SetBranchAddress("puWeightDown", &puWeightDown);
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

	histograms["mjjsr_nom"] = new TH1F("mjjsrnom", "mjj from two leading jets in event for signal region ", 120,0, 1200);	
	histograms["mjjcr1_nom"] = new TH1F("mjjcr1nom", "mjj from two leading jets in event for control region 1", 120,0, 1200);	
	histograms["mjjcr2_nom"] = new TH1F("mjjcr2nom", "mjj from two leading jets in event for control region 2", 120,0, 1200);		
	histograms["mjjcr3_nom"] = new TH1F("mjjcr3nom", "mjj from two leading jets in event for control region 3", 120,0, 1200);	
	histograms["mjjsr"] = new TH1F("mjjsr", "mjj from two leading jets in event for signal region ", 120,0, 1200);	
	histograms["mjjcr1"] = new TH1F("mjjcr1", "mjj from two leading jets in event for control region 1", 120,0, 1200);	
	histograms["mjjcr2"] = new TH1F("mjjcr2", "mjj from two leading jets in event for control region 2", 120,0, 1200);		
	histograms["mjjcr3"] = new TH1F("mjjcr3", "mjj from two leading jets in event for control region 3", 120,0, 1200);			
	histograms["mjjcr0p"] = new TH1F("mjjcr0p", "mjj frommjjcr3 two leading jets in event for control region 0'", 120,0, 1200);	
	histograms["mjjcr1p"] = new TH1F("mjjcr1p", "mjj from two leading jets in event for control region 1'", 120,0, 1200);	
	histograms["mjjcr2p"] = new TH1F("mjjcr2p", "mjj from two leading jets in event for control region 2'", 120,0, 1200);		
	histograms["mjjcr3p"] = new TH1F("mjjcr3p", "mjj from two leading jets in event for control region 3'", 120,0, 1200);

	histograms["mjjsr_jup"] = new TH1F("mjjsr_jup", "mjj from two leading jets in event for signal region +7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1_jup"] = new TH1F("mjjcr1_jup", "mjj from two leading jets in event for control region 1 +7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2_jup"] = new TH1F("mjjcr2_jup", "mjj from two leading jets in event for control region 2 +7p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3_jup"] = new TH1F("mjjcr3_jup", "mjj from two leading jets in event for control region 3 +7p minbias xsec", 120,0, 1200);			
	histograms["mjjcr0p_jup"] = new TH1F("mjjcr0p_jup", "mjj frommjjcr3 two leading jets in event for control region 0' +7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1p_jup"] = new TH1F("mjjcr1p_jup", "mjj from two leading jets in event for control region 1' +7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2p_jup"] = new TH1F("mjjcr2p_jup", "mjj from two leading jets in event for control region 2' +7p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3p_jup"] = new TH1F("mjjcr3p_jup", "mjj from two leading jets in event for control region 3' +7p minbias xsec", 120,0, 1200);

	histograms["mjjsr_jdown"] = new TH1F("mjjsr_jdown", "mjj from two leading jets in event for signal region -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1_jdown"] = new TH1F("mjjcr1_jdown", "mjj from two leading jets in event for control region 1 -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2_jdown"] = new TH1F("mjjcr2_jdown", "mjj from two leading jets in event for control region 2 -7p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3_jdown"] = new TH1F("mjjcr3_jdown", "mjj from two leading jets in event for control region 3 -7p minbias xsec", 120,0, 1200);			
	histograms["mjjcr0p_jdown"] = new TH1F("mjjcr0p_jdown", "mjj frommjjcr3 two leading jets in event for control region 0' -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr1p_jdown"] = new TH1F("mjjcr1p_jdown", "mjj from two leading jets in event for control region 1' -7p minbias xsec", 120,0, 1200);	
	histograms["mjjcr2p_jdown"] = new TH1F("mjjcr2p_jdown", "mjj from two leading jets in event for control region 2' -7p minbias xsec", 120,0, 1200);		
	histograms["mjjcr3p_jdown"] = new TH1F("mjjcr3p_jdown", "mjj from two leading jets in event for control region 3' -7p minbias xsec", 120,0, 1200);

	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << name << endl;
	cout<< "Number of events in t: " << chain->GetEntries() << endl;
	Int_t nPass = 0;	
	Int_t nPass_up = 0;
	Int_t nPass_down = 0;
	Int_t nPass_nom = 0;

	Int_t tot = 0;
	Int_t tot_up = 0;
	Int_t tot_down = 0;			

for (Int_t indx=0;indx<chain->GetEntries();indx++){ 			
	//for (Int_t indx=0;indx<10;indx++){ 
		chain->GetEntry(indx);

		tot = tot + puWeight;
		tot_up = tot_up + puWeightUp;
		tot_down = tot_down + puWeightDown;

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

//============================= mjj Calculation ===============================
		
		if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v==1){						
			TLorentzVector j1;
			TLorentzVector j2;

			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);

			TLorentzVector s = j1+j2;

			//=================  Plots ====================
			if(passT4){
				histograms["mjjsr_nom"]->Fill(s.M());
				histograms["mjjsr"]->Fill(s.M(),puWeight);
				histograms["mjjsr_jup"]->Fill(s.M(), puWeightUp);
				histograms["mjjsr_jdown"]->Fill(s.M(), puWeightDown);
				nPass = nPass + puWeight;
				nPass_up = nPass_up + puWeightUp;
				nPass_down = nPass_down + puWeightDown;
				nPass_nom++;
			} 
			else{
				if(passT2){
					histograms["mjjcr1_nom"]->Fill(s.M());
					histograms["mjjcr1"]->Fill(s.M(), puWeight);
					histograms["mjjcr1_jup"]->Fill(s.M(), puWeightUp);
					histograms["mjjcr1_jdown"]->Fill(s.M(), puWeightDown);
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
			TLorentzVector j1;
			TLorentzVector j2;
			
			j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
			j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);

			TLorentzVector s = j1+j2;
			
			if(passT4){
				histograms["mjjcr2_nom"]->Fill(s.M());
				histograms["mjjcr2"]->Fill(s.M(),puWeight);
				histograms["mjjcr2_jdown"]->Fill(s.M(), puWeightDown);
				histograms["mjjcr2_jup"]->Fill(s.M(), puWeightUp);
				//histograms["test"]->Fill(indx,2.);					
			} 
			else{
				if(passT2){
					histograms["mjjcr3_nom"]->Fill(s.M());
					histograms["mjjcr3"]->Fill(s.M(),puWeight);
					histograms["mjjcr3_jdown"]->Fill(s.M(),puWeightDown);
					histograms["mjjcr3_jup"]->Fill(s.M(),puWeightUp);
					//histograms["test"]->Fill(indx,3.);						
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

	Double_t nom_acc = 0.;
	Double_t acc_up = 0.;
	Double_t acc_down = 0.;

	nom_acc = nPass*100./tot;
	acc_up = nPass_up*100./tot_up;
	acc_down = nPass_down*100./tot_down;

	std::cout << "Events passing selection: " << nPass_nom << " E: " << chain->GetEntries() <<" A: " << nPass_nom*100./chain->GetEntries() << " +- "<< 100.*TMath::Sqrt(nPass_nom*(1.0 - (nPass_nom/chain->GetEntries())))/chain->GetEntries() << endl;
	std::cout << "PU weight: " << nPass<< " E: " << tot << " A: " << nPass*100./tot << " +- "<< 100.*TMath::Sqrt(nPass*(1.0 - (nPass/tot)))/tot<< endl;	
	std::cout << "PU Up: " << nPass_up << " E: " << tot_up << " A: " << nPass_up*100./tot_up<< " +- "<<100.*TMath::Sqrt(nPass_up*(1.0 - (nPass_up/tot_up)))/tot_up<<endl;
	std::cout << "PU Down: " <<nPass_down << " E: " << tot_down << " A: " << nPass_down*100./tot_down<< " +- "<<100.*TMath::Sqrt(nPass_down*(1.0 - (nPass_down/tot_down)))/tot_down<<endl;
	std::cout << "% Variation up: " << 100.*(nom_acc - acc_up)/nom_acc << endl;
	std::cout << "% Variation down: " << 100.*(nom_acc - acc_down)/nom_acc << endl;

	//Save invariant mass histogram to a file
	TFile* outfile = new TFile("../data/mjj+PU/mjj_"+name+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;

}