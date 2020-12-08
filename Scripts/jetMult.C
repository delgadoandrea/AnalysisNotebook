//Script that studies changes in jet multiplicity shapes after PU weight up/down variations

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

void jetMult(TString fname, TString outname){
	TFile *f = new TFile(fname, "READ");
	TTree *t = (TTree*) f->Get("tree");

	Int_t nJet;
	Float_t puWeight;
	Float_t puWeightUp;
	Float_t puWeightDown;

	t->SetBranchAddress("puWeight", &puWeight);
	t->SetBranchAddress("puWeightUp", &puWeightUp);
	t->SetBranchAddress("puWeightDown", &puWeightDown);
	t->SetBranchAddress("nJet",&nJet);	

	histograms["njet_nom"] = new TH1F("njet_nom", "Jet multiplicity", 10, 0, 10);
	histograms["njet_pu"] = new TH1F("njet_pu", "Jet multiplicity + PU weight", 10, 0, 10);
	histograms["njet_puup"] = new TH1F("njet_puup", "Jet multiplicity + PU weight up", 10, 0, 10);
	histograms["njet_pudown"] = new TH1F("njet_pudown", "Jet multiplicity + PU weight down", 10, 0, 10);

	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Sumw2();
	}

	cout<< "Analyzing sample: " << outname << endl;
	cout<< "Number of events in t: " << t->GetEntries() << endl;

	for (Int_t indx=0;indx<t->GetEntries();indx++){ 			
	//for (Int_t indx=0;indx<10;indx++){ 
		t->GetEntry(indx);
		histograms["njet_nom"]->Fill(nJet);
		histograms["njet_pu"]->Fill(nJet, puWeight);
		histograms["njet_puup"]->Fill(nJet, puWeightUp);
		histograms["njet_pudown"]->Fill(nJet, puWeightDown);

	}

	TFile* outfile = new TFile("../data/jetmult/njet_"+outname+".root","RECREATE");
	for(std::map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); ++it){
		histograms[it->first]->Write();	
	}
	outfile->Close();
	delete outfile;
}