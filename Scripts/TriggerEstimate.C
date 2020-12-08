#define TriggerEstimate_cxx
#include "TriggerEstimate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
#include <TH1.h>

#include <vector>
#include <map>
#include <iterator>
#include <iostream>

//--------------------------------------------------------------------
int Match(vector<TLorentzVector> ToMatch, TLorentzVector comparison){
  float smallestDeltaR=0.1;
  int Matched=-1;
  for(unsigned i=0; i<ToMatch.size(); ++i){
    const float DeltaR=comparison.DeltaR(ToMatch[i]);
    if(DeltaR<smallestDeltaR){
      smallestDeltaR=DeltaR;
      Matched=i;
    }
  }
  return Matched;
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
  weight = passed->GetBinContent(rBin)/total->GetBinContent(rBin);

  delete total;
  delete passed;
  return weight;
}
//--------------------------------------------------------------------

map<string, TEfficiency*> histograms;
map <string, TH1*> histogramsVal;

void TriggerEstimate::Loop(){
  if (fChain == 0) return;
  bool applyWeights = true;

  fChain->SetBranchStatus("*",0);  // disable all branches
   //jet block
  fChain->SetBranchStatus("Jet_pt",1);  // activate branchname
  fChain->SetBranchStatus("Jet_eta",1);  // activate branchname
  fChain->SetBranchStatus("Jet_phi",1);  // activate branchname
  fChain->SetBranchStatus("Jet_mass",1);  // activate branchname
  fChain->SetBranchStatus("Jet_btagCSV",1);  // activate branchname
  fChain->SetBranchStatus("Jet_btagCMVA",1);  // activate branchname   
  fChain->SetBranchStatus("nJet",1);  // activate branchname
  fChain->SetBranchStatus("Jet_btagDeepCSVb",1);
  fChain->SetBranchStatus("Jet_btagDeepCSVbb",1);
  fChain->SetBranchStatus("Jet_chHEF",1);
  fChain->SetBranchStatus("Jet_neHEF",1);
  fChain->SetBranchStatus("Jet_chEmEF",1);
  fChain->SetBranchStatus("Jet_neEmEF",1);
  fChain->SetBranchStatus("Jet_numberOfDaughters",1);
  fChain->SetBranchStatus("Jet_muEF",1);
  fChain->SetBranchStatus("Jet_chMult",1);
  fChain->SetBranchStatus("Jet_nhMult",1);
   //Quad45 block
  fChain->SetBranchStatus("ntrgObjects_hltL1sQuadJetCIorTripleJetVBFIorHTT",1);
  fChain->SetBranchStatus("trgObjects_hltL1sQuadJetCIorTripleJetVBFIorHTT_pt",1);
  fChain->SetBranchStatus("trgObjects_hltL1sQuadJetCIorTripleJetVBFIorHTT_eta",1);
  fChain->SetBranchStatus("trgObjects_hltL1sQuadJetCIorTripleJetVBFIorHTT_phi",1);
  fChain->SetBranchStatus("trgObjects_hltL1sQuadJetCIorTripleJetVBFIorHTT_mass",1);

  fChain->SetBranchStatus("ntrgObjects_hltQuadCentralJet45",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet45_pt",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet45_eta",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet45_phi",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet45_mass",1);

  fChain->SetBranchStatus("ntrgObjects_hltQuadPFCentralJetLooseID45",1);
  fChain->SetBranchStatus("trgObjects_hltQuadPFCentralJetLooseID45_pt",1);
  fChain->SetBranchStatus("trgObjects_hltQuadPFCentralJetLooseID45_eta",1);
  fChain->SetBranchStatus("trgObjects_hltQuadPFCentralJetLooseID45_phi",1);
  fChain->SetBranchStatus("trgObjects_hltQuadPFCentralJetLooseID45_mass",1);

  fChain->SetBranchStatus("HLT_BIT_HLT_QuadJet45_TripleBTagCSV_p087_v",1);
  fChain->SetBranchStatus("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV_p087_v",1);
   //Di90Di30 block
  fChain->SetBranchStatus("ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet",1);
  fChain->SetBranchStatus("ntrgObjects_hltQuadCentralJet30",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet30_pt",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet30_eta",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet30_phi",1);
  fChain->SetBranchStatus("trgObjects_hltQuadCentralJet30_mass",1);

  fChain->SetBranchStatus("ntrgObjects_hltQuadPFCentralJetLooseID30",1);
  fChain->SetBranchStatus("trgObjects_hltQuadPFCentralJetLooseID30_pt",1);
  fChain->SetBranchStatus("ntrgObjects_hltDoubleCentralJet90",1);
  fChain->SetBranchStatus("trgObjects_hltDoubleCentralJet90_pt",1);
  fChain->SetBranchStatus("ntrgObjects_hltDoublePFCentralJetLooseID90",1);
  fChain->SetBranchStatus("trgObjects_hltDoublePFCentralJetLooseID90_pt",1);
  fChain->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v",1);  // activate branchname
  fChain->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v",1);  // activate branchname
   //general trigger block
  fChain->SetBranchStatus("ntrgObjects_caloJets",1);
  fChain->SetBranchStatus("trgObjects_caloJets_pt",1);

  fChain->SetBranchStatus("ntrgObjects_pfJets",1);
  fChain->SetBranchStatus("trgObjects_pfJets_pt",1);

  fChain->SetBranchStatus("HLT_BIT_HLT_IsoTkMu27_v",1);
  fChain->SetBranchStatus("HLT_BIT_HLT_IsoMu24_v",1);
  fChain->SetBranchStatus("ntrgObjects_hltBTagCaloCSVp087Triple",1);
  fChain->SetBranchStatus("trgObjects_hltBTagCaloCSVp087Triple_pt",1);
  fChain->SetBranchStatus("trgObjects_hltBTagCaloCSVp087Triple_eta",1);
  fChain->SetBranchStatus("trgObjects_hltBTagCaloCSVp087Triple_phi",1);
  fChain->SetBranchStatus("trgObjects_hltBTagCaloCSVp087Triple_mass",1);
   //lepton block
  fChain->SetBranchStatus("nvLeptons",1);
  fChain->SetBranchStatus("vLeptons_isGlobalMuon",1);
  fChain->SetBranchStatus("vLeptons_mediumMuonId",1);
  fChain->SetBranchStatus("vLeptons_pt",1);
  fChain->SetBranchStatus("vLeptons_eta",1);
  fChain->SetBranchStatus("vLeptons_phi",1);
  fChain->SetBranchStatus("vLeptons_mass",1);
  fChain->SetBranchStatus("nselLeptons",1);
  fChain->SetBranchStatus("selLeptons_pt",1);         
  fChain->SetBranchStatus("selLeptons_eta",1);         
  fChain->SetBranchStatus("selLeptons_phi",1);         
  fChain->SetBranchStatus("selLeptons_mass",1);
  fChain->SetBranchStatus("selLeptons_isGlobalMuon",1);
  fChain->SetBranchStatus("selLeptons_mediumMuonId",1); 

  //----------------------------Define histograms here  
  TFile *w = TFile::Open("/fdata/hepx/store/user/delgado_andrea/94X/updatedNtuplesJune19/weightHistos_B_v2.root");

  TEfficiency *l1 = (TEfficiency*) w->Get("L1Quad30_caloT");
  TEfficiency *t1 = (TEfficiency*) w->Get("Quad30Calo_caloT");
  TEfficiency *t2 = (TEfficiency*) w->Get("Di90Calo_caloT");
  TEfficiency *t3 = (TEfficiency*) w->Get("Di90Di30LeadCSV3_min_PF");
  TEfficiency *t4 = (TEfficiency*) w->Get("Quad30PF_PfT");  
  TEfficiency *t5 = (TEfficiency*) w->Get("Di90PF_PfT");

  Double_t xbins[34] = {100., 120., 140., 160., 180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 425., 450., 475., 500., 525., 550., 575.,600., 650., 700., 750., 800., 850., 900., 950.,1000., 1100., 1200.};

  if(!applyWeights){                           
    //TH2F *necessaryNumberD90=new TH2F("necessaryNumberD90","trigger object number for fired HLT D90;trigger objects;subtrigger index", 11, -0.5, 10.5, 7, -0.5, 6.5);
    histograms["CutflowDi90Di30"] = new TEfficiency("CutflowDi90Di30","IsoMu24, GlobalMuon pT>54, 4pF jets>30, L1, Calo4, PF4, Calo2, PF2, CSV;category;acceptance",9, -0.5,8.5);
    histograms["L1Quad30_caloT"] = new TEfficiency("L1Quad30_caloT","L1 turnon Quad30 calo trg objects;#Sigma p_{T} 0-3 [GeV];efficiency",100,100.,600.); 
    histograms["L1Quad30_PfT"] = new TEfficiency("L1Quad30_PfT","L1 turnon Quad30 PF trg objects;#Sigma p_{T} 0-3 [GeV];efficiency",100,100.,600.);   
    histograms["L1Quad30_PF"] = new TEfficiency("L1Quad30_PF","L1 turnon Quad30 PF;#Sigma p_{T} 0-3 [GeV];efficiency",100,100.,600.);

    histograms["Quad30Calo_caloT"] = new TEfficiency("Quad30Calo_caloT","HLT Calo30 turnon calo trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);  
    histograms["Quad30Calo_PfT"] = new TEfficiency("Quad30Calo_PfT","HLT Calo30 turnon PF trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);   
    histograms["Quad30Calo_PF"] = new TEfficiency("Quad30Calo_PF","HLT Calo30 turnon unmatched PF;p_{T} 4 [GeV];efficiency",60,20.,140.);     

    histograms["Quad30PF_PfT"] = new TEfficiency("Quad30PF_PfT","HLT PF30 turnon PF trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);     
    histograms["Quad30PF_PF"] = new TEfficiency("Quad30PF_PF","HLT PF30 turnon unmatched PF;p_{T} 4 [GeV];efficiency",60,20.,140.);  

    histograms["Di90Calo_caloT"] = new TEfficiency("Di90Calo_caloT","HLT Calo30 turnon calo trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);  
    histograms["Di90Calo_PfT"]=new TEfficiency("Di90Calo_PfT","HLT Calo30 turnon PF trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);     
    histograms["Di90Calo_PF"]=new TEfficiency("Di90Calo_PF","HLT Calo30 turnon unmatched PF;p_{T} 4 [GeV];efficiency",60,20.,140.);     

    histograms["Di90PF_PfT"]=new TEfficiency("Di90PF_PfT","HLT PF30 turnon PF trg objects;p_{T} 4 [GeV];efficiency",60,20.,140.);       
    histograms["Di90PF_PF"]=new TEfficiency("Di90PF_PF","HLT PF30 turnon unmatched PF;p_{T} 4 [GeV];efficiency",60,20.,140.);       

    histograms["Di90Di30LeadCSV3_max_mTrgPF"] = new TEfficiency("Di90Di30LeadCSV3_max_mTrgPF","HLT leading CSV turnon trg match to PF;max CSV discriminator value;efficiency",100.,0.,1.);  
    histograms["Di90Di30LeadCSV3_max_PF"] = new TEfficiency("Di90Di30LeadCSV3_max_PF","HLT leading CSV turnon unmatched PF;max CSV discriminator value;efficiency",100.,0.,1.);   
    histograms["Di90Di30LeadCSV3_min_mTrgPF"] = new TEfficiency("Di90Di30LeadCSV3_min_mTrgPF","HLT leading CSV turnon trg match to PF;min_3 CSV discriminator value;efficiency",100.,0.,1.);
    histograms["Di90Di30LeadCSV3_min_PF"] = new TEfficiency("Di90Di30LeadCSV3_min_PF","HLT leading CSV turnon unmatched PF;min_3 CSV discriminator value;efficiency",100.,0.,1.);   
    histograms["Di90Di30LeadCSV3_min2_mTrgPF"] = new TEfficiency("Di90Di30LeadCSV3_min2_mTrgPF","HLT leading CSV turnon trg match to PF;min_2 CSV discriminator value;efficiency",100.,0.,1.);
    histograms["Di90Di30LeadCSV3_min2_PF"] = new TEfficiency("Di90Di30LeadCSV3_min2_PF","HLT leading CSV turnon unmatched PF;min_2 CSV discriminator value;efficiency",100.,0.,1.);   
  }
  else{

    histogramsVal["trigEta1"] = new TH1D("triggeredEventsvsEta1", "Triggered Events vs Eta of leading jet", 40, -4.0, 4.0);
    histogramsVal["trigEta2"] = new TH1D("triggeredEventsvsEta2", "Triggered Events vs Eta of 2nd leading jet", 40, -4.0, 4.0);  
    histogramsVal["trigEta3"] = new TH1D("triggeredEventsvsEta3", "Triggered Events vs Eta of PF selection leading jet", 40, -4.0, 4.0);   
    histogramsVal["weigthEta1"] = new TH1D("weightedEventsvsEta1", "Weighted Events vs Eta of leading jet", 40, -4.0, 4.0);  
    histogramsVal["weigthEta2"] = new TH1D("weightedEventsvsEta2", "Weighted Events vs Eta of 2nd leading jet", 40, -4.0, 4.0);    
    histogramsVal["weigthEta3"] = new TH1D("weightedEventsvsEta3", "Weighted Events vs Eta of PF selection leading jet", 40, -4.0, 4.0); 

    histogramsVal["trigPt1"] = new TH1D("triggeredEventsvsPt1", "Triggered Events vs Pt of leading jet", 160, 0, 800);
    histogramsVal["trigPt2"] = new TH1D("triggeredEventsvsPt2", "Triggered Events vs Pt of 2nd leading jet", 160, 0, 800);
    histogramsVal["trigPt3"] = new TH1D("triggeredEventsvsPt3", "Triggered Events vs Pt of PF selection  4th leading jet", 40, 0, 200);  
    histogramsVal["trigPt4"] = new TH1D("triggeredEventsvsPt4", "Triggered Events vs Pt of 4th leading jet", 40, 0, 200);  
    
    histogramsVal["weightPt1"] = new TH1D("weightededEventsvsPt1", "Weighted Events vs Pt of leading jet", 160, 0, 800);  
    histogramsVal["weightPt2"] = new TH1D("weightededEventsvsPt2", "Weighted Events vs Pt of 2nd leading jet", 160, 0, 800);  
    histogramsVal["weightPt3"] = new TH1D("weightededEventsvsPt3", "Weighted Events vs Pt of PF selection 4th leading jet", 40, 0, 200); 
    histogramsVal["weightPt4"] = new TH1D("weightededEventsvsPt4", "Weighted Events vs Pt of 4th leading jet", 40, 0, 200);  

    histogramsVal["trigCSV1"] = new TH1D("triggeredEventsvsCSV1", "Triggered Events vs CSV of leading jet", 20, 0, 1);
    histogramsVal["trigCSV2"] = new TH1D("triggeredEventsvsCSV2", "Triggered Events vs CSV of 2nd leading jet", 20, 0, 1);
    histogramsVal["trigCSV3"] = new TH1D("triggeredEventsvsCSV3", "Triggered Events vs CSV of PF selection 3rd leading jet", 20, 0, 1);
    histogramsVal["trigCSV4"] = new TH1D("triggeredEventsvsCSV4", "Triggered Events vs CSV of 3rd leading jet", 20, 0, 1);
        
    histogramsVal["weightCSV1"] = new TH1D("weightedEventsvsCSV1", "Weighted Events vs CSV of leading jet", 20, 0, 1); 
    histogramsVal["weightCSV2"] = new TH1D("weightedEventsvsCSV2", "Weighted Events vs CSV of 2nd leading jet", 20, 0, 1); 
    histogramsVal["weightCSV3"] = new TH1D("weightedEventsvsCSV3", "Weighted Events vs CSV of PF selection 3rd leading jet", 20, 0, 1);  
    histogramsVal["weightCSV4"] = new TH1D("weightedEventsvsCSV4", "Weighted Events vs CSV of 3rd leading jet", 20, 0, 1);       

    /*histogramsVal["trigmjjc1"] = new TH1D("trigmjjc1", "Triggered Events vs mjj in 2T", 120, 0, 1200);
    histogramsVal["trigmjjc2"] = new TH1D("trigmjjc2", "Triggered Events vs mjj in 2T!(3T||4T)", 120, 0, 1200);
    histogramsVal["trigmjjc3"] = new TH1D("trigmjjc3", "Triggered Events vs mjj in 3T!4T", 120, 0, 1200);
    histogramsVal["trigmjjc4"] = new TH1D("trigmjjc4", "Triggered Events vs mjj in 4T ", 120, 0, 1200);

    histogramsVal["weightmjjc1"] = new TH1D("weightmjjc1", "Weighted Events vs mjj in 2T", 120, 0, 1200);
    histogramsVal["weightmjjc2"] = new TH1D("weightmjjc2", "Weighted Events vs mjj in 2T!(3T||4T)", 120, 0, 1200);
    histogramsVal["weightmjjc3"] = new TH1D("weightmjjc3", "Weighted Events vs mjj in 3T!4T", 120, 0, 1200);
    histogramsVal["weightmjjc4"] = new TH1D("weightmjjc4", "Weighted Events vs mjj in 4T ", 120, 0, 1200); */

    histogramsVal["trigmjjc1"] = new TH1D("trigmjjc1", "Triggered Events vs mjj in 2T", 33, xbins);
    histogramsVal["trigmjjc2"] = new TH1D("trigmjjc2", "Triggered Events vs mjj in 2T!(3T||4T)", 33, xbins);
    histogramsVal["trigmjjc3"] = new TH1D("trigmjjc3", "Triggered Events vs mjj in 3T!4T", 33, xbins);
    histogramsVal["trigmjjc4"] = new TH1D("trigmjjc4", "Triggered Events vs mjj in 4T ", 33, xbins);

    histogramsVal["weightmjjc1"] = new TH1D("weightmjjc1", "Weighted Events vs mjj in 2T", 33, xbins);
    histogramsVal["weightmjjc2"] = new TH1D("weightmjjc2", "Weighted Events vs mjj in 2T!(3T||4T)", 33, xbins);
    histogramsVal["weightmjjc3"] = new TH1D("weightmjjc3", "Weighted Events vs mjj in 3T!4T", 33, xbins);
    histogramsVal["weightmjjc4"] = new TH1D("weightmjjc4", "Weighted Events vs mjj in 4T ", 33, xbins); 

  }

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries2 = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<100;jentry++) {    
//    if(jentry%10==0) cout<<jentry<<"/"<<nentries2<<endl;
    if(jentry%100000==0) cout<<jentry<<"/"<<nentries2<<endl;    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //-------------------------verbose switch/bool initialization
    bool verbose=false;

    bool passT2 = false; 
    bool passT3 = false;                
    bool passT4 = false;                

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

    //--------------------------muon good event selection simple
    if(!applyWeights) histograms["CutflowDi90Di30"]->Fill(HLT_BIT_HLT_IsoMu24_v,0);
    if(!HLT_BIT_HLT_IsoMu24_v) continue;            
    if(verbose) cout<<"passed IsoMu24"<<endl;
    bool passMuon=false;

    for(unsigned l=0; l<nselLeptons; ++l){
      if(!(selLeptons_isGlobalMuon[l] && selLeptons_mediumMuonId[l])) continue;
      if(selLeptons_pt[l]>54) passMuon=true;
    }
      
    if(!applyWeights) histograms["CutflowDi90Di30"]->Fill(passMuon,1);
    if(!passMuon) continue;      
    if(verbose) cout<<"passed muon requirement"<<endl;

    //----------------------------require 4 PF jets >30
    float MaxBtag=0.;
    float MinBtag=1.;
    float PreMinBtag=1.;
    vector<TLorentzVector> PFjets;
    vector<float> CSV;
    for(unsigned j=0; j<nJet; ++j){
      if(Jet_pt[j]<30 || fabs(Jet_eta[j])>2.6) continue;
      TLorentzVector jet;
      const float DeepCSV=Jet_btagDeepCSVb[j]+Jet_btagDeepCSVbb[j];
      jet.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],Jet_mass[j]);
      if(fabs(Jet_eta[j])<2.7){
        if(Jet_neHEF[j]<0.9 && Jet_neHEF[j]<0.9 && Jet_numberOfDaughters[j]>1 && Jet_muEF[j]<0.8){
          if(fabs(Jet_eta[j])<2.4){
            if(Jet_chHEF[j]>0 && Jet_chMult[j]>0 && Jet_chEmEF[j] < 0.9){
              if(PFjets.size()<4 && DeepCSV>MaxBtag) MaxBtag=DeepCSV;                     
              PFjets.push_back(jet);
              CSV.push_back(DeepCSV);                     
            }
          }
          else{ PFjets.push_back(jet);CSV.push_back(DeepCSV);}               
        }
      }//else if eta < 2.7
      else if(fabs(Jet_eta[j])<3.0){
        if(Jet_neEmEF[j]>0.01 && Jet_neHEF[j]<0.98 && Jet_nhMult[j]>2) {PFjets.push_back(jet);CSV.push_back(DeepCSV);}            
      }
      else if(fabs(Jet_eta[j])>=3.0 && Jet_neEmEF[j]<0.9 && Jet_nhMult[j]>10) {PFjets.push_back(jet);CSV.push_back(DeepCSV);}         
    }//loop over jets

    //-----------------------------ensure the 3rd highest of the first 4 leading jets in btag value is chosen as minimum
    vector<float> CSV2=CSV;
    std::sort(CSV2.begin(),CSV2.end(),[](float a, float b){return a>b;});
    if(CSV2.size()){
      MinBtag=CSV2[min((int)CSV2.size()-1,2)];
      PreMinBtag=CSV2[min((int)CSV2.size()-1,1)];
    }

    //---------------------------- lepton-jet veto
    for(unsigned i=0; i<nselLeptons; i++){
      TLorentzVector lepton;
      lepton.SetPtEtaPhiM(selLeptons_pt[i], selLeptons_eta[i], selLeptons_phi[i], selLeptons_mass[i]);
      for(unsigned int j=0; j<PFjets.size(); j++){
        if(lepton.DeltaR(PFjets[j])<0.4) PFjets.erase(PFjets.begin() + j);    
      }  
    }
    if(verbose) cout<<"passed lepton-jet veto requirement"<<endl;
    //----------------------------jet selection test
    if(!applyWeights) histograms["CutflowDi90Di30"]->Fill(PFjets.size()>=4,2);
    if(PFjets.size()<4) continue;
    if(verbose) cout<<"passed jet requirement"<<endl;


    //--------------------------b-tagging selection
    if(Jet_pt[0]>30 && fabs(Jet_eta[0])<2.6){
      const float DeepCSV=Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0];
      if(DeepCSV>0.8958)       FiT=true;
      else if (DeepCSV>0.6324) FiM=true;
      else if (DeepCSV>0.2219) FiL=true;
    }
    if(Jet_pt[1]>30 && fabs(Jet_eta[1])<2.6){
      const float DeepCSV=Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1];
      if(DeepCSV>0.8958)       SeT=true;
      else if (DeepCSV>0.6324) SeM=true;
      else if (DeepCSV>0.2219) SeL=true;
    }

    if(Jet_pt[2]>30 && fabs(Jet_eta[2])<2.6){
      const float DeepCSV=Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2];
      if(DeepCSV>0.8958)       ThT=true;
      else if (DeepCSV>0.6324) ThM=true;
      else if (DeepCSV>0.2219) ThL=true;
    }
    if(Jet_pt[3]>30 && fabs(Jet_eta[3])<2.6){
      const float DeepCSV=Jet_btagDeepCSVb[3]+Jet_btagDeepCSVbb[3];
      if(DeepCSV>0.8958)       FoT=true;
      else if (DeepCSV>0.6324) FoM=true;
      else if (DeepCSV>0.2219) FoL=true;
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

    if(!passT2) continue;
    if(verbose) cout<<"passed T2 requirement"<<endl;
    //-------------------------------match PF jets to objects triggering CSVp087 trigger
    std::vector<float> CSVvals;
    for(unsigned i=0; i<ntrgObjects_hltBTagCaloCSVp087Triple; ++i){
      TLorentzVector bTrigger;
      bTrigger.SetPtEtaPhiM(trgObjects_hltBTagCaloCSVp087Triple_pt[i],trgObjects_hltBTagCaloCSVp087Triple_eta[i],trgObjects_hltBTagCaloCSVp087Triple_phi[i],trgObjects_hltBTagCaloCSVp087Triple_mass[i]);
      int indexB=Match(PFjets,bTrigger);
      if(indexB!=-1) CSVvals.push_back(CSV[indexB]);
    }
      
    std::sort(CSVvals.begin(),CSVvals.end(),[](float a, float b){return a>b;});
    if(verbose) cout<<"sorted CSV - matched"<<endl;
    //------------------------------D90D30 trigger
    float SumCalo=0.;
    for(unsigned c=0; c<min(ntrgObjects_caloJets,4); ++c) SumCalo+=trgObjects_caloJets_pt[c];
    
    if(verbose) cout<<"sumcalo computed"<<endl;
    
    if(applyWeights){
      Double_t totalWeight = 0;
      vector<double> weights;

      weights.push_back(weight(l1, SumCalo));
      weights.push_back(weight(t1, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]));
      weights.push_back(weight(t2, trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]));
      weights.push_back(weight(t3, MinBtag));
      weights.push_back(weight(t4, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]));
      weights.push_back(weight(t5, trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]));
    
      for (int i = 0; i < weights.size(); ++i){
        if(TMath::IsNaN(weights[i])) weights[i] = 0.0;
      }
      totalWeight = weights[0] *weights[1] *weights[2] *weights[3] *weights[4] *weights[5];
   
      TLorentzVector j1;
      TLorentzVector j2;
      j1.SetPtEtaPhiM(Jet_pt[0], Jet_eta[0], Jet_phi[0], Jet_mass[0]);
      j2.SetPtEtaPhiM(Jet_pt[1], Jet_eta[1], Jet_phi[1], Jet_mass[1]);  
      TLorentzVector s = j1+j2;

      if(HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v){
        histogramsVal["trigEta1"]->Fill(Jet_eta[0]);
        histogramsVal["trigEta2"]->Fill(Jet_eta[1]);
        histogramsVal["trigEta3"]->Fill(PFjets[0].Eta());

        histogramsVal["trigPt1"]->Fill(Jet_pt[0]);
        histogramsVal["trigPt2"]->Fill(Jet_pt[1]);
        histogramsVal["trigPt3"]->Fill(PFjets[3].Pt());
        histogramsVal["trigPt4"]->Fill(Jet_pt[3]);

        histogramsVal["trigCSV1"]->Fill(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0]);                            
        histogramsVal["trigCSV2"]->Fill(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1]);                            
        histogramsVal["trigCSV3"]->Fill(MinBtag);                             
        histogramsVal["trigCSV4"]->Fill(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2]);             

        histogramsVal["trigmjjc1"]->Fill(s.M());
        if(!passT3  && !passT4) histogramsVal["trigmjjc2"]->Fill(s.M());
        if(passT3 && !passT4)   histogramsVal["trigmjjc3"]->Fill(s.M());  
        if(passT4)              histogramsVal["trigmjjc4"]->Fill(s.M());                         
      }
      
      histogramsVal["weigthEta1"]->Fill(Jet_eta[0], totalWeight);
      histogramsVal["weigthEta2"]->Fill(Jet_eta[1], totalWeight);
      histogramsVal["weigthEta3"]->Fill(PFjets[0].Eta(), totalWeight);
                                    
      histogramsVal["weightPt1"]->Fill(Jet_pt[0], totalWeight);
      histogramsVal["weightPt2"]->Fill(Jet_pt[1], totalWeight);
      histogramsVal["weightPt3"]->Fill(PFjets[3].Pt(), totalWeight);
      histogramsVal["weightPt4"]->Fill(Jet_pt[3], totalWeight);

      histogramsVal["weightCSV1"]->Fill(Jet_btagDeepCSVb[0]+Jet_btagDeepCSVbb[0], totalWeight);
      histogramsVal["weightCSV2"]->Fill(Jet_btagDeepCSVb[1]+Jet_btagDeepCSVbb[1], totalWeight);
      histogramsVal["weightCSV3"]->Fill(MinBtag, totalWeight);
      histogramsVal["weightCSV4"]->Fill(Jet_btagDeepCSVb[2]+Jet_btagDeepCSVbb[2], totalWeight);            

      histogramsVal["weightmjjc1"]->Fill(s.M(),totalWeight);
      if(!passT3  && !passT4) histogramsVal["weightmjjc2"]->Fill(s.M(),totalWeight);
      if(passT3 && !passT4)   histogramsVal["weightmjjc3"]->Fill(s.M(),totalWeight);  
      if(passT4)              histogramsVal["weightmjjc4"]->Fill(s.M(),totalWeight);  
    }

    else{
      histograms["CutflowDi90Di30"]->Fill(ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet,3);
      histograms["L1Quad30_caloT"]->Fill(ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet>=1, SumCalo);
      
      SumCalo=0.;
      for(unsigned c=0; c<min(ntrgObjects_pfJets,4); ++c) SumCalo+=trgObjects_pfJets_pt[c];
      histograms["L1Quad30_PfT"]->Fill(ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet>=1, SumCalo); 

      SumCalo=0.;
      for(unsigned c=0; c<min((int)PFjets.size(),4); ++c) SumCalo+=PFjets[c].Pt();
      histograms["L1Quad30_PF"]->Fill(ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet>=1, SumCalo);
      
      if(verbose) cout<<"L1 histos filled"<<endl;
      
      if(ntrgObjects_hltL1sTripleJetVBFIorHTTIorDoubleJetCIorSingleJet>=1){
      //Quad30Calo
        histograms["Quad30Calo_caloT"]->Fill(ntrgObjects_hltQuadCentralJet30>=4,trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,3)]);
        histograms["Quad30Calo_PfT"]->Fill(ntrgObjects_hltQuadCentralJet30>=4,trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]);
        histograms["Quad30Calo_PF"]->Fill(ntrgObjects_hltQuadCentralJet30>=4,PFjets[3].Pt());
        if(verbose) cout<<"Quad30 histos filled"<<endl;
        
        if(ntrgObjects_hltQuadCentralJet30>=4){
        //Di90Calo
          histograms["Di90Calo_caloT"]->Fill(ntrgObjects_hltDoubleCentralJet90>=2,trgObjects_caloJets_pt[min(ntrgObjects_caloJets-1,1)]);
          histograms["Di90Calo_PfT"]->Fill(ntrgObjects_hltDoubleCentralJet90>=2,trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]);
          histograms["Di90Calo_PF"]->Fill(ntrgObjects_hltDoubleCentralJet90>=2,PFjets[1].Pt());
          if(verbose) cout<<"Di90 histos filled"<<endl;
          
          if(ntrgObjects_hltDoubleCentralJet90>=2){
            /*std::vector<float> CSVvals;
            for(unsigned i=0; i<ntrgObjects_hltBTagCaloCSVp087Triple; ++i){
              TLorentzVector bTrigger;
              bTrigger.SetPtEtaPhiM(trgObjects_hltBTagCaloCSVp087Triple_pt[i],trgObjects_hltBTagCaloCSVp087Triple_eta[i],trgObjects_hltBTagCaloCSVp087Triple_phi[i],trgObjects_hltBTagCaloCSVp087Triple_mass[i]);
              int indexB=Match(PFjets,bTrigger);
              if(indexB!=-1) CSVvals.push_back(CSV[indexB]);
            }
      
            std::sort(CSVvals.begin(),CSVvals.end(),[](float a, float b){return a>b;});*/

            if(CSVvals.size()){
              if(verbose) cout<<"CSVvalssize: " << CSVvals.size()<<endl;              
              histograms["Di90Di30LeadCSV3_max_mTrgPF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,CSVvals[0]);
              histograms["Di90Di30LeadCSV3_min_mTrgPF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,CSVvals[min(2,(int)CSVvals.size()-1)]);
              histograms["Di90Di30LeadCSV3_min2_mTrgPF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,CSVvals[min(1,(int)CSVvals.size()-1)]);
            }
      
            histograms["Di90Di30LeadCSV3_max_PF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,MaxBtag);
            histograms["Di90Di30LeadCSV3_min_PF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,MinBtag);
            histograms["Di90Di30LeadCSV3_min2_PF"]->Fill(ntrgObjects_hltBTagCaloCSVp087Triple>=3,PreMinBtag);
            if(verbose) cout<<" CSV histos filled"<<endl;
            if(ntrgObjects_hltBTagCaloCSVp087Triple>=3){
            //Quad30PF
              histograms["Quad30PF_PfT"]->Fill(ntrgObjects_hltQuadPFCentralJetLooseID30>=4,trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,3)]);
              histograms["Quad30PF_PF"]->Fill(ntrgObjects_hltQuadPFCentralJetLooseID30>=4,PFjets[3].Pt());
            if(verbose) cout<<"Quad30PF histos filled"<<endl;
              if(ntrgObjects_hltQuadPFCentralJetLooseID30>=4){
              //Di90PF
                if(verbose) cout<<"pt of trigger obj: "<<trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]<<endl;                
                if(verbose) cout<<"pt of PF: "<<PFjets[1].Pt()<<endl;  
                if(verbose) cout<<"ntrgObjects_hltDoublePFCentralJetLooseID90: "<<ntrgObjects_hltDoublePFCentralJetLooseID90<<endl;                                                              
                histograms["Di90PF_PfT"]->Fill(ntrgObjects_hltDoublePFCentralJetLooseID90>=2,trgObjects_pfJets_pt[min(ntrgObjects_pfJets-1,1)]);
                //histograms["Di90PF_Pf"]->Fill(ntrgObjects_hltDoublePFCentralJetLooseID90>=2,PFjets[1].Pt());
                if(verbose) cout<<"Di90 PF histos filled"<<endl;
              }
            }
          }
        }
      }
    }//Loop over !applyWeights
  }//Loop over entries

  TFile *savefile;

  if(applyWeights){
    savefile = new TFile("triggerClosure"+name+"_v2.root","RECREATE");
    for(std::map<string,TH1*>::iterator it=histogramsVal.begin(); it!=histogramsVal.end(); ++it){
      histogramsVal[it->first]->Write(); 
    }
  }   
  else{
    savefile = new TFile("weightHistos_"+name+"_v2.root","RECREATE");
    histograms["CutflowDi90Di30"]->Write();
    histograms["L1Quad30_caloT"]->Write(); 
    histograms["L1Quad30_PfT"]->Write(); 
    histograms["L1Quad30_PF"]->Write(); 
    histograms["Quad30Calo_caloT"]->Write();
    histograms["Quad30Calo_PfT"]->Write(); 
    histograms["Quad30Calo_PF"]->Write();
    histograms["Quad30PF_PfT"]->Write();
    histograms["Quad30PF_PF"]->Write();
    histograms["Di90Calo_caloT"]->Write();
    histograms["Di90Calo_PfT"]->Write();
    histograms["Di90Calo_PF"]->Write();
    histograms["Di90PF_PfT"]->Write();
    histograms["Di90PF_PF"]->Write();
    histograms["Di90Di30LeadCSV3_max_mTrgPF"]->Write();
    histograms["Di90Di30LeadCSV3_max_PF"]->Write();
    histograms["Di90Di30LeadCSV3_min_mTrgPF"]->Write();
    histograms["Di90Di30LeadCSV3_min_PF"]->Write();
    histograms["Di90Di30LeadCSV3_min2_mTrgPF"]->Write();
    histograms["Di90Di30LeadCSV3_min2_PF"]->Write(); 
  }
  
  savefile->Close();
}
