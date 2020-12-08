//BFF 4b channel trigger reweight fitting tool
//Author: Denis Rathjens
//Date of creation: 02/11/2020

#include <TROOT.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TBinomialEfficiencyFitter.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <vector>
#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <cmath>
#include <TMath.h>

TF1 *central, *uncertainty;
double Add(double *x, double *par){
  double result=central->EvalPar(x,par)+uncertainty->EvalPar(x,par);
  if(result<=1.) return result;
  else return 1.;
}

double Substract(double *x, double *par){
  double result=central->EvalPar(x,par)-uncertainty->EvalPar(x,par);
  if(result>=0.) return result;
  else return 0.;
}

void TriggerFits4b(){
  TFile *infile=new TFile("weightHistosAll.root","READ");
  //histograms to be read and fitted with error functions
  vector<TString> histogramsErf={"L1Quad30_caloT","Quad30Calo_caloT","Di90Calo_caloT","Quad30PF_PfT","Di90PF_PfT"};
  //starting parameters for each fit, so convergence works
  vector<double> p0Erf={215., 30., 90., 30., 88.};
  vector<double> p1Erf={ 44.,  3.,  1.,  1.,  5.};
  vector<double> p2Erf={  0.,0.05,0.01,0.05, 0.1};
 
  //savefile
  TFile *outfile=new TFile("ConfInterval.root","RECREATE");

  //loop over error function fits
  for(unsigned i=0; i<histogramsErf.size(); ++i){
    //fit definition and parameter setting
    central=new TF1("central_"+histogramsErf[i],"TMath::Erf((x-[0])/[1])/(2.+[2])+1./(2.+[2])",0.,1000.);
    central->SetParameters(p0Erf[i],p1Erf[i],p2Erf[i]);

    //converting the TEfficiency into something ROOT can get the covariance matrix out of
    TEfficiency *EffIn=(TEfficiency*) infile->Get(histogramsErf[i]);
    TH1F *TotalHist=(TH1F*) EffIn->GetTotalHistogram();
    TH1F *PassedHist=(TH1F*) EffIn->GetPassedHistogram();
    TBinomialEfficiencyFitter *fitter=new TBinomialEfficiencyFitter(PassedHist,TotalHist);
    //fit and covariance extraction
    fitter->Fit(central,"S");
    ROOT::Fit::Fitter *_fitter=fitter->GetFitter();
    ROOT::Fit::FitResult Result=_fitter->Result();
    float s00=Result.CovMatrix(0,0);
    float s10=Result.CovMatrix(1,0);
    float s01=Result.CovMatrix(0,1);
    float s11=Result.CovMatrix(1,1);
    float s20=Result.CovMatrix(2,0);
    float s02=Result.CovMatrix(0,2);
    float s12=Result.CovMatrix(2,1);
    float s21=Result.CovMatrix(1,2);
    float s22=Result.CovMatrix(2,2);
  
    //uncertainty function (symmetric)
    uncertainty=new TF1("uncertainty_"+histogramsErf[i],"TMath::Sqrt( ([0] -([6]-x)/[7]*[1] +pow(([6]-x)/[7],2)*[2])*pow(2*TMath::Exp(-pow(([6]-x)/[7],2))/TMath::Sqrt(TMath::Pi())/[7]/([8]+2),2)+([3] +([6]-x)/[7]*[4])*2*TMath::Exp(-pow(([6]-x)/[7],2))/TMath::Sqrt(TMath::Pi())/[7]/([8]+2) *(TMath::Erf([6]-x/[7])-1)/pow([8]+2,2)+ [5]*pow((TMath::Erf([6]-x/[7])-1)/pow([8]+2,2),2))",0.,1000.);
    uncertainty->SetLineColor(3);
    uncertainty->SetParameters(s00, s01+s10, s11, -(s02+s20), s12+s21, s22, central->GetParameter(0), central->GetParameter(1), central->GetParameter(2));

    //define upper and lower variations
    TF1 *upper=new TF1("upper_"+histogramsErf[i],Add,0.,1000.);
    upper->SetLineColor(4);
    TF1 *lower=new TF1("lower_"+histogramsErf[i],Substract,0.,1000.);
    lower->SetLineColor(8);

    //explicit test function
    for(float i=10.; i<300.; i+=5.) std::cout<<i<<": "<<central->Eval(i)<<"+/-"<<uncertainty->Eval(i)<<" = "<<upper->Eval(i)<<" to "<<lower->Eval(i)<<std::endl;

    //drawing and writing
    TCanvas *Canv=new TCanvas("Fit_"+histogramsErf[i],histogramsErf[i],500,500);
    EffIn->Draw();
    central->Draw("same");
    uncertainty->Draw("same");
    upper->Draw("same");
    lower->Draw("same");
    Canv->Write();
    Canv->SaveAs(histogramsErf[i]+".pdf");
    central->Write();
    uncertainty->Write();
    upper->Write();
    lower->Write();
    EffIn->Write();
  }

  //linear fit section
  TEfficiency *EffIn=(TEfficiency*) infile->Get("Di90Di30LeadCSV3_min_PF");
  TH1F *TotalHist=(TH1F*) EffIn->GetTotalHistogram();
  TH1F *PassedHist=(TH1F*) EffIn->GetPassedHistogram();
  TBinomialEfficiencyFitter *fitter=new TBinomialEfficiencyFitter(PassedHist,TotalHist);
  central=new TF1("central_Di90Di30LeadCSV3_min_PF","[0]*x+[1]",0.,1.);
  fitter->Fit(central,"S");
  ROOT::Fit::Fitter *_fitter=fitter->GetFitter();
  ROOT::Fit::FitResult Result=_fitter->Result();
  float s00=Result.CovMatrix(0,0);
  float s10=Result.CovMatrix(1,0);
  float s01=Result.CovMatrix(0,1);
  float s11=Result.CovMatrix(1,1);
  uncertainty=new TF1("uncertainty_Di90Di30LeadCSV3_min_PF","TMath::Sqrt([0]*x*x+[1]*x+[2])",0.,1.);
  uncertainty->SetParameters(s00,s10+s01,s11);
  uncertainty->SetLineColor(3);
  TF1 *upper=new TF1("upper_Di90Di30LeadCSV3_min_PF",Add,0.,1.);
  upper->SetLineColor(4);
  TF1 *lower=new TF1("lower_Di90Di30LeadCSV3_min_PF",Substract,0.,1.);
  lower->SetLineColor(8);
  TCanvas *Canv=new TCanvas("Fit_Di90Di30LeadCSV3_min_PF","Di90Di30LeadCSV3_min_PF",500,500);
  EffIn->Draw();
  central->Draw("same");
  uncertainty->Draw("same");
  upper->Draw("same");
  lower->Draw("same");
  Canv->Write();
  central->Write();
  uncertainty->Write();
  upper->Write();
  lower->Write();
  EffIn->Write();

  //file closing and cleanup
  infile->Close();
  outfile->Close();
}
