#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

void njetPlot(TString fname, TString outname){

	auto c1 = new TCanvas("c1","c1",600,500);
   	gStyle->SetOptStat(0);

	TFile *f = new TFile(fname, "READ");
	TH1F *h = (TH1F*) f->Get("njet_nom");
	TH1F *hpu = (TH1F*) f->Get("njet_pu");
	TH1F *hup = (TH1F*) f->Get("njet_puup");
	TH1F *hdown = (TH1F*) f->Get("njet_pudown");

	h->SetMarkerStyle(20);
	h->SetLineColor(kRed);
	h->SetMarkerColor(kRed);
	hpu->SetMarkerStyle(20);
	hpu->SetLineColor(kBlack);
	hpu->SetMarkerColor(kBlack);
	hup->SetMarkerStyle(20);
	hup->SetLineColor(kGreen);
	hup->SetMarkerColor(kGreen);
	hdown->SetMarkerStyle(20);
	hdown->SetLineColor(kCyan);
	hdown->SetMarkerColor(kCyan);
	hup->SetTitle("Jet multiplicity (Normalized)");
	hup->GetXaxis()->SetTitle("nJet");
	hup->GetYaxis()->SetTitle("A.U.");	

	hup->DrawNormalized();
	hpu->DrawNormalized("same");
	hdown->DrawNormalized("same");
	//h->DrawNormalized("same");

   auto legend = new TLegend(0.0986622,0.765823,0.434783,0.900844);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   //legend->AddEntry(h,"No weight","apl");
   legend->AddEntry(hpu,"Nominal PU weight","apl");
   legend->AddEntry(hup,"PU weight up","apl");
   legend->AddEntry(hdown,"PU weight down","apl");
   legend->Draw("same");

   //c1->SaveAs("../plots/jetmult/njet_Zprime_m500.png");
   c1->SaveAs(outname);

}