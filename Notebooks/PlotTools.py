from ROOT import TCanvas, TFile, gRandom, TH1F, TColor, TGaxis, TPad, TLegend, TLine
from ROOT import kBlack, kBlue, kRed, kGreen
import ROOT

def createRatio(h1, h2, name):
    h3 = h1.Clone("h3")
    h3.SetLineColor(kBlack)
    h3.SetMarkerColor(kBlack)
    h3.SetMarkerStyle(21)
    h3.SetTitle("")
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)

    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle(name)
    y.SetNdivisions(505)
    y.SetTitleSize(18)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)

    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)

    return h3

def createCanvasPads():
    c = TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()

    return c, pad1, pad2

def make_plots(h1,h1_up, h1_down, name):
    h1.GetYaxis().SetTitle('Events')
    h1.GetXaxis().SetTitle('m_{jj}[GeV]')
    h1.GetYaxis().SetLabelSize(0.02)
    h1.GetYaxis().SetTitleOffset(1.2)
    h1.GetXaxis().SetLabelSize(0.02)
    h1.GetXaxis().SetTitleOffset(1.2)
    h1.SetStats(0)

    h1.SetMarkerStyle(21)
    h1_up.SetMarkerStyle(21)
    h1_down.SetMarkerStyle(21)
    h1_up.SetLineColor(kRed)
    h1_up.SetMarkerColor(kRed)
    h1_down.SetLineColor(kGreen)
    h1_down.SetMarkerColor(kGreen)
    h1.SetTitle(name)
    
    leg = TLegend(.53,.82,.87,.53)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(h1,"JEC","APL")
    leg.AddEntry(h1_up,"JEC +1#sigma","APL")
    leg.AddEntry(h1_down,"JEC -1#sigma","APL")
    
    return h1, h1_up, h1_down, leg