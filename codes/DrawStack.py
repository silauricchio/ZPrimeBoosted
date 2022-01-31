import ROOT
import matplotlib.pyplot as plt
import numpy as np

var_list = [
#"MET_Et", 
"ttbar_mass"
#"ttbar_ptimbalance"
]

masses = [
#"400","500","750","1000","1250", 
"1500"
#, "2500"
]
ZPsColors = [
ROOT.kBlue
#, ROOT.kAzure-2, ROOT.kViolet, ROOT.kMagenta, ROOT.kOrange-1, ROOT.kYellow-2
]

f = ROOT.TFile.Open("../run/file_histos.root")


# Get histos                                                                                                                                                                    

for var in var_list:                  
    print(var)

    data = f.Get(var + "_data")
    tt = f.Get(var + "_tt")
    Vjets = f.Get(var + "_Vjets")
    st = f.Get(var + "_single_top")
    dib = f.Get(var + "_Diboson")
    ZPs = []
    h_signals = []
    for mass in masses:
        h_signal = f.Get(var + "_ZPrime"+mass)
        h_signals.append(h_signal)
        ZPs.append(h_signal)
    
    
    # Set styles
    ROOT.gROOT.SetStyle("ATLAS")
    
    # Create canvas
    c = ROOT.TCanvas("c", "", 600, 600)
    c.SetTickx(0)
    c.SetTicky(0)
    #c.SetLogy()
    #c.SetLogx()
    
    stack = ROOT.THStack()
    for h, color in zip([dib, Vjets, st, tt],
            [(208, 240, 193), (155, 152, 204), (248, 206, 104), (222, 90, 106)]):
        h.SetLineWidth(1)
        h.SetLineColor(1)
        h.SetFillColor(ROOT.TColor.GetColor(*color))
        #h.GetXaxis().SetLimits(0,2000)
        stack.Add(h)
    stack.Draw("HIST")
    stack.GetXaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetTitleSize(0.045)
    stack.GetXaxis().SetTitleOffset(1.3)
    stack.GetXaxis().SetTitle(var)
    #stack.GetYaxis().SetTitle("Events/bin")
    stack.GetYaxis().SetTitle("Events/60 GeV") 
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetYaxis().SetTitleSize(0.045)
    stack.SetMaximum(36)
    stack.SetMinimum(0)
        
    # Draw data
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    #data.GetXaxis().SetLimits(0,2000)
    data.Draw("E SAME")
    
    # Draw signals
    for h, color in zip(ZPs, ZPsColors): 
        h.SetLineWidth(1)
        h.SetLineStyle(9)
        h.SetLineColor(color)
        h.Scale(100)
        #h.Scale(tt.Integral("width")/h.Integral("width"))
        #h.GetXaxis().SetLimits(0,2000)
        h.Draw("HIST SAME")
        
    #Define the legend
    legend = ROOT.TLegend(0.60, 0.65, 0.92, 0.92)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.02)
    legend.SetTextAlign(32)
    legend.AddEntry(data, "Data" ,"lep")
    legend.AddEntry(tt, "tt", "f")
    legend.AddEntry(st, "Single top", "f")
    legend.AddEntry(Vjets, "V+jets", "f")
    legend.AddEntry(dib, "Diboson", "f")
    for h_signal, mass in zip(h_signals, masses):
        #legend.AddEntry(h_signal, "Z' " + mass,"f")
        legend.AddEntry(h_signal, "Z' " + mass + " x100","f") 
            
    legend.Draw()
    
    ratio_data_mc = data.Integral("width")/(st.Integral("width") + tt.Integral("width") + dib.Integral("width") + Vjets.Integral("width"))

    # Add ATLAS label
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.045)
    text.DrawLatex(0.21, 0.86, "ATLAS")
    text.SetTextFont(42)
    text.DrawLatex(0.21 + 0.16, 0.86, "Open Data")
    text.SetTextSize(0.04)
    text.DrawLatex(0.21, 0.80, "#sqrt{{s}} = 13 TeV, {:.0f} fb^{{-1}}".format(10))
    text.SetTextSize(0.02)
    text.DrawLatex(0.76, 0.58, "Data/MC: {:.4f}" .format(ratio_data_mc))
    #text.DrawLatex(0.7, 0.53, "Signals normalized to tt yields")
    
    # Save the plot
    c.SaveAs("../run/Stack_"+ var+".pdf")
        

    print(ratio_data_mc)


