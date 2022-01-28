import ROOT
import matplotlib.pyplot as plt
import numpy as np

masses = ["400","500","750","1000","1250","1500"]
ZPsColors = [ROOT.kBlue, ROOT.kAzure-2, ROOT.kViolet, ROOT.kMagenta, ROOT.kPink-1, ROOT.kYellow-2]

f = ROOT.TFile.Open("file_histos.root")


# Get histos                                                                                                                                                                                        
data = f.Get("ttbar_mass_data")
tt = f.Get("ttbar_mass_tt")
Vjets = f.Get("ttbar_mass_Vjets")
st = f.Get("ttbar_mass_single_top")
dib = f.Get("ttbar_mass_Diboson")
ZPs = []
h_signals = []
for mass in masses:
    h_signal = f.Get("ttbar_mass_ZPrime"+mass)
    h_signals.append(h_signal)
    ZPs.append(h_signal)

# Set styles
ROOT.gROOT.SetStyle("ATLAS")
 
# Create canvas
c = ROOT.TCanvas("c", "", 600, 600)
c.SetTickx(0)
c.SetTicky(0)
#c.SetLogy()
 

stack = ROOT.THStack()
for h, color in zip(
        [dib, Vjets, st, tt],
        [(208, 240, 193), (155, 152, 204), (248, 206, 104), (222, 90, 106)]):
    h.SetLineWidth(1)
    h.SetLineColor(1)
    h.SetFillColor(ROOT.TColor.GetColor(*color))
    h.Scale(0.930433686828)
    stack.Add(h)
stack.Draw("HIST")
stack.GetXaxis().SetLabelSize(0.04)
stack.GetXaxis().SetTitleSize(0.045)
stack.GetXaxis().SetTitleOffset(1.3)
stack.GetXaxis().SetTitle("m^{Z->tt}_{visible} [GeV]")
stack.GetYaxis().SetTitle("Events/60 GeV")
stack.GetYaxis().SetLabelSize(0.04)
stack.GetYaxis().SetTitleSize(0.045)
#stack.SetMaximum(200)
stack.SetMinimum(1)

# Draw data
data.SetMarkerStyle(20)
data.SetMarkerSize(1.2)
data.SetLineWidth(2)
data.SetLineColor(ROOT.kBlack)
data.Draw("E SAME")

# Draw signals
for h, color in zip(ZPs, ZPsColors): 
    h.SetLineWidth(1)
    h.SetLineColor(1)
    h.SetFillColor(color)
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
    legend.AddEntry(h_signal, "Z' " + mass ,"f")

legend.Draw()

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
 
# Save the plot
c.SaveAs("Stack_Zprime_norm.pdf")


print(data.Integral("width")/(st.Integral("width") + tt.Integral("width") + dib.Integral("width") + Vjets.Integral("width")))


