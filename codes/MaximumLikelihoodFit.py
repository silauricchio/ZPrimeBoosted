import ROOT
import sys


Zprime_masses =["400", "500", "750", "1000", "1250", "1500", "1750", "2000", "2250", "2500", "2750", "3000"]


# Open files
f = ROOT.TFile.Open("../run/file_histos.root")

#x Axis variable
var_mass = ROOT.RooRealVar("var_mass","mass plot", 0, 3300,"GeV")

# Get histos
data = f.Get("ttbar_mass_data")
tt = f.Get("ttbar_mass_tt")
Vj = f.Get("ttbar_mass_Vjets")
st = f.Get("ttbar_mass_single_top")
dib = f.Get("ttbar_mass_Diboson")

bkg = tt.Clone()
bkg.Add(Vj)
bkg.Add(st)
bkg.Add(dib)

#data hist
data_hist = ROOT.RooDataHist("data_hist","data_hist",ROOT.RooArgList(var_mass),data)
bkg_hist = ROOT.RooDataHist("bkg_hist","bkg_hist",ROOT.RooArgList(var_mass),bkg)

Zprimes = []
Zprimes_hist = []
for mass in Zprime_masses:
        Zprime = f.Get("ttbar_mass_ZPrime"+mass)
        Zprimes.append(Zprime)
        Zprime_hist = ROOT.RooDataHist("Zprime_" + mass, "Zprime_" + mass, ROOT.RooArgList(var_mass), Zprime)
        Zprimes_hist.append(Zprime_hist)

#pdfs
bkg_pdf = ROOT.RooHistPdf("bkg_pdf","bkg_pdf",ROOT.RooArgSet(var_mass),bkg_hist,0)

Zprimes_pdf = []
i = int(0)
for mass in Zprime_masses:
        Zprimes_pdf.append(ROOT.RooHistPdf("Zprime_" + mass + "_pdf", "Zprime_" + mass + "_pdf", ROOT.RooArgSet(var_mass), Zprimes_hist[i],0))
        i = i+1

#bkg yields
nbkg = ROOT.RooRealVar("nbkg","n bkg events",0, bkg.Integral("width")*2)
nbkg.setVal(bkg.Integral())
nbkg.setConstant(ROOT.kTRUE)

#signal yields
nZprimes = []
i = 0
for mass in Zprime_masses:
        nZprime = ROOT.RooRealVar("nZp_" + mass, "Zp " + mass + " events", 0, Zprimes[i].Integral("width")*2)
        nZprime.setVal(Zprimes[i].Integral())
        nZprime.setConstant(ROOT.kTRUE)
        nZprimes.append(nZprime)
        i = i+1

# sig strenght
mu = ROOT.RooRealVar("mu","mu", 0 , 2)

#total bkg normalization factor
k = ROOT.RooRealVar("k", "k", 0.0001, 10)
#k.setVal(1.)
#k.setConstant(ROOT.kTRUE)

#exp signal yields                                                                                                                                                                              
i = 0
exp_signals_yields = []
for mass in Zprime_masses:
        nZp = nZprimes[i]
        exp_signals_yields.append(ROOT.RooFormulaVar("exp_signal_yields_" + mass, "mu*nZp_" + mass , ROOT.RooArgList(mu,nZp)))
        i = i+1

#exp bkgs yields                                                                                                                                                                                   
exp_bkg_yields = ROOT.RooFormulaVar("exp_bkg_yields", "k*nbkg",  ROOT.RooArgList(k,nbkg))



# Construct extended composite models (per signal model)                                                                                                                                                 
# -------------------------------------------------------------------                                                                                                                
# Sum the composite signal and background into an extended pdf                                                                                                                       
# exp_signal_yields+exp_bkg_yields                                                                                                                                                               

i = 0
for mass in Zprime_masses:

        Zprime_pdf = Zprimes_pdf[i]
        exp_signal_yields = exp_signals_yields[i]

        model = ROOT.RooAddPdf("model_" + mass, "Total s(" + mass + ") +b pdf", ROOT.RooArgList(Zprime_pdf, bkg_pdf), ROOT.RooArgList(exp_signal_yields, exp_bkg_yields))

        # construct (binned) negative log likelihood
        nll = model.createNLL(data_hist, ROOT.RooFit.NumCPU(2))

        # Minimize likelihood w.r.t all parameters before making plots
        ROOT.RooMinimizer(nll).migrad()

        # Plot likelihood scan for mu 
        frame1 = mu.frame(ROOT.RooFit.Bins(10),ROOT.RooFit.Range(0, 1.5), ROOT.RooFit.Title("negative LL and profileLL in mu for Z' " + mass + "GeV")) 
        nll.plotOn(frame1,ROOT.RooFit.ShiftToZero(),ROOT.RooFit.Name("NLL_mu"))

        # Plot likelihood scan for k                                                                                                                                                                 
        frame2 = k.frame(ROOT.RooFit.Bins(10), ROOT.RooFit.Range(0.0001, 3), ROOT.RooFit.Title("negative LL and profileLL in k for Z' " + mass + "GeV"))
        nll.plotOn(frame2,ROOT.RooFit.ShiftToZero(), ROOT.RooFit.Name("NLL_k"))
        


        # C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   m u
        # -----------------------------------------------------------------------
        # The profile likelihood estimator on nll for mu will minimize nll w.r.t
        # all floating parameters except mu for each evaluation
        pll_mu = nll.createProfile(ROOT.RooArgSet(mu))
        
        #Plot the profile likelihood in mu
        pll_mu.plotOn(frame1, ROOT.RooFit.Name("NprofileLL_mu"), ROOT.RooFit.LineColor(ROOT.kRed)) 
        
        # Adjust frame maximum for visual clarity
        #frame1.SetMinimum(0) 
        #frame1.SetMaximum(3) 
        
        
        # C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   k                                                                                                                              
        # -----------------------------------------------------------------------                                                                                                                          
        # The profile likelihood estimator on nll for k will minimize nll w.r.t                                                                                                                           
        # all floating parameters except k for each evaluation                                                                                                                                           
        pll_k = nll.createProfile(ROOT.RooArgSet(k))
        #Plot the profile likelihood in k                                                                                                                                                               
        pll_k.plotOn(frame2, ROOT.RooFit.Name("NprofileLL_k"), ROOT.RooFit.LineColor(ROOT.kRed)) 
        # Adjust frame maximum for visual clarity                                                                                                                                                          
        #frame1.SetMinimum(0)                                                                                                                                                                            
        #frame1.SetMaximum(3)  
        
        
        #Make canvas and draw RooPlots
        c = ROOT.TCanvas("profilell","profilell", 800, 400)
        c.Divide(2)
        c.cd(1) 
        #leg1 = ROOT.TLegend()
        #leg1.SetFillStyle(0)
        #leg1.AddEntry("NLL_mu", "NLL", "L")
        #leg1.AddEntry("NprofileLL_mu", "NPLL", "L")
        ROOT.gPad.SetLeftMargin(0.15) 
        frame1.GetYaxis().SetTitleOffset(1.4)
        frame1.Draw() 
        #leg1.Draw()
        
        c.cd(2)
        #leg2 = ROOT.TLegend()
        #leg2.SetFillStyle(0)
        #leg2.AddEntry("NLL_k", "NLL", "L")
        #leg2.AddEntry("NprofileLL_k", "NPLL", "L")
        ROOT.gPad.SetLeftMargin(0.15)
        frame2.GetYaxis().SetTitleOffset(1.4)
        frame2.Draw() 
        #leg2.Draw()
        
        c.SaveAs("../run/ProfileLL_" + mass + ".pdf")
        
        print(mu.getValV())
        
        i = i+1


        #--------------------------------------------------------------------------
        # -------- create the workspace and save pdfs and composite models --------
        #--------------------------------------------------------------------------

        w = ROOT.RooWorkspace("w_" + mass)
        getattr(w,'import')(model)
        getattr(w,'import')(bkg_pdf)
        getattr(w,'import')(Zprime_pdf)
        getattr(w,'import')(data_hist)
        
        #-------------------------------------------------------
        # -------- Set the Model Configuration and Save --------
        #-------------------------------------------------------  
        
        mc = ROOT.RooStats.ModelConfig("ModelConfig_" + mass, w)
        mc.SetPdf(model)
        mc.SetParametersOfInterest(ROOT.RooArgSet(w.var("mu")))
        mc.SetObservables(ROOT.RooArgSet(w.var("var_mass")))
        mc.SetNuisanceParameters(ROOT.RooArgSet(w.var("k"))) #? k nuis par?

        getattr(w,'import')(mc)
        w.writeToFile("../run/Model_" + mass + ".root", True)

        del w
