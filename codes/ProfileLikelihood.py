import ROOT
import sys

def DoFit(stringa):

	# Open files
	f = ROOT.TFile.Open("../run/file_histos.root")

	#x Axis variable
	mass = ROOT.RooRealVar("mass","mass plot", 0, 3300,"GeV")

	# Get histos
	data = f.Get("ttbar_mass_data")
	tt = f.Get("ttbar_mass_tt")
	Vj = f.Get("ttbar_mass_Vjets")
	st = f.Get("ttbar_mass_single_top")
	dib = f.Get("ttbar_mass_Diboson")
	Zp_mass_point = f.Get("ttbar_mass_ZPrime"+stringa)

        bkg = tt.Clone()
        bkg.Add(Vj)
        bkg.Add(st)
        bkg.Add(dib)

	#data hist
	data_hist = ROOT.RooDataHist("data_hist","data_hist",ROOT.RooArgList(mass),data)
	Zp_mass_point_hist = ROOT.RooDataHist("Zp_mass_point_hist","Zp_mass_point_hist",ROOT.RooArgList(mass),Zp_mass_point)
	bkg_hist = ROOT.RooDataHist("bkg_hist","bkg_hist",ROOT.RooArgList(mass),bkg)
        

	#pdfs
	bkg_pdf = ROOT.RooHistPdf("bkg_pdf","bkg_pdf",ROOT.RooArgSet(mass),bkg_hist,0)
	Zp_mass_point_pdf = ROOT.RooHistPdf("Zp_mass_point_pdf","Zp_mass_point_pdf",ROOT.RooArgSet(mass),Zp_mass_point_hist,0)


        #bkg yields
        nbkg = ROOT.RooRealVar("nbkg","n bkg events",0, bkg.Integral("width")*2)
        #nbkg.setVal(bkg.Integral("width"))
        nbkg.setVal(bkg.Integral())
        nbkg.setConstant(ROOT.kTRUE)

        #signal yields
	nZp = ROOT.RooRealVar("nZp","Zp_mass_point events",0,Zp_mass_point.Integral("width")*2)
	#nZp.setVal(Zp_mass_point.Integral("width"))
        nZp.setVal(Zp_mass_point.Integral())
	nZp.setConstant(ROOT.kTRUE)

	# sig strenght
        mu = ROOT.RooRealVar("mu","mu", 0 , 2)

        #total bkg normalization factor
        k = ROOT.RooRealVar("k", "k", 0.0001, 10)
        #k.setVal(1.)
        #k.setConstant(ROOT.kTRUE)

        #exp signal yields                                                                                                                                                                              
        exp_signal_yields = ROOT.RooFormulaVar("exp_signal_yields", "mu*nZp", ROOT.RooArgList(mu,nZp))

        #exp bkgs yields                                                                                                                                                                                   
        exp_bkg_yields = ROOT.RooFormulaVar("exp_bkg_yields", "k*nbkg",  ROOT.RooArgList(k,nbkg))

        # Construct extended composite model                                                                                                                                                                       # -------------------------------------------------------------------                                                                                                                                      # Sum the composite signal and background into an extended pdf                                                                                                                                             # exp_signal_yields+exp_bkg_yields                                                                                                                                                               
        model = ROOT.RooAddPdf("model","Total s+b pdf", ROOT.RooArgList(Zp_mass_point_pdf, bkg_pdf), ROOT.RooArgList(exp_signal_yields, exp_bkg_yields))

        # construct (binned) negative log likelihood
        nll = model.createNLL(data_hist, ROOT.RooFit.NumCPU(2))

        # Minimize likelihood w.r.t all parameters before making plots
        ROOT.RooMinimizer(nll).migrad()
        
        # Plot likelihood scan for mu 
        frame1 = mu.frame(ROOT.RooFit.Bins(10),ROOT.RooFit.Range(0,1), ROOT.RooFit.Title("negative LL and profileLL in mu for Z' " + stringa + "GeV")) 
        nll.plotOn(frame1,ROOT.RooFit.ShiftToZero(),ROOT.RooFit.Name("NLL_mu"))
        
        # Plot likelihood scan for k                                                                                                                                                                 
        frame2 = k.frame(ROOT.RooFit.Bins(10), ROOT.RooFit.Range(0.0001, 3), ROOT.RooFit.Title("negative LL and profileLL in k for Z' " + stringa + "GeV"))
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
        
        c.SaveAs("../run/ProfileLL_" + stringa + ".pdf")


        #calculate Profiled Likelihood Ratio for mu=0
        #idea: pllratio_mu = pll_mu(mu=0)/pll_min
        #then get p-value, then calculate significance
        #then run for all mass points and create a graph

        return mu.getValV()

lista_m =[
#"400","500","750","1000",
"2500",
#"1500"
]
lista_mu =[]
for i in lista_m: lista_mu.append(DoFit(i))
print(lista_m, lista_mu)
