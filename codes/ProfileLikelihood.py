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

	#data hist
	data_hist = ROOT.RooDataHist("data_hist","data_hist",ROOT.RooArgList(mass),data)
	tt_hist = ROOT.RooDataHist("tt_hist","tt_hist",ROOT.RooArgList(mass),tt)
	Vj_hist = ROOT.RooDataHist("Vj_hist","Vj_hist",ROOT.RooArgList(mass),Vj)
	st_hist = ROOT.RooDataHist("st_hist","st_hist",ROOT.RooArgList(mass),st)
	dib_hist = ROOT.RooDataHist("dib_hist","dib_hist",ROOT.RooArgList(mass),dib)
	Zp_mass_point_hist = ROOT.RooDataHist("Zp_mass_point_hist","Zp_mass_point_hist",ROOT.RooArgList(mass),Zp_mass_point)
	
	#pdfs
	tt_pdf = ROOT.RooHistPdf("tt_pdf","tt_pdf",ROOT.RooArgSet(mass),tt_hist,0)
	Vj_pdf = ROOT.RooHistPdf("Vj_pdf","Vj_pdf",ROOT.RooArgSet(mass),Vj_hist,0)
	st_pdf = ROOT.RooHistPdf("st_pdf","st_pdf",ROOT.RooArgSet(mass),st_hist,0)
	dib_pdf = ROOT.RooHistPdf("dib_pdf","dib_pdf",ROOT.RooArgSet(mass),dib_hist,0)
	Zp_mass_point_pdf = ROOT.RooHistPdf("Zp_mass_point_pdf","Zp_mass_point_pdf",ROOT.RooArgSet(mass),Zp_mass_point_hist,0)

	# tot numb events observed
	N_obs = ROOT.RooRealVar("n","n events",0,data.Integral("width")*2)
	N_obs.setVal(data.Integral("width"))
	N_obs.setConstant(ROOT.kTRUE)

	# sig and bkg events observed
	ntt = ROOT.RooRealVar("ntt","tt events",0,tt.Integral("width")*2)
	ntt.setVal(tt.Integral("width"))
	ntt.setConstant(ROOT.kTRUE)

	nVj = ROOT.RooRealVar("nVj","Vj events",0,Vj.Integral("width")*2)
	nVj.setVal(Vj.Integral("width"))
	nVj.setConstant(ROOT.kTRUE)

	nst = ROOT.RooRealVar("nst","st events",0,st.Integral("width")*2)
	nst.setVal(st.Integral("width"))
	nst.setConstant(ROOT.kTRUE)

	ndib = ROOT.RooRealVar("ndib","dib events",0,dib.Integral("width")*2)
	ndib.setVal(dib.Integral("width"))
	ndib.setConstant(ROOT.kTRUE)

        print(ndib, nst, nVj, ntt)
        
        #signal yields
	nZp = ROOT.RooRealVar("nZp","Zp_mass_point events",0,Zp_mass_point.Integral("width")*2)
	nZp.setVal(Zp_mass_point.Integral("width"))
	nZp.setConstant(ROOT.kTRUE)

	# sig strenght
        mu = ROOT.RooRealVar("mu","mu", 0 , 10)

        #total bkg normalization factor
        k = ROOT.RooRealVar("k", "k", 0, 100)
        #k.setVal(1.)
        #k.setConstant(ROOT.kTRUE)

        #exp signal yields                                                                                                                                                                              
        exp_signal_yields = ROOT.RooFormulaVar("exp_signal_yields", "mu*nZp", ROOT.RooArgList(mu,nZp))

        #exp bkgs yields                                                                                                                                                                                   
        exp_bkg_dib = ROOT.RooFormulaVar("exp_bkg_dib", "k*ndib",  ROOT.RooArgList(k,ndib))
        exp_bkg_st = ROOT.RooFormulaVar("exp_bkg_st", "k*nst",  ROOT.RooArgList(k,nst))
        exp_bkg_Vj = ROOT.RooFormulaVar("exp_bkg_Vj", "k*nVj",  ROOT.RooArgList(k,nVj))
        exp_bkg_tt = ROOT.RooFormulaVar("exp_bkg_tt", "k*ntt",  ROOT.RooArgList(k,ntt))

        # Construct extended composite model                                                                                                                                                                       # -------------------------------------------------------------------                                                                                                                                      # Sum the composite signal and background into an extended pdf                                                                                                                                             # exp_signal_yields+exp_bkg_yields                                                                                                                                                               
        model = ROOT.RooAddPdf("model","Total s+b pdf", ROOT.RooArgList(Zp_mass_point_pdf, dib_pdf, st_pdf, Vj_pdf, tt_pdf), ROOT.RooArgList(exp_signal_yields, exp_bkg_dib, exp_bkg_st, exp_bkg_Vj, exp_bkg_tt))

        # construct (binned) negative log likelihood
        nll = model.createNLL(data_hist, ROOT.RooFit.NumCPU(2))

        # Minimize likelihood w.r.t all parameters before making plots
        ROOT.RooMinimizer(nll).migrad()
        
        
        # Plot likelihood scan for mu 
        frame1 = mu.frame(ROOT.RooFit.Bins(10),ROOT.RooFit.Range(0,1), ROOT.RooFit.Title("negative LL and profileLL in mu")) 
        nll.plotOn(frame1,ROOT.RooFit.ShiftToZero(),ROOT.RooFit.Name("NLL_mu"))
        
        # Plot likelihood scan for k                                                                                                                                                                 
        frame2 = k.frame(ROOT.RooFit.Bins(10), ROOT.RooFit.Range(0,0.001), ROOT.RooFit.Title("negative LL and profileLL in k"))
        nll.plotOn(frame2,ROOT.RooFit.ShiftToZero(), ROOT.RooFit.Name("NLL_k"))
        
        '''
        
        # C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   m u
        # -----------------------------------------------------------------------
        # The profile likelihood estimator on nll for mu will minimize nll w.r.t
        # all floating parameters except mu for each evaluation
        pll_mu = nll.createProfile(mu)
        #Plot the profile likelihood in mu
        pll_mu.plotOn(frame1, ROOT.RooFit.Name("NprofileLL_mu"), ROOT.RooFit.LineColor(kRed)) 
        # Adjust frame maximum for visual clarity
        #frame1.SetMinimum(0) 
        #frame1.SetMaximum(3) 

        
        # C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   k                                                                                                                              
        # -----------------------------------------------------------------------                                                                                                                          
        # The profile likelihood estimator on nll for k will minimize nll w.r.t                                                                                                                           
        # all floating parameters except k for each evaluation                                                                                                                                           
        pll_k = nll.createProfile({k})
        #Plot the profile likelihood in k                                                                                                                                                               
        pll_k.plotOn(frame2, ROOT.RooFit.Name("NprofileLL_k"), ROOT.RooFit.LineColor(kRed)) 
        # Adjust frame maximum for visual clarity                                                                                                                                                          
        #frame1.SetMinimum(0)                                                                                                                                                                            
        #frame1.SetMaximum(3)  
        '''


        #Make canvas and draw RooPlots
        c = ROOT.TCanvas("profilell","profilell", 800, 400)
        c.Divide(2)
        c.cd(1) 
        leg1 = ROOT.TLegend()
        leg1.AddEntry("NLL_mu", "NLL")
        #leg1.AddEntry("NprofileLL_mu", "NPLL")
        ROOT.gPad.SetLeftMargin(0.15) 
        frame1.GetYaxis().SetTitleOffset(1.4)
        frame1.Draw() 
        leg1.Draw()
        
        c.cd(2)
        leg2 = ROOT.TLegend()
        leg2.AddEntry("NLL_k", "NLL")
        #leg2.AddEntry("NprofileLL_k", "NPLL")
        ROOT.gPad.SetLeftMargin(0.15)
        frame2.GetYaxis().SetTitleOffset(1.4)
        frame2.Draw() 
        leg2.Draw()
        
        c.SaveAs("../run/ProfileLL.pdf")

        return mu.getValV()

lista_m =[
#"400","500","750","1000",
"2500",
#"1500"
]
lista_mu =[]
for i in lista_m: lista_mu.append(DoFit(i))
print(lista_m, lista_mu)
