import ROOT
import sys

def rooPart1(stringa):

	# Open files
	f = ROOT.TFile.Open("../run/file_histos.root")

	#x Axis variable
	mass = ROOT.RooRealVar("mass","mass plot",400,1600,"GeV")

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

	# sig and bkcg evts oberved
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

	nZp_mass_point = ROOT.RooRealVar("nZp_mass_point","Zp_mass_point events",0,Zp_mass_point.Integral("width")*2)
	nZp_mass_point.setVal(Zp_mass_point.Integral("width"))
	nZp_mass_point.setConstant(ROOT.kTRUE)


	# sig strenght
	#mu = ROOT.RooRealVar("mu","mu",1,0,2)
        mu = ROOT.RooRealVar("mu","mu", 0 , 100)

        #total bkg normalization factor
        k = ROOT.RooRealVar("k", "k", 0, 10.)
        k.setVal(0.930433686828)
        #k.setVal(1.)
        k.setConstant(ROOT.kTRUE)

	#signal yields *mu
	signal = ROOT.RooFormulaVar("signal", "mu*nZp_mass_point", ROOT.RooArgList(mu,nZp_mass_point))

        #bkgs yields
        bkg_dib = ROOT.RooFormulaVar("bkg_dib", "k*ndib",  ROOT.RooArgList(k,ndib))
        bkg_st = ROOT.RooFormulaVar("bkg_st", "k*nst",  ROOT.RooArgList(k,nst))
        bkg_Vj = ROOT.RooFormulaVar("bkg_Vj", "k*nVj",  ROOT.RooArgList(k,nVj))
        bkg_tt = ROOT.RooFormulaVar("bkg\_tt", "k*ntt",  ROOT.RooArgList(k,ntt))

	#expected mu*s+b(*k)
	total_events = ROOT.RooFormulaVar("total_events", "mu*nZp_mass_point + k*(ndib+nst+nVj+ntt)", ROOT.RooArgList(mu,nZp_mass_point, k, ndib,nst,nVj,ntt))
        #TotEvents_pdf
	TotEvents_pdf = ROOT.RooPoisson("pois","Total number of events",N_obs,total_events)

        #total pdf 
	sum_pdf = ROOT.RooAddPdf("sumPDF","Total s+b PDF",ROOT.RooArgList(Zp_mass_point_pdf,dib_pdf,st_pdf,Vj_pdf,tt_pdf),ROOT.RooArgList(signal, bkg_dib, bkg_st, bkg_Vj, bkg_tt))
	prod_sum_pdf = ROOT.RooProdPdf("extPDF","Extended Maximum Likelihood",ROOT.RooArgList(sum_pdf, TotEvents_pdf))

	#fit
	prod_sum_pdf.fitTo(data_hist)# ROOT.RooFit.Extended(1))

	
	#c1=ROOT.TCanvas()

	#frame = mass.frame(ROOT.RooFit.Title("mass distribution"))
	#tt_arg = ROOT.RooArgSet(tt_pdf)
        #Vj_arg =  ROOT.RooArgSet(Vj_pdf)
        #st_arg = ROOT.RooArgSet(st_pdf)
        #dib_arg = ROOT.RooArgSet(dib_pdf)
	#Zp_arg = ROOT.RooArgSet(Zp_mass_point_pdf)
	#tt_pdf.plotOn(frame)
	#Vj_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed))
	#st_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray))
	#dib_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen))
	#Zp_mass_point_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kMagenta))
	#prod_sum_pdf.plotOn(frame,ROOT.RooFit.Components(tt_arg),ROOT.RooFit.LineColor(ROOT.kBlack))
	#prod_sum_pdf.plotOn(frame,ROOT.RooFit.Components(Zp_arg),ROOT.RooFit.LineColor(ROOT.kBlue))
	##frame.Draw()
	
	#data.Draw("Pesame")
	#c1.SaveAs("pdfs.png")
	
	return mu.getValV()

lista_m =["400","500","750","1000","1250","1500"]
lista_mu =[]
for i in lista_m: lista_mu.append(rooPart1(i))
print(lista_m, lista_mu)
