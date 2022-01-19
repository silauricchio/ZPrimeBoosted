import ROOT
import sys

def rooPart1(stringa):

	# Open files
	f = ROOT.TFile.Open("file_histos.root")

	#x Axis variable
	mass = ROOT.RooRealVar("mass","mass plot",400,1600,"GeV")

	# Get histos
	data = f.Get("ttbar_mass_data")
	tt = f.Get("ttbar_mass_tt")
	Vj = f.Get("ttbar_mass_Vjets")
	st = f.Get("ttbar_mass_single_top")
	dib = f.Get("ttbar_mass_Diboson")
	Zp_400 = f.Get("ttbar_mass_ZPrime"+stringa)

	#data hist
	data_hist = ROOT.RooDataHist("data_hist","data_hist",ROOT.RooArgList(mass),data)
	tt_hist = ROOT.RooDataHist("tt_hist","tt_hist",ROOT.RooArgList(mass),tt)
	Vj_hist = ROOT.RooDataHist("Vj_hist","Vj_hist",ROOT.RooArgList(mass),Vj)
	st_hist = ROOT.RooDataHist("st_hist","st_hist",ROOT.RooArgList(mass),st)
	dib_hist = ROOT.RooDataHist("dib_hist","dib_hist",ROOT.RooArgList(mass),dib)
	Zp_400_hist = ROOT.RooDataHist("Zp_400_hist","Zp_400_hist",ROOT.RooArgList(mass),Zp_400)
	
	#pdfs
	tt_pdf = ROOT.RooHistPdf("tt_pdf","tt_pdf",ROOT.RooArgSet(mass),tt_hist,0)
	Vj_pdf = ROOT.RooHistPdf("Vj_pdf","Vj_pdf",ROOT.RooArgSet(mass),Vj_hist,0)
	st_pdf = ROOT.RooHistPdf("st_pdf","st_pdf",ROOT.RooArgSet(mass),st_hist,0)
	dib_pdf = ROOT.RooHistPdf("dib_pdf","dib_pdf",ROOT.RooArgSet(mass),dib_hist,0)
	Zp_400_pdf = ROOT.RooHistPdf("Zp_400_pdf","Zp_400_pdf",ROOT.RooArgSet(mass),Zp_400_hist,0)



	# tot numb events observed
	nev = ROOT.RooRealVar("n","n events",0,500)
	nev.setVal(data.Integral())
	nev.setConstant(ROOT.kTRUE)

	# sig and bkcg evts oberved
	ntt = ROOT.RooRealVar("ntt","tt events",0,tt.Integral()*2)
	ntt.setVal(tt.Integral())
	#ntt.setConstant(ROOT.kTRUE)

	nVj = ROOT.RooRealVar("nVj","Vj events",0,Vj.Integral()*2)
	nVj.setVal(Vj.Integral())
	nVj.setConstant(ROOT.kTRUE)

	nst = ROOT.RooRealVar("nst","st events",0,st.Integral()*2)
	nst.setVal(st.Integral())
	nst.setConstant(ROOT.kTRUE)

	ndib = ROOT.RooRealVar("ndib","dib events",0,dib.Integral()*2)
	ndib.setVal(dib.Integral())
	ndib.setConstant(ROOT.kTRUE)

	nZp_400 = ROOT.RooRealVar("nZp_400","Zp_400 events",0,Zp_400.Integral()*2)
	nZp_400.setVal(Zp_400.Integral())
	nZp_400.setConstant(ROOT.kTRUE)


	# sig strenght
	mu = ROOT.RooRealVar("mu","mu",1,0,2)

	#signal
	signal = ROOT.RooFormulaVar("signal", "mu*nZp_400", ROOT.RooArgList(mu,nZp_400))

	#expected mu*s+b
	total_events = ROOT.RooFormulaVar("total_events", "mu*nZp_400+ndib+nst+nVj+ntt", ROOT.RooArgList(mu,nZp_400,ndib,nst,nVj,ntt))
	ExpEvents = ROOT.RooPoisson("pois","Total number of events",nev,total_events)


	sum_pdf = ROOT.RooAddPdf("sumPDF","Total s+b PDF",ROOT.RooArgList(Zp_400_pdf,dib_pdf,st_pdf,Vj_pdf,tt_pdf),ROOT.RooArgList(signal,ndib,nst,nVj,ntt))
	prod_sum_pdf = ROOT.RooProdPdf("extPDF","Total s+b PDF with nevents",ROOT.RooArgList(sum_pdf,ExpEvents))
	#fit
	prod_sum_pdf.fitTo(data_hist)#, ROOT.RooFit.Extended(1))

	"""
	c1=ROOT.TCanvas()

	frame = mass.frame(ROOT.RooFit.Title("mass distribution"))
	tt_arg = ROOT.RooArgSet(tt_pdf)
	Zp_arg = ROOT.RooArgSet(Zp_400_pdf)
	#tt_pdf.plotOn(frame)
	#Vj_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed))
	#st_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray))
	#dib_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen))
	#Zp_400_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kMagenta))
	prod_sum_pdf.plotOn(frame,ROOT.RooFit.Components(tt_arg),ROOT.RooFit.LineColor(ROOT.kBlack))
	prod_sum_pdf.plotOn(frame,ROOT.RooFit.Components(Zp_arg),ROOT.RooFit.LineColor(ROOT.kBlue))
	frame.Draw()
	
	data.Draw("Pesame")
	c1.SaveAs("pdfs.png")
	"""
	return mu.getValV()

lista_m =["400","500","750","1000","1250","1500"]
lista_mu =[]
for i in lista_m: lista_mu.append(rooPart1(i))
print(lista_mu)