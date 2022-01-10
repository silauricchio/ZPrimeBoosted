import ROOT
import sys

def rooPart1():

	mass = ROOT.RooRealVar("mass","mass plot",0,1000,"GeV")
	#nev da 0 500 hist Data ha circa 200 entries
	#real value n_OBS
	nev = ROOT.RooRealVar("n","n events",0,500)

	nev.setVal(198)
	# valori plausibili da h->GetIntegral()
	#s
	nsignal = 3
	#b
	nbackground = 190
	mu = 1


	signal = ROOT.RooRealVar("signal","signal events",nsignal,0,20)
	background = ROOT.RooRealVar("background","signal events",nbackground,0,500)
	mu = ROOT.RooRealVar("mu","mu",mu,0,2) 	
	#expected mu*s+b
	total_events = ROOT.RooFormulaVar("total_events", "mu*signal+background", ROOT.RooArgList(mu,signal,background))
	ExpEvents = ROOT.RooPoisson("pois","Total number of events",nev,total_events)
	#prod_sum_pdf = ROOT.RooProdPdf("extPDF","Total s+b PDF with nevents",ROOT.RooArgList(ExpEvents))

	

rooPart1()