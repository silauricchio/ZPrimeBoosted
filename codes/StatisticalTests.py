import ROOT
import sys

run_Asimov = True
run_Toys = False

poiValueForBackground = 0

Zprime_masses = ["400", "500", "750", "1000", "1250", "1500", "1750", "2000", "2250", "2500", "2750", "3000"]
#Zprime_masses = ["3000"]


AsymLimits = {}
ExpAsymLimits = {}
ExpAsymLimits_m1s = {}
ExpAsymLimits_p1s = {}
ExpAsymLimits_m2s = {}
ExpAsymLimits_p2s = {}

ToysLimits = {}
ExpToysLimits = {}
ExpToysLimits_m1s = {}
ExpToysLimits_p1s = {}
ExpToysLimits_m2s = {}
ExpToysLimits_p2s = {}


for mass in Zprime_masses:

    # Open models files                                                                                                                                                                                   
    f = ROOT.TFile.Open("../run/Model_" + mass + ".root")
    
    #retrieve the workspace
    
    w = f.Get("w_" + mass)
    w.Print() 

    data = w.data("data_hist") 
    mc = w.obj("ModelConfig_" + mass)
    model = w.pdf("model_" + mass)

    # make sure ingredients are found
    if not data or not mc or not model:   
        print("data or ROOT.RooStats.ModelConfig or model was not found")
        exit

    POI = mc.GetParametersOfInterest().first()
    NP = mc.GetNuisanceParameters().first()

    
    # ---------------------------------------------------
    # ----- Configure a ProfileLikelihoodCalculator -----
    # ---------------------------------------------------
    

    sbModel = mc.Clone()                                                                                                                                                                                  
    sbModel.SetName('sbmodel')   
    poi = sbModel.GetParametersOfInterest().first()                                                                                                                                                        
    poi.setVal(1.)                                                                                                                                                                                         
    sbModel.SetSnapshot(ROOT.RooArgSet(poi))                                                                                                                                                               
    
    # Clone S+B model, set POI to zero and set as B-Only model                                                                                                                                             
    bModel = sbModel.Clone()                                                                                                                                                                               
    bModel.SetName('bmodel')                                                                                                                                                                               
    poi_null = poi.Clone()                                                                                                                                                                                 
    poi_null.setVal(0.)               
    bModel.SetSnapshot(ROOT.RooArgSet(poi_null))         
    

    if run_Asimov:

        print('')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('ASYMPTOTIC CALCULATION')
        print('')
        
        AsymCalc = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel, False)  
        AsymCalc.SetOneSided(True)
        
        #Return the hypothesis test result obtained from the likelihood ratio of the maximum likelihood value with the null parameters fixed to their values, with respect keeping all parameters floating
        testresult = AsymCalc.GetHypoTest()
        significance = testresult.Significance()    
        p_value = testresult.NullPValue() 
        
        print("Observed significance with asym formulae")
        print("For Z'" + mass + " p value --> " + str(p_value) + " and Z --> " + str(significance))
        
        
        fit_params = AsymCalc.GetBestFitParams() 
        print(fit_params)
        mu_best = AsymCalc.GetMuHat() 
        
        
        # HypoTestInverter #per l'estrazione dei limiti                                                                                                                                                    

        calc = ROOT.RooStats.HypoTestInverter(AsymCalc)                                                                                                                                                    
        calc.SetConfidenceLevel(0.95)                                                                                                                                                                      
        calc.UseCLs(True)                                                                                                                                                                                  
        calc.SetVerbose(False)                                                                                                                                                                             
        calc.SetAutoScan()                                                                                                                                                                                 
        hypotestresult = calc.GetInterval()   
        
        ObsUpperLimit = hypotestresult.UpperLimit()
        ExpUpperLimit = hypotestresult.GetExpectedUpperLimit(0)
        ExpUpperLimit_m1s = hypotestresult.GetExpectedUpperLimit(-1)
        ExpUpperLimit_p1s = hypotestresult.GetExpectedUpperLimit(1)
        ExpUpperLimit_m2s = hypotestresult.GetExpectedUpperLimit(-2)
        ExpUpperLimit_p2s = hypotestresult.GetExpectedUpperLimit(2)
        
        
        print(mass)
        print("MU up obs = " + str(ObsUpperLimit))
        print("MU up exp = " + str(ExpUpperLimit))
        print('') 
        
        AsymLimits['mass'] = ObsUpperLimit
        ExpAsymLimits['mass'] = ExpUpperLimit
        ExpAsymLimits_m1s['mass'] = ExpUpperLimit_m1s
        ExpAsymLimits_p1s['mass'] = ExpUpperLimit_p1s
        ExpAsymLimits_m2s['mass'] = ExpUpperLimit_m2s
        ExpAsymLimits_p2s['mass'] = ExpUpperLimit_p2s
        
        del AsymCalc
    
    
    # ---------------------------------------
    # ----- repeat the test with Toys MC ----
    # ---------------------------------------
    
    if run_Toys:
        print('')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('TOYS CALCULATION')
        print('')
        
        
        NTOYS = 10000
        
        hc = ROOT.RooStats.FrequentistCalculator(data, bModel, sbModel) 
        hc.SetToys(int(NTOYS), int(NTOYS))
        hc.StoreFitInfo(True)
        hc.UseSameAltToys()
        
        # Test statistics: profile likelihood
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(bModel.GetPdf())
        profll.EnableDetailedOutput()
        profll.SetLOffset(True)
        profll.SetMinimizer('Minuit2')
        profll.SetOneSided(True)
        profll.SetPrintLevel(0)
        profll.SetStrategy(2)
        profll.SetAlwaysReuseNLL(True)
        profll.SetReuseNLL(True)
        
        toymcs = hc.GetTestStatSampler()
        toymcs.SetTestStatistic(profll)
        toymcs_testresult = hc.GetHypoTest()
        toymcs_significance = toymcs_testresult.Significance()
        toymcs_p_value = toymcs_testresult.NullPValue()
        
        
        print("Observed significance with Toys MC")
        print("For Z'" + mass + " p value --> " + str(toymcs_p_value) + " and Z --> " + str(toymcs_significance))
        exit()
        
        # HypoTestInverter #per l'estrazione dei limiti
        toymc_calc = ROOT.RooStats.HypoTestInverter(hc)
        toymc_calc.SetConfidenceLevel(0.95)
        toymc_calc.UseCLs(True)
        toymc_calc.SetVerbose(False)
        toymc_calc.SetAutoScan()
        toymc_hypotestresult = toymc_calc.GetInterval()
        
        toymc_ObsUpperLimit = toymc_hypotestresult.UpperLimit()
        toymc_ExpUpperLimit = toymc_hypotestresult.GetExpectedUpperLimit(0)
        toymc_ExpUpperLimit_m1s = toymc_hypotestresult.GetExpectedUpperLimit(-1)
        toymc_ExpUpperLimit_p1s = toymc_hypotestresult.GetExpectedUpperLimit(1)
        toymc_ExpUpperLimit_m2s = toymc_hypotestresult.GetExpectedUpperLimit(-2)
        toymc_ExpUpperLimit_p2s = toymc_hypotestresult.GetExpectedUpperLimit(2)
        
        print(mass)
        print("MU up obs = " + str(toymc_ObsUpperLimit))
        print("MU up exp = " + str(toymc_ExpUpperLimit))
        print('')
        
        ToysLimits['mass'] = toymc_ObsUpperLimit
        ExpToysLimits['mass'] = toymc_ExpUpperLimit
        ExpToysLimits_m1s['mass'] = toymc_ExpUpperLimit_m1s
        ExpToysLimits_p1s['mass'] = toymc_ExpUpperLimit_p1s
        ExpToysLimits_m2s['mass'] = toymc_ExpUpperLimit_m2s
        ExpToysLimits_p2s['mass'] = toymc_ExpUpperLimit_p2s
        
        del hc
    


#p_value_graph_asym = ROOT.TGraph("p value", "p value", masses, p_values)
#significance_graph_asym = 
#p_value_graph_toys =
#significance_graph_toys =

#upper_limits_asym= brazilian plot
#upper_limits_toys=

'''
# Plot HypoTestInverter result
c1 = ROOT.TCanvas()
plot = ROOT.RooStats.HypoTestInverterPlot('hypotest_inverter', 'hypotest inverter', hypotestresult)
plot.MakePlot()
plot.MakeExpectedPlot()
plot.Draw('CLB2CL')
c1.SaveAs("../run/boh_" + mass + ".pdf")
'''

