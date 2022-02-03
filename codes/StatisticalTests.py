import ROOT
import sys


poiValueForBackground = 0

Zprime_masses = ["400", "500", "750", "1000", "1250", "1500", "1750", "2000", "2250", "2500", "2750", "3000"]
#Zprime_masses = [
#"1000"
#, "1250", "1500"
#]


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

    # ---------------------------------------------------
    # ----- Configure a ProfileLikelihoodCalculator -----
    # ---------------------------------------------------

    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, model, ROOT.RooArgSet(POI)) #this class uses the Wilk's theorem and calculate p value with asymptotic formulae

    '''
    plc.SetConfidenceLevel(0.95)
    interval = plc.GetInterval()
    mu_best = interval.GetBestFitParameters().first()

    lowerLimit = interval.LowerLimit(POI)
    upperLimit = interval.UpperLimit(POI)
    print("MU = " + str(mu_best))
    

    #Make a plot
    c = ROOT.TCanvas()
    plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
    plot.SetLineColor(ROOT.kWhite)
    #plot.SetNPoints(50)  #do not use too many points, it could become very slow for some models
    plot.Draw()
    #plot.Draw("tf1") #use option TF1 if too slow (plot.Draw("tf1"))
    c.SaveAs("../run/PLRScan_mu_Zp" + mass + ".pdf")
    '''

    # Get the significance using the GetHypoTest function:
    # create a copy of the POI parameters to set the values to zero

    nullparams = mc.GetParametersOfInterest().Clone()
    nullparams.first().setVal(0)
    plc.SetNullParameters(nullparams)

    del nullparams
    
    #Return the hypothesis test result obtained from the likelihood ratio of the maximum likelihood value with the null parameters fixed to their values, with respect keeping all parameters floating
    testresult = plc.GetHypoTest()
    significance = testresult.Significance()    
    p_value = testresult.NullPValue() 
    

    print("Observed significance with asym formulae")
    print("For Z'" + mass + " p value --> " + str(p_value) + " and Z --> " + str(significance))

    del plc

    
    # ---------------------------------------
    # ----- repeat the test with Toys MC ----
    # ---------------------------------------

    NTOYS = 1000

    # Get RooStats ModelConfig from workspace and save it as Signal+Background model
    sbModel = mc.Clone()
    poi = sbModel.GetParametersOfInterest().first()
    poi.setVal(1.)
    sbModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    # Clone S+B model, set POI to zero and set as B-Only model
    bModel = sbModel.Clone()
    bModel.SetName('bmodel')
    poi.setVal(0.)
    bModel.SetSnapshot(ROOT.RooArgSet(poi))
    
    hc = ROOT.RooStats.FrequentistCalculator(data, bModel, sbModel)
    hc.SetToys(int(NTOYS), int(NTOYS/2.))
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

        
    # HypoTestInverter #per l'estrazione dei limiti
    calc = ROOT.RooStats.HypoTestInverter(hc)
    calc.SetConfidenceLevel(0.95)
    calc.UseCLs(True)
    calc.SetVerbose(False)
    calc.SetAutoScan()
    #calc.SetFixedScan(6, 0., 10., False)
    hypotestresult = calc.GetInterval()
    
    # Plot HypoTestInverter result
    c1 = ROOT.TCanvas()
    plot = ROOT.RooStats.HypoTestInverterPlot('hypotest_inverter', 'hypotest inverter', hypotestresult)
    plot.MakePlot()
    plot.MakeExpectedPlot()
    plot.Draw('CLB2CL')
    c1.SaveAs("../run/boh_" + mass + ".pdf")

    # Print information
    hc.GetFitInfo().Print('v')


    del hc

    
