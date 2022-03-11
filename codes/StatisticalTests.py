import ROOT
import sys
import matplotlib.pyplot as plt
import numpy as np
import collections
ROOT.gROOT.SetBatch(ROOT.kTRUE)

run_Asimov = False
run_Toys = True
NTOYS = 7000

poiValueForBackground = 0

Zprime_masses = ["400", "500", "750", "1000", "1250", "1500", "1750", "2000", "2250", "2500", "2750", "3000"] #GeV
#Zprime_masses = ["2000"]

Zprime_xsections = [8.9856996, 8.7384996, 3.1201000, 1.1260999, 0.4598099, 0.2068500, 0.1001600, 0.0513460, 0.0274810, 0.015226, 0.0086883, 0.0050842] #fb
#Zprime_xsections = [0.0513460]

Asym_pvalue = {}
Asym_sig = {}

Asym_exp_pvalue = {}
Asym_exp_pvalue_1s = {}
Asym_exp_pvalue_2s = {}
Asym_exp_pvalue_m1s = {}
Asym_exp_pvalue_m2s = {}

Asym_exp_sig = {}
Asym_exp_sig_1s = {}
Asym_exp_sig_2s = {}
Asym_exp_sig_m1s = {}
Asym_exp_sig_m2s = {}

Toys_pvalue = {}
Toys_sig = {}

Toys_exp_pvalue = {}
Toys_exp_pvalue_1s = {}
Toys_exp_pvalue_2s = {}
Toys_exp_pvalue_m1s = {}
Toys_exp_pvalue_m2s = {}

Toys_exp_sig = {}
Toys_exp_sig_1s = {}
Toys_exp_sig_2s = {}
Toys_exp_sig_m1s = {}
Toys_exp_sig_m2s = {}


AsymLimits = {}
ExpAsymLimits = {}
ExpAsymLimits_1s = {}
ExpAsymLimits_2s = {}
ExpAsymLimits_m1s = {}
ExpAsymLimits_m2s = {}

ToysLimits = {}
ExpToysLimits = {}
ExpToysLimits_1s = {}
ExpToysLimits_2s = {}
ExpToysLimits_m1s = {}
ExpToysLimits_m2s = {}


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
        
        AsymCalc = ROOT.RooStats.AsymptoticCalculator(data, sbModel, bModel, False) #False if for Asimov data with fitted k
        AsymCalc.SetOneSidedDiscovery(True) 

        #Return the hypothesis test result obtained from the likelihood ratio of the maximum likelihood value with mu=0 (null hypothesis), with respect keeping all parameters floating
        testresult = AsymCalc.GetHypoTest()
        significance = testresult.Significance()    
        p_value = testresult.NullPValue() 
        
        '''
        can = ROOT.TCanvas()
        plot = ROOT.RooStats.HypoTestPlot(testresult)
        plot.Draw()
        can.SaveAs('../run/TestStat_Distr_Asym_' + mass +'.pdf')
        '''
        print("Observed null significance with asym formulae")
        print("For Z' " + mass + " p value --> " + str(p_value) + " and Z --> " + str(significance))
        print("Observed null significance (in s+b assumption) with asym formulae")
        print("For Z' " + mass + " p value --> " + str(testresult.AlternatePValue()) + " and Z --> " + str(ROOT.RooStats.PValueToSignificance(testresult.AlternatePValue())))
  
        Exp_pvalue = AsymCalc.GetExpectedPValues(p_value, testresult.AlternatePValue(), 0,False)
        Exp_pvalue_1s = AsymCalc.GetExpectedPValues(p_value, testresult.AlternatePValue(), 1, False)
        Exp_pvalue_2s =AsymCalc.GetExpectedPValues(p_value, testresult.AlternatePValue(), 2, False)
        Exp_pvalue_m1s = AsymCalc.GetExpectedPValues(p_value, testresult.AlternatePValue(), -1, False)
        Exp_pvalue_m2s =AsymCalc.GetExpectedPValues(p_value, testresult.AlternatePValue(), -2, False)

        Exp_significance = ROOT.RooStats.PValueToSignificance(Exp_pvalue)
        Exp_significance_1s = ROOT.RooStats.PValueToSignificance(Exp_pvalue_1s)
        Exp_significance_2s = ROOT.RooStats.PValueToSignificance(Exp_pvalue_2s)
        Exp_significance_m1s = ROOT.RooStats.PValueToSignificance(Exp_pvalue_m1s)
        Exp_significance_m2s = ROOT.RooStats.PValueToSignificance(Exp_pvalue_m2s)
        
        Asym_pvalue[int(mass)] = p_value
        Asym_sig[int(mass)] = significance

        Asym_exp_pvalue[int(mass)] = Exp_pvalue
        Asym_exp_pvalue_1s[int(mass)] = Exp_pvalue_1s
        Asym_exp_pvalue_2s[int(mass)] = Exp_pvalue_2s
        Asym_exp_pvalue_m1s[int(mass)] = Exp_pvalue_m1s
        Asym_exp_pvalue_m2s[int(mass)] = Exp_pvalue_m2s


        Asym_exp_sig[int(mass)] = Exp_significance
        Asym_exp_sig_1s[int(mass)] = Exp_significance_1s
        Asym_exp_sig_2s[int(mass)] = Exp_significance_2s
        Asym_exp_sig_m1s[int(mass)] = Exp_significance_m1s
        Asym_exp_sig_m2s[int(mass)] = Exp_significance_m2s

        fit_params = AsymCalc.GetBestFitParams() 
        print(fit_params)
        mu_best = AsymCalc.GetMuHat() 
        
        del AsymCalc
        
        # HypoTestInverter for limits extraction                                                                                                                                                           
        AsymCalc = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel, False)
        AsymCalc.SetOneSided(True)
    
        calc = ROOT.RooStats.HypoTestInverter(AsymCalc)                                                                                                                                                    
        calc.SetConfidenceLevel(0.95)                                                                                                                                                                      
        calc.UseCLs(True)                                                                                                                                                                                  
        calc.SetVerbose(True)                                                                                                                                                                             
        calc.SetAutoScan()                                                                                                                                                                                 
        hypotestresult = calc.GetInterval()   
        
        ObsUpperLimit = hypotestresult.UpperLimit()
        ExpUpperLimit = hypotestresult.GetExpectedUpperLimit(0)
        ExpUpperLimit_1s = hypotestresult.GetExpectedUpperLimit(1)
        ExpUpperLimit_2s = hypotestresult.GetExpectedUpperLimit(2)
        ExpUpperLimit_m1s = hypotestresult.GetExpectedUpperLimit(-1)
        ExpUpperLimit_m2s = hypotestresult.GetExpectedUpperLimit(-2)

        
        print(mass)
        print("MU up obs = " + str(ObsUpperLimit))
        print("MU up exp = " + str(ExpUpperLimit))
        print('') 
        
        AsymLimits[int(mass)] = ObsUpperLimit
        ExpAsymLimits[int(mass)] = ExpUpperLimit
        ExpAsymLimits_1s[int(mass)] = ExpUpperLimit_1s
        ExpAsymLimits_2s[int(mass)] = ExpUpperLimit_2s
        ExpAsymLimits_m1s[int(mass)] = ExpUpperLimit_m1s
        ExpAsymLimits_m2s[int(mass)] = ExpUpperLimit_m2s

        del AsymCalc
    
    
    # ---------------------------------------
    # ----- repeat the test with Toys MC ----
    # ---------------------------------------
    
    if run_Toys:
        print('')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('TOYS CALCULATION')
        print('')
        
        
        hc = ROOT.RooStats.FrequentistCalculator(data, sbModel, bModel) 
        hc.SetToys(int(NTOYS), int(NTOYS))
        hc.StoreFitInfo(True)
        hc.UseSameAltToys()
        
        # Test statistics: profile likelihood
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(bModel.GetPdf())
        profll.SetMinimizer('Minuit2')
        
        profll.SetOneSidedDiscovery(True)
        profll.SetOneSided(False)
        profll.SetAlwaysReuseNLL(True)
        profll.SetReuseNLL(True)
        
        toymcs = hc.GetTestStatSampler()
        toymcs.SetTestStatistic(profll)
        toymcs_testresult = hc.GetHypoTest()
        toymcs_significance = toymcs_testresult.Significance()
        toymcs_p_value = toymcs_testresult.NullPValue()

        Toys_pvalue[int(mass)] = toymcs_p_value
        Toys_sig[int(mass)] = toymcs_significance

        '''
        Exp_toymcs_pvalue = hc.GetExpectedPValues(toymcs_p_value, toymcs_testresult.AlternatePValue(), 0,False)
        Exp_toymcs_pvalue_1s = hc.GetExpectedPValues(toymcs_p_value, toymcs_testresult.AlternatePValue(), 1,False)
        Exp_toymcs_pvalue_2s = hc.GetExpectedPValues(toymcs_p_value, toymcs_testresult.AlternatePValue(), 2,False)
        Exp_toymcs_pvalue_m1s = hc.GetExpectedPValues(toymcs_p_value, toymcs_testresult.AlternatePValue(), -1,False)                                                                                       
        Exp_toymcs_pvalue_m2s = hc.GetExpectedPValues(toymcs_p_value, toymcs_testresult.AlternatePValue(), -2,False)  

        Exp_toymcs_significance = ROOT.RooStats.PValueToSignificance(Exp_toymcs_pvalue)
        Exp_toymcs_significance_1s = ROOT.RooStats.PValueToSignificance(Exp_toymcs_pvalue_1s)
        Exp_toymcs_significance_2s = ROOT.RooStats.PValueToSignificance(Exp_toymcs_pvalue_2s)
        Exp_toymcs_significance_m1s = ROOT.RooStats.PValueToSignificance(Exp_toymcs_pvalue_m1s) 
        Exp_toymcs_significance_m2s = ROOT.RooStats.PValueToSignificance(Exp_toymcs_pvalue_m2s) 
        
        Toys_exp_pvalue[int(mass)] = Exp_toymcs_pvalue
        Toys_exp_pvalue_1s[int(mass)] = Exp_toymcs_pvalue_1s
        Toys_exp_pvalue_2s[int(mass)] = Exp_toymcs_pvalue_2s
        Toys_exp_pvalue_m1s[int(mass)] = Exp_toymcs_pvalue_m1s   
        Toys_exp_pvalue_m2s[int(mass)] = Exp_toymcs_pvalue_m2s  

        Toys_exp_sig[int(mass)] = Exp_toymcs_significance
        Toys_exp_sig_1s[int(mass)] = Exp_toymcs_significance_1s
        Toys_exp_sig_2s[int(mass)] = Exp_toymcs_significance_2s
        Toys_exp_sig_m1s[int(mass)] = Exp_toymcs_significance_m1s 
        Toys_exp_sig_m2s[int(mass)] = Exp_toymcs_significance_m2s
        '''
            
        print("Observed significance with Toys MC")
        print("For Z'" + mass + " p value --> " + str(toymcs_p_value) + " and Z --> " + str(toymcs_significance))

        del profll
        del hc

        hc = ROOT.RooStats.FrequentistCalculator(data, bModel, sbModel)
        hc.SetToys(int(NTOYS), int(NTOYS))
        hc.StoreFitInfo(True)
        hc.UseSameAltToys()

        # HypoTestInverter for limits extraction
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(sbModel.GetPdf())
        profll.SetMinimizer('Minuit2')
        profll.SetOneSidedDiscovery(False)
        profll.SetAlwaysReuseNLL(True)
        profll.SetReuseNLL(True)

        profll.SetOneSided(True)
        toymcs.SetTestStatistic(profll)
        toymc_calc = ROOT.RooStats.HypoTestInverter(hc)
        toymc_calc.SetConfidenceLevel(0.95)
        toymc_calc.UseCLs(True)
        toymc_calc.SetVerbose(False)
        toymc_calc.SetAutoScan()
        toymc_hypotestresult = toymc_calc.GetInterval()
        
        toymc_ObsUpperLimit = toymc_hypotestresult.UpperLimit()
        toymc_ExpUpperLimit = toymc_hypotestresult.GetExpectedUpperLimit(0)
        toymc_ExpUpperLimit_1s = toymc_hypotestresult.GetExpectedUpperLimit(1)
        toymc_ExpUpperLimit_2s = toymc_hypotestresult.GetExpectedUpperLimit(2)
        toymc_ExpUpperLimit_m1s = toymc_hypotestresult.GetExpectedUpperLimit(-1)
        toymc_ExpUpperLimit_m2s = toymc_hypotestresult.GetExpectedUpperLimit(-2)


        print(mass)
        print("MU up obs = " + str(toymc_ObsUpperLimit))
        print("MU up exp = " + str(toymc_ExpUpperLimit))
        print('')
        
        ToysLimits[int(mass)] = toymc_ObsUpperLimit
        ExpToysLimits[int(mass)] = toymc_ExpUpperLimit
        ExpToysLimits_1s[int(mass)] = toymc_ExpUpperLimit_1s
        ExpToysLimits_2s[int(mass)] = toymc_ExpUpperLimit_2s
        ExpToysLimits_m1s[int(mass)] = toymc_ExpUpperLimit_m1s
        ExpToysLimits_m2s[int(mass)] = toymc_ExpUpperLimit_m2s
        
        del profll
        del hc
     

#order dictionary by masses
Asym_pvalue = collections.OrderedDict(sorted(Asym_pvalue.items()))
Asym_exp_pvalue = collections.OrderedDict(sorted(Asym_exp_pvalue.items()))
Asym_exp_pvalue_1s = collections.OrderedDict(sorted(Asym_exp_pvalue_1s.items()))
Asym_exp_pvalue_2s = collections.OrderedDict(sorted(Asym_exp_pvalue_2s.items()))
Asym_exp_pvalue_m1s = collections.OrderedDict(sorted(Asym_exp_pvalue_m1s.items()))
Asym_exp_pvalue_m2s = collections.OrderedDict(sorted(Asym_exp_pvalue_m2s.items()))

Asym_sig = collections.OrderedDict(sorted(Asym_sig.items()))
Asym_exp_sig = collections.OrderedDict(sorted(Asym_exp_sig.items()))
Asym_exp_sig_1s = collections.OrderedDict(sorted(Asym_exp_sig_1s.items()))
Asym_exp_sig_2s = collections.OrderedDict(sorted(Asym_exp_sig_2s.items()))
Asym_exp_sig_m1s = collections.OrderedDict(sorted(Asym_exp_sig_m1s.items()))
Asym_exp_sig_m2s = collections.OrderedDict(sorted(Asym_exp_sig_m2s.items()))

AsymLimits = collections.OrderedDict(sorted(AsymLimits.items()))
ExpAsymLimits = collections.OrderedDict(sorted(ExpAsymLimits.items()))
ExpAsymLimits_1s = collections.OrderedDict(sorted(ExpAsymLimits_1s.items()))
ExpAsymLimits_2s = collections.OrderedDict(sorted(ExpAsymLimits_2s.items()))
ExpAsymLimits_m1s = collections.OrderedDict(sorted(ExpAsymLimits_m1s.items()))
ExpAsymLimits_m2s = collections.OrderedDict(sorted(ExpAsymLimits_m2s.items()))

Toys_sig = collections.OrderedDict(sorted(Toys_sig.items()))
Toys_pvalue = collections.OrderedDict(sorted(Toys_pvalue.items()))
ToysLimits = collections.OrderedDict(sorted(ToysLimits.items()))
ExpToysLimits = collections.OrderedDict(sorted(ExpToysLimits.items()))
ExpToysLimits_1s = collections.OrderedDict(sorted(ExpToysLimits_1s.items()))
ExpToysLimits_2s = collections.OrderedDict(sorted(ExpToysLimits_2s.items()))
ExpToysLimits_m1s = collections.OrderedDict(sorted(ExpToysLimits_m1s.items()))
ExpToysLimits_m2s = collections.OrderedDict(sorted(ExpToysLimits_m2s.items()))


masses = np.array(Zprime_masses)
xsecs = np.array(Zprime_xsections)

if run_Asimov:
    sigs_asym = np.array(Asym_sig.values())
    pvalues_asym = np.array(Asym_pvalue.values())

    pvalues_asym_exp = np.array(Asym_exp_pvalue.values())
    pvalues_asym_exp_1s = np.array(Asym_exp_pvalue_1s.values())
    pvalues_asym_exp_2s = np.array(Asym_exp_pvalue_2s.values())
    pvalues_asym_exp_m1s = np.array(Asym_exp_pvalue_m1s.values())
    pvalues_asym_exp_m2s = np.array(Asym_exp_pvalue_m2s.values())

    sigs_asym_exp = np.array(Asym_exp_sig.values())
    sigs_asym_exp_1s = np.array(Asym_exp_sig_1s.values())
    sigs_asym_exp_2s = np.array(Asym_exp_sig_2s.values())
    sigs_asym_exp_m1s = np.array(Asym_exp_sig_m1s.values())
    sigs_asym_exp_m2s = np.array(Asym_exp_sig_m2s.values())

    lim_asym = np.array(AsymLimits.values())
 
    exp_lim_asym = np.array(ExpAsymLimits.values())
    exp_lim_1s_asym = np.array(ExpAsymLimits_1s.values())
    exp_lim_2s_asym = np.array(ExpAsymLimits_2s.values())
    exp_lim_m1s_asym = np.array(ExpAsymLimits_m1s.values())
    exp_lim_m2s_asym = np.array(ExpAsymLimits_m2s.values())

    xsec_lim_asym = lim_asym * xsecs
    xsec_exp_lim_asym = exp_lim_asym * xsecs
    xsec_exp_lim_1s_asym = exp_lim_1s_asym * xsecs
    xsec_exp_lim_2s_asym = exp_lim_2s_asym * xsecs
    xsec_exp_lim_m1s_asym = exp_lim_m1s_asym * xsecs
    xsec_exp_lim_m2s_asym = exp_lim_m2s_asym * xsecs


    print('masses')
    print(masses)

    print('Asymptotic formulae')
    print('Observed pvalue')
    print(pvalues_asym)
    print('Expected pvalue')
    print(pvalues_asym_exp)
    print('Observed Significance')
    print(sigs_asym)
    print('Expecteded Significance')
    print(sigs_asym_exp)
    print('Observed mu up')
    print(lim_asym)
    print('Expected mu up')
    print(exp_lim_asym)
    print('Expected mu up + 1 sigma')
    print(exp_lim_1s_asym)
    print('Expected mu up + 2 sigma')
    print(exp_lim_2s_asym)
    print('Expected mu up - 1 sigma')
    print(exp_lim_m1s_asym)
    print('Expected mu up - 2 sigma')
    print(exp_lim_m2s_asym)

    
    plt.fill_between(masses, pvalues_asym_exp_2s, pvalues_asym_exp_m2s, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')     
    plt.fill_between(masses, pvalues_asym_exp_1s, pvalues_asym_exp_m1s, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')
    plt.plot(masses, pvalues_asym_exp, linestyle='dashed', linewidth=0.5, label='Exp')
    plt.plot(masses, pvalues_asym, marker='o', linestyle='dashed', markerfacecolor='black', label='Obs')
    plt.title('Asymptotic approximation')
    plt.ylabel('Null p value')
    plt.legend(loc="upper right")
    plt.savefig("../run/pvalues_Asym.pdf")
    plt.clf()


    plt.fill_between(masses, sigs_asym_exp_2s, sigs_asym_exp_m2s, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')
    plt.fill_between(masses, sigs_asym_exp_1s, sigs_asym_exp_m1s, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')
    plt.plot(masses, sigs_asym_exp, linestyle='dashed', linewidth=0.5, label='Exp')
    plt.plot(masses, sigs_asym, marker='o', linestyle='dashed', markerfacecolor='black', label='Obs')
    plt.title('Asymptotic approximation')
    plt.ylabel('Gaussian significance')
    #plt.yscale('log')
    #plt.ylim(bottom= -0.5)
    plt.legend(loc="upper right")
    plt.savefig("../run/significances_Asym.pdf")
    plt.clf()

    plt.fill_between(masses, xsec_exp_lim_2s_asym, xsec_exp_lim_m2s_asym, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')           
    plt.fill_between(masses, xsec_exp_lim_1s_asym, xsec_exp_lim_m1s_asym, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')                                      
    plt.plot(masses, xsecs, color='red', label='Theory')
    plt.plot(masses, xsec_exp_lim_asym, linestyle='dashed', color='lightslategrey', label='Exp', linewidth=0.5)
    plt.plot(masses, xsec_lim_asym, marker='o', linestyle='dashed', markerfacecolor='black', color='black', label='Obs', markersize=3)
    plt.title('95CLs Upper limits with asymptotic formulae')
    plt.ylabel('Zprime -> tt Xsection [fb]')
    plt.legend(loc="upper right")
    plt.yscale('log')
    #plt.ylim(0, 10000)
    plt.savefig("../run/Limits_Asym.pdf")
    plt.clf()


if run_Toys:
    sigs_toys = np.array(Toys_sig.values())
    pvalues_toys = np.array(Toys_pvalue.values())

    '''
    pvalues_toys_exp = np.array(Toys_exp_pvalue.values())
    pvalues_toys_exp_1s = np.array(Toys_exp_pvalue_1s.values())
    pvalues_toys_exp_2s = np.array(Toys_exp_pvalue_2s.values())
    pvalues_toys_exp_m1s = np.array(Toys_exp_pvalue_m1s.values()) 
    pvalues_toys_exp_m2s = np.array(Toys_exp_pvalue_m2s.values())  

    sigs_toys_exp = np.array(Toys_exp_sig.values())
    sigs_toys_exp_1s = np.array(Toys_exp_sig_1s.values())
    sigs_toys_exp_2s = np.array(Toys_exp_sig_2s.values())
    sigs_toys_exp_m1s = np.array(Toys_exp_sig_m1s.values()) 
    sigs_toys_exp_m2s = np.array(Toys_exp_sig_m2s.values())    
    '''

    lim_toys = np.array(ToysLimits.values())
    exp_lim_toys = np.array(ExpToysLimits.values())
    exp_lim_1s_toys = np.array(ExpToysLimits_1s.values())
    exp_lim_2s_toys = np.array(ExpToysLimits_2s.values())
    exp_lim_m1s_toys = np.array(ExpToysLimits_m1s.values())
    exp_lim_m2s_toys = np.array(ExpToysLimits_m2s.values())

    xsec_lim_toys = lim_toys * xsecs
    xsec_exp_lim_toys = exp_lim_toys * xsecs
    xsec_exp_lim_1s_toys = exp_lim_1s_toys * xsecs
    xsec_exp_lim_2s_toys = exp_lim_2s_toys * xsecs
    xsec_exp_lim_m1s_toys = exp_lim_m1s_toys * xsecs
    xsec_exp_lim_m2s_toys = exp_lim_m2s_toys * xsecs

 
    print('MC Toys')
    print('Observed pvalue')
    print(pvalues_toys)
    print('Observed Significance')
    print(sigs_toys)
    print('Observed mu up')
    print(lim_toys)
    print('Expected mu up')
    print(exp_lim_toys)
    print('Expected mu up + 1 sigma')
    print(exp_lim_1s_toys)
    print('Expected mu up + 2 sigma')
    print(exp_lim_2s_toys)
    print('Expected mu up - 1 sigma')
    print(exp_lim_m1s_toys)
    print('Expected mu up - 2 sigma')
    print(exp_lim_m2s_toys)


    #plt.fill_between(masses, pvalues_toys_exp_2s, pvalues_toys_exp_m2s, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')
    #plt.fill_between(masses, pvalues_toys_exp_1s, pvalues_toys_exp_m1s, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')
    #plt.plot(masses, pvalues_toys_exp, linestyle='dashed', linewidth=0.5, label='Exp')
    plt.plot(masses, pvalues_toys, marker='o', linestyle='dashed', markerfacecolor='black', label='Obs')
    plt.title('Toys MC')
    plt.ylabel('Null p value')
    plt.legend(loc="upper right")
    plt.savefig("../run/pvalues_Toys.pdf")
    plt.clf()

    #plt.fill_between(masses, sigs_toys_exp_2s, sigs_toys_exp_m2s, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')
    #plt.fill_between(masses, sigs_toys_exp_1s, sigs_toys_exp_m1s, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')
    #plt.plot(masses, sigs_toys_exp, linestyle='dashed', linewidth=0.5, label='Exp')
    plt.plot(masses, sigs_toys, marker='o', linestyle='dashed', markerfacecolor='black', label='Obs')
    plt.title('Toys MC')
    plt.ylabel('Gaussian significance')
    #plt.yscale('log')
    #plt.ylim(bottom=0.0001)
    plt.legend(loc="upper right")
    plt.savefig("../run/significances_Toys.pdf")
    plt.clf()

    plt.fill_between(masses, xsec_exp_lim_2s_toys, xsec_exp_lim_m2s_toys, linewidth=0, color='yellow', label='$Exp \pm 2\sigma$')           
    plt.fill_between(masses, xsec_exp_lim_1s_toys, xsec_exp_lim_m1s_toys, linewidth=0, color='lime', label='$Exp \pm 1\sigma$')                                     
    plt.plot(masses, xsecs, color='red', label='Theory')
    plt.plot(masses, xsec_exp_lim_toys, linestyle='dashed', color='lightslategrey', label='Exp',  linewidth=0.5)
    plt.plot(masses, xsec_lim_toys, marker='o', markerfacecolor='black', color='black', label='Obs', markersize=3)
    plt.title('95CLs Upper limits with MC Toys')
    plt.ylabel('Zprime -> tt Xsection [fb]')
    plt.legend(loc="upper right")
    plt.yscale('log')
    plt.savefig("../run/Limits_Toys.pdf")
    plt.clf()



