import ROOT
import sys


poiValueForBackground = 0

Zprime_masses = ["400", "500", "750", "1000", "1250", "1500", "1750", "2000", "2250", "2500", "2750", "3000"]
#Zprime_masses = ["1000", "1250", "1500"]


for mass in Zprime_masses:

    # Open models files                                                                                                                                                                                   
    f = ROOT.TFile.Open("../run/Model_" + mass + ".root")
    
    #retrieve the workspace
    
    w = f.Get("w_" + mass)
    w.Print() 

    data = w.data("data_hist") 
    mc = w.obj("ModelConfig_" + mass)
    model = w.pdf("model " + mass)

    # make sure ingredients are found
    if not data or not mc or not model:   
        print("data or ROOT.RooStats.ModelConfig or model was not found")
        exit

    POI = mc.GetParametersOfInterest().first()

    # -------------------------------------------------
    # ----- Configure a ProfileLikelihoodCalculator -----
    # -------------------------------------------------

    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, model, ROOT.RooArgSet(POI))

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

    print("For Z'" + mass + " p value --> " + str(p_value) + " and Z --> " + str(significance))

    del plc
