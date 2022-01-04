import ROOT

#path_trees = '/data2/sauricchio/corso_AnalisiDati_PhD/Data/1largeRjet1lep/' 
path_trees = 'https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/1largeRjet1lep/'

samples = ['data', 'tt', 'Vjets', 'single_top', 'Diboson', 'ZPrime400', 'ZPrime500', 'ZPrime750', 'ZPrime1000', 'ZPrime1250', 'ZPrime1500', 'ZPrime1750', 'ZPrime2000', 'ZPrime2250', 'ZPrime2500', 'ZPrime2750', 'ZPrime3000']
#samples = ['single_top']

Lumi = 10000 #pb-1

ROOT.gInterpreter.Declare(
    """
    float Define_Weight(float xsection, float sumweights, float mcweight, float SF_ele, float SF_muon, float SF_lepTrigg, float SF_pu, float SF_btag) 
    {      
        float SF = SF_ele*SF_muon*SF_lepTrigg*SF_pu*SF_btag;
        float weight = SF*mcweight/sumweights;
        weight = weight*xsection;
        
        return weight;                                                                                                                                           
    }
    """
    )


ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;
    using Vec_unint = const ROOT::VecOps::RVec<unsigned int>;
    using Vec_bool = const ROOT::VecOps::RVec<bool>;
    int GoodLep(UInt_t lepsize, Vec_t& pt, Vec_t& eta, Vec_t& phi, Vec_t& e, Vec_t& ptcone30, Vec_t& etcone20, Vec_t& trackd0pvunbiased, Vec_t& tracksigd0pvunbiased, Vec_t& z0, Vec_unint& leptype, Vec_bool& istight) 
    {
        int n_goodleps = 0;
        int index = -9;
        
        for(int i=0; i<lepsize; i++)
        {
            ROOT::Math::PtEtaPhiEVector temp_lep(pt[i]/1000., eta[i], phi[i], e[i]/1000.);
            if(istight[i]) 
            {
                if((pt[i]>30000.) && ((ptcone30[i]/pt[i])<0.15) && ((etcone20[i]/pt[i])<0.15))
                {
                    if((leptype[i]==11) && (TMath::Abs(eta[i])<2.47) && ((TMath::Abs(eta[i])<1.37) || (TMath::Abs(eta[i])>1.52)))
                    {
                        if((TMath::Abs(trackd0pvunbiased[i])/tracksigd0pvunbiased[i]<5) && (TMath::Abs(z0[i]*TMath::Sin(temp_lep.Theta()))<0.5))
                        {
                            index = i;
                            n_goodleps++;
                            }
                        }

                    if((leptype[i]==13) && (TMath::Abs(eta[i])<2.5))
                    {
                        if((TMath::Abs(trackd0pvunbiased[i])/tracksigd0pvunbiased[i]<3) && ((TMath::Abs(z0[i])*TMath::Sin(temp_lep.Theta()))<0.5))
                        {
                            index=i;
                            n_goodleps++;
                            }
                         }
                     }
                  }
            }

        if(n_goodleps!=1) return -9; //save only index for events with exactly 1 good lepton
        return index;
         }
    """
 )

ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;  
    float Calculate_Wtransvmass(int goodlep_index, Vec_t& lep_pt, Vec_t& lep_eta, Vec_t& lep_phi, Vec_t& lep_e, Float_t met_et, Float_t  met_phi)
    {
        ROOT::Math::PtEtaPhiEVector lepton(lep_pt[goodlep_index], lep_eta[goodlep_index], lep_phi[goodlep_index], lep_e[goodlep_index]);
        ROOT::Math::PtEtaPhiEVector MeT(met_et, 0, met_phi, met_et);
        float mtw = TMath::Sqrt(2*lepton.Pt()*met_et*(1-TMath::Cos(lepton.Phi() - MeT.Phi())));
        return mtw;
        }
    """
)


ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;                                                                                                                                                         
    using Vec_i = ROOT::VecOps::RVec<int>; 
    auto GoodBJets(Vec_t& jet_pt, Vec_t& jet_eta, Vec_t& jet_jvt, Vec_t& jet_MV2c10)                                                                                                       
    {                                                                                                                                                                                                      
        int nbjets=0;
        Vec_i goodbjets_index;

        for(int i=0; i<jet_pt.size(); i++)                                                                                                                                                                         
        {                                                                                                                                                                                                  
            if((jet_pt[i]>30000.) && (TMath::Abs(jet_eta[i])<2.5))                                                                                                                                         
            {                                                                                                                                                                                              
                bool jvt_pass=true;                                                                                                                                                                        
                if ((jet_pt[i]<60000.) && (TMath::Abs(jet_eta[i])<2.4) && (jet_jvt[i]<0.59)) jvt_pass=false;                                                                                               
                if (jvt_pass)                                                                                                                                                                              
                {                                                                                                                                                                                          
                    if(jet_MV2c10[i]>0.8244273)
                    {
                        goodbjets_index.push_back(i);
                        nbjets++;
                        }                                                                                                                                                  
                                                                                                                                                                                                          
                    }                                                                                                                                                                                      
                }                                                                                                                                                                                 
             }         
         return goodbjets_index;                                                                                                                                                                   
        
}                                                                                                                                                                                                  
   """
)


ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;  
    using Vec_i = ROOT::VecOps::RVec<int>; 

    auto GoodJet(UInt_t jet_n, Vec_t& jet_pt, Vec_t& jet_eta, Vec_t& jet_phi, Vec_t& jet_E, Vec_t& jet_jvt, int goodlep_index, Vec_t& lep_pt, Vec_t& lep_eta, Vec_t& lep_phi, Vec_t& lep_e)                     {
        ROOT::Math::PtEtaPhiEVector Lepton(lep_pt[goodlep_index], lep_eta[goodlep_index], lep_phi[goodlep_index], lep_e[goodlep_index]);   

        Vec_i GoodJet_index;
        int NGoodjets=0;
        
        for(int i=0; i<jet_n; i++)                                                                                                                                                                         
        {                                                                                                                                                                                                  
            if((jet_pt[i]>30000.) && (TMath::Abs(jet_eta[i])<2.5))                                                                                                                                         
            {                                                                                                                                                                                              
                bool jvt_pass=true;                                                                                                                                                                        
                if ((jet_pt[i]<60000.) && (TMath::Abs(jet_eta[i])<2.4) && (jet_jvt[i]<0.59)) jvt_pass=false;                                                                                               
                if (jvt_pass)                                                                                                                                                                              
                {                                                                                                                                                                                          
                    ROOT::Math::PtEtaPhiEVector temp_jet(jet_pt[i], jet_eta[i], jet_phi[i], jet_E[i]);
                    float dR_lep_jet = sqrt((Lepton.Eta() - temp_jet.Eta())*(Lepton.Eta() - temp_jet.Eta()) + (Lepton.Phi() - temp_jet.Phi())*(Lepton.Phi() - temp_jet.Phi()));
                    if(dR_lep_jet<2.0) 
                    {
                        GoodJet_index.push_back(i);
                        NGoodjets++;
                        }
                    }                                                                                                                                                                                      
                }    
            }
    
        if(NGoodjets==0) GoodJet_index.push_back(-9);
        return GoodJet_index[0];
        }
    """
    )
ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>; 
    using Vec_int = const ROOT::VecOps::RVec<int>; 
    int GoodTopLRJet(int goodlep_index, Vec_t& lep_pt, Vec_t& lep_eta, Vec_t& lep_phi, Vec_t& lep_e, int goodjet_index, Vec_t& jet_pt, Vec_t& jet_eta, Vec_t& jet_phi, Vec_t& jet_E, UInt_t NlargeRjets, Vec_t& LRjet_m, Vec_t& LRjet_pt, Vec_t& LRjet_eta, Vec_t& LRjet_phi, Vec_t& LRjet_E, Vec_t& LRjet_tau32, int Ngoodbjets, Vec_int& goodbjets_indexes)
    { 

        ROOT::Math::PtEtaPhiEVector Lepton(lep_pt[goodlep_index], lep_eta[goodlep_index], lep_phi[goodlep_index], lep_e[goodlep_index]);  
        ROOT::Math::PtEtaPhiEVector goodsmallRjet(jet_pt[goodjet_index], jet_eta[goodjet_index], jet_phi[goodjet_index], jet_E[goodjet_index]);

        int NTopLRjets=0;
        int TopLRjets_index=-9;

        for(int l=0; l<NlargeRjets; l++)
        {
            if (((LRjet_m[l]/1000.)>50) && ((LRjet_pt[l]/1000.)>300) && ((LRjet_pt[l]/1000.)<1500) && (TMath::Abs(LRjet_eta[l])<2.0) && (LRjet_tau32[l]<0.75))
            {
                ROOT::Math::PtEtaPhiEVector TopLRjet(LRjet_pt[l], LRjet_eta[l], LRjet_phi[l], LRjet_E[l]);
                float dPhi = TMath::Abs(Lepton.Phi() - TopLRjet.Phi());
                dPhi  = dPhi < TMath::Pi() ? dPhi : 2*TMath::Pi() - dPhi;
                
                float dR = sqrt((goodsmallRjet.Eta() - TopLRjet.Eta())*(goodsmallRjet.Eta() - TopLRjet.Eta()) + (goodsmallRjet.Phi() - TopLRjet.Phi())*(goodsmallRjet.Phi() - TopLRjet.Phi()));
                
                bool btag_is_within_TopLR = false;
                bool btag_is_smallRjet = false; 
                
                for(int b=0; b<Ngoodbjets; b++)
                {
                    ROOT::Math::PtEtaPhiEVector bjet(jet_pt[goodbjets_indexes[b]], jet_eta[goodbjets_indexes[b]], jet_phi[goodbjets_indexes[b]], jet_E[goodbjets_indexes[b]]);
                    float dR_bjet_TopLR = sqrt((bjet.Eta() - TopLRjet.Eta())*(bjet.Eta() - TopLRjet.Eta()) + (bjet.Phi() - TopLRjet.Phi())*(bjet.Phi() - TopLRjet.Phi()));               
                    float dR_bjet_smallRjet = sqrt((bjet.Eta() - goodsmallRjet.Eta())*(bjet.Eta() - goodsmallRjet.Eta()) + (bjet.Phi() - goodsmallRjet.Phi())*(bjet.Phi() - goodsmallRjet.Phi()));   
                    
                    if(dR_bjet_TopLR<1.0) btag_is_within_TopLR=true;     
                    if(dR_bjet_smallRjet<0.01) btag_is_smallRjet=true;
                    
                    }
              
                if((dR>1.5) && (dPhi>1.0) && (btag_is_within_TopLR || btag_is_smallRjet))
                {
                    NTopLRjets++;
                    TopLRjets_index=l;
                    }
                
                }
            
            }
        if (NTopLRjets==1) return TopLRjets_index; //return the index of good top large-R jet only for events with exactly 1 good top large-R jet
        else return -9;
        }
    """
)


ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;
    float FinalInvMass(int goodlep_index, Vec_t& lep_pt, Vec_t& lep_eta, Vec_t& lep_phi, Vec_t& lep_e, int goodjet_index, Vec_t& jet_pt, Vec_t& jet_eta, Vec_t& jet_phi, Vec_t& jet_E, int TopLRjets_index, Vec_t& LRjet_pt, Vec_t& LRjet_eta, Vec_t& LRjet_phi, Vec_t& LRjet_E, Vec_t& jet_MV2c10)
    {   
        ROOT::Math::PtEtaPhiEVector Lepton(lep_pt[goodlep_index], lep_eta[goodlep_index], lep_phi[goodlep_index], lep_e[goodlep_index]);                                                                   
        ROOT::Math::PtEtaPhiEVector goodsmallRjet(jet_pt[goodjet_index], jet_eta[goodjet_index], jet_phi[goodjet_index], jet_E[goodjet_index]);     

        ROOT::Math::PtEtaPhiEVector TopHad(LRjet_pt[TopLRjets_index], LRjet_eta[TopLRjets_index], LRjet_phi[TopLRjets_index], LRjet_E[TopLRjets_index]);
        ROOT::Math::PtEtaPhiEVector TopLep, TTbar;

        float m_tt=-99.;

        if(jet_MV2c10[goodjet_index]>0.8244273) 
        {
            TopLep = goodsmallRjet + Lepton;
            TTbar = TopHad + TopLep; 
            m_tt = TTbar.M()/1000.;
            }
        return m_tt;

        }
    """
    )
    

file_histos =  ROOT.TFile('file_histos.root', 'RECREATE')

for s in samples:
    ntuples = ROOT.TChain('mini')

    if s=='data':
        path = path_trees + 'Data/'
        ntuples_list = ['data_A.1largeRjet1lep.root', 'data_B.1largeRjet1lep.root', 'data_C.1largeRjet1lep.root', 'data_D.1largeRjet1lep.root']

    else:
        path = path_trees + 'MC/'
        
        if s=='tt':
            ntuples_list = ['mc_410000.ttbar_lep.1largeRjet1lep.root']
        if s=='Vjets':
            ntuples_list = ['mc_361106.Zee.1largeRjet1lep.root', 'mc_361107.Zmumu.1largeRjet1lep.root', 'mc_361108.Ztautau.1largeRjet1lep.root', 
                            'mc_364156.Wmunu_PTV0_70_CVetoBVeto.1largeRjet1lep.root', 'mc_364158.Wmunu_PTV0_70_BFilter.1largeRjet1lep.root', 'mc_364157.Wmunu_PTV0_70_CFilterBVeto.1largeRjet1lep.root', 'mc_364159.Wmunu_PTV70_140_CVetoBVeto.1largeRjet1lep.root', 'mc_364160.Wmunu_PTV70_140_CFilterBVeto.1largeRjet1lep.root', 'mc_364161.Wmunu_PTV70_140_BFilter.1largeRjet1lep.root', 'mc_364162.Wmunu_PTV140_280_CVetoBVeto.1largeRjet1lep.root', 'mc_364163.Wmunu_PTV140_280_CFilterBVeto.1largeRjet1lep.root', 'mc_364164.Wmunu_PTV140_280_BFilter.1largeRjet1lep.root',  'mc_364165.Wmunu_PTV280_500_CVetoBVeto.1largeRjet1lep.root', 'mc_364166.Wmunu_PTV280_500_CFilterBVeto.1largeRjet1lep.root', 'mc_364167.Wmunu_PTV280_500_BFilter.1largeRjet1lep.root', 'mc_364168.Wmunu_PTV500_1000.1largeRjet1lep.root', 'mc_364169.Wmunu_PTV1000_E_CMS.1largeRjet1lep.root',
                            'mc_364170.Wenu_PTV0_70_CVetoBVeto.1largeRjet1lep.root', 'mc_364171.Wenu_PTV0_70_CFilterBVeto.1largeRjet1lep.root', 'mc_364172.Wenu_PTV0_70_BFilter.1largeRjet1lep.root', 'mc_364175.Wenu_PTV70_140_BFilter.1largeRjet1lep.root', 'mc_364173.Wenu_PTV70_140_CVetoBVeto.1largeRjet1lep.root', 'mc_364174.Wenu_PTV70_140_CFilterBVeto.1largeRjet1lep.root', 'mc_364176.Wenu_PTV140_280_CVetoBVeto.1largeRjet1lep.root', 'mc_364177.Wenu_PTV140_280_CFilterBVeto.1largeRjet1lep.root', 'mc_364178.Wenu_PTV140_280_BFilter.1largeRjet1lep.root', 'mc_364179.Wenu_PTV280_500_CVetoBVeto.1largeRjet1lep.root', 'mc_364177.Wenu_PTV140_280_CFilterBVeto.1largeRjet1lep.root', 'mc_364178.Wenu_PTV140_280_BFilter.1largeRjet1lep.root', 'mc_364182.Wenu_PTV500_1000.1largeRjet1lep.root', 'mc_364183.Wenu_PTV1000_E_CMS.1largeRjet1lep.root', 
                            'mc_364184.Wtaunu_PTV0_70_CVetoBVeto.1largeRjet1lep.root', 'mc_364186.Wtaunu_PTV0_70_BFilter.1largeRjet1lep.root', 'mc_364185.Wtaunu_PTV0_70_CFilterBVeto.1largeRjet1lep.root', 'mc_364187.Wtaunu_PTV70_140_CVetoBVeto.1largeRjet1lep.root', 'mc_364188.Wtaunu_PTV70_140_CFilterBVeto.1largeRjet1lep.root', 'mc_364189.Wtaunu_PTV70_140_BFilter.1largeRjet1lep.root', 'mc_364190.Wtaunu_PTV140_280_CVetoBVeto.1largeRjet1lep.root', 'mc_364192.Wtaunu_PTV140_280_BFilter.1largeRjet1lep.root', 'mc_364191.Wtaunu_PTV140_280_CFilterBVeto.1largeRjet1lep.root', 'mc_364193.Wtaunu_PTV280_500_CVetoBVeto.1largeRjet1lep.root', 'mc_364194.Wtaunu_PTV280_500_CFilterBVeto.1largeRjet1lep.root', 'mc_364195.Wtaunu_PTV280_500_BFilter.1largeRjet1lep.root', 'mc_364196.Wtaunu_PTV500_1000.1largeRjet1lep.root', 'mc_364197.Wtaunu_PTV1000_E_CMS.1largeRjet1lep.root', 
            ]
        if s=='single_top':
            ntuples_list = ['mc_410011.single_top_tchan.1largeRjet1lep.root', 'mc_410012.single_antitop_tchan.1largeRjet1lep.root', 'mc_410025.single_top_schan.1largeRjet1lep.root', 'mc_410026.single_antitop_schan.1largeRjet1lep.root', 'mc_410013.single_top_wtchan.1largeRjet1lep.root', 'mc_410014.single_antitop_wtchan.1largeRjet1lep.root']
        if s=='Diboson':
            ntuples_list = ['mc_363356.ZqqZll.1largeRjet1lep.root', 'mc_363358.WqqZll.1largeRjet1lep.root', 'mc_363360.WplvWmqq.1largeRjet1lep.root', 'mc_363489.WlvZqq.1largeRjet1lep.root', 'mc_363491.lllv.1largeRjet1lep.root', 'mc_363492.llvv.1largeRjet1lep.root', 'mc_363493.lvvv.1largeRjet1lep.root']
        if s=='ZPrime400':
            ntuples_list = ['mc_301322.ZPrime400_tt.1largeRjet1lep.root']
        if s=='ZPrime500':
            ntuples_list = ['mc_301323.ZPrime500_tt.1largeRjet1lep.root']
        if s=='ZPrime750':
            ntuples_list = ['mc_301324.ZPrime750_tt.1largeRjet1lep.root']
        if s=='ZPrime1000':
            ntuples_list = ['mc_301325.ZPrime1000_tt.1largeRjet1lep.root']
        if s=='ZPrime1250':
            ntuples_list = ['mc_301326.ZPrime1250_tt.1largeRjet1lep.root']
        if s=='ZPrime1500':
            ntuples_list = ['mc_301327.ZPrime1500_tt.1largeRjet1lep.root']
        if s=='ZPrime1750':
            ntuples_list = ['mc_301328.ZPrime1750_tt.1largeRjet1lep.root']
        if s=='ZPrime2000':
            ntuples_list = ['mc_301329.ZPrime2000_tt.1largeRjet1lep.root']
        if s=='ZPrime2250':
            ntuples_list = ['mc_301330.ZPrime2250_tt.1largeRjet1lep.root']
        if s=='ZPrime2500':
            ntuples_list = ['mc_301331.ZPrime2500_tt.1largeRjet1lep.root']
        if s=='ZPrime2750':
            ntuples_list = ['mc_301332.ZPrime2750_tt.1largeRjet1lep.root']
        if s=='ZPrime3000':
            ntuples_list = ['mc_301333.ZPrime3000_tt.1largeRjet1lep.root']

            

    for n in ntuples_list:
        ntuples.Add(path + n)

    print(' ') 
    print(s)
    #print(ntuples.Print())
    print(' ')

    initial_df = ROOT.RDataFrame(ntuples)

    
    ####~~~~~  Selections  ~~~~~####


    #if 'single' in s:
    #    initial_df = initial_df.Define('mcweight', "return (mcWeight/TMath::Abs(mcWeight));")    
    #    weighted_df = initial_df.Define('weight', 'Define_Weight(XSection, SumWeights, mcweight, scaleFactor_ELE, scaleFactor_MUON, scaleFactor_LepTRIGGER, scaleFactor_PILEUP, scaleFactor_BTAG)')
    #else:
    weighted_df = initial_df.Define('weight', 'Define_Weight(XSection, SumWeights, mcWeight, scaleFactor_ELE, scaleFactor_MUON, scaleFactor_LepTRIGGER, scaleFactor_PILEUP, scaleFactor_BTAG)')

    #apply trigger request
    #select only events with at least one large_r jet and apply the cut MET>20GeV

    pre_sel_df = weighted_df.Filter("((trigE==true) || (trigM==true)) && (largeRjet_n>=1) && (met_et > 20000.)")


    #take events with ONLY 1 good lepton passing quality and isolation criteria
    good_lep_df = pre_sel_df.Define('goodLep_index', 'GoodLep(lep_n, lep_pt, lep_eta, lep_phi, lep_E, lep_ptcone30, lep_etcone20, lep_trackd0pvunbiased, lep_tracksigd0pvunbiased, lep_z0, lep_type, lep_isTightID)').Filter("goodLep_index!=(-9)")
    
    #MET + mTW > 60 GeV
    met_df = good_lep_df.Define('MTW', 'Calculate_Wtransvmass(goodLep_index, lep_pt, lep_eta, lep_phi, lep_E, met_et, met_phi)').Filter("(met_et+MTW)>60000.")
    
    #define goodbjets 
    good_bjets_df = met_df.Define('goodbjets_indexes', 'GoodBJets(jet_pt, jet_eta, jet_jvt, jet_MV2c10)')
    good_bjets_df = good_bjets_df.Define('Ngoodbjets', "return goodbjets_indexes.size();")
    
    #good small-r jets close to the lepton (the leading one is the b-tagged candidate from t->bW decay)
    good_bjets_df = good_bjets_df.Define('GoodJet_index', 'GoodJet(jet_n, jet_pt, jet_eta, jet_phi, jet_E, jet_jvt, goodLep_index, lep_pt, lep_eta, lep_phi, lep_E)')

    #select only events with >=1 b tagged jet and >=1 good small-r jet 
    good_bjets_df = good_bjets_df.Filter("(Ngoodbjets>0) && (GoodJet_index!=(-9))")
    
    #find the large-R jet from the top hadronic decay
    goodTopjet_df = good_bjets_df.Define('TopLRjets_index', 'GoodTopLRJet(goodLep_index, lep_pt, lep_eta, lep_phi, lep_E, GoodJet_index, jet_pt, jet_eta, jet_phi, jet_E, largeRjet_n, largeRjet_m, largeRjet_pt, largeRjet_eta, largeRjet_phi, largeRjet_E, largeRjet_tau32, Ngoodbjets, goodbjets_indexes)')

    #require exactly one good top large-R jet 
    goodTopjet_df = goodTopjet_df.Filter("TopLRjets_index!=-9")

    #Construct the invariant mass with top-tagged large-R jet + small-R jet (b-tagged) + lepton
    final_df = goodTopjet_df.Define('TTbar_M', 'FinalInvMass(goodLep_index, lep_pt, lep_eta, lep_phi, lep_E, GoodJet_index, jet_pt, jet_eta, jet_phi, jet_E, TopLRjets_index, largeRjet_pt, largeRjet_eta, largeRjet_phi, largeRjet_E, jet_MV2c10)')

    #remove events where small-R jet is not b-tagged 
    final_df = final_df.Filter('TTbar_M>0')

    nbins = 20
    a = 400
    b = 1600
    
    if s!='data':
        h_mass = final_df.Histo1D(ROOT.RDF.TH1DModel('ttbar_mass_' + s, 'ttbar_mass', nbins, a, b), "TTbar_M", "weight")
        h_mass.Scale(Lumi, "width")
    else:
        h_mass = final_df.Histo1D(ROOT.RDF.TH1DModel('ttbar_mass_' + s, 'ttbar_mass', nbins, a, b), "TTbar_M")
        h_mass.Scale(1, "width")

    h_mass.GetXaxis().SetTitle("m^{Z->tt}_{visible} [GeV]")
    h_mass.GetYaxis().SetTitle("Entries/60 GeV")
    
    print(h_mass.Integral("width"))
    
    #h_mass.Draw()
    h_mass.Write()
    
file_histos.Close()


#to draw with entries/bin_size do h_mass.Scale(1, "width")
#save DFs in files and save histograms
#another python script which open the histos, draw TTbar_M for data, signal and bkgs. With bks create the stack plot and draw data/MC, then create the full_bkg histogram. Then save full_bkg, data and signals in a file
#create code for statistical analysis
