import ROOT
import os

class BTagSFHandle:
    def __init__(self) :
        self.initialized=False
btagSFhandle = BTagSFHandle()

def initBTagSF () :
# load the BTagCalibrationStandalone.cc macro from https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
    btagSFhandle.initialized=True
    csvpath = "../../HeppyNtuplizer/CMSSW_8_0_26/src/VHbbAnalysis/Heppy/data/csv/"
    ROOT.gSystem.Load(csvpath+'/BTagCalibrationStandalone.so')

# DeepCSV
    btagSFhandle.calib_csv = ROOT.BTagCalibration("DeepCSV", "../data/DeepCSV_2016LegacySF_WP_V1.csv")

# map between algo/flavour and measurement type
    btagSFhandle.sf_type_map = {
    "DeepCSV" : {
        "file" : btagSFhandle.calib_csv,
        "bc" : "comb",
        "l" : "incl",
        },
    }
# map of calibrators. E.g. btag_calibrators["CSVM_nominal_bc"], btag_calibrators["CSVM_up_l"], ...
    btagSFhandle.btag_calibrators = {}
    for algo in ["DeepCSV"]:
        for wp in [ [0, "L"],[1, "M"], [2,"T"] ]:
            for syst in ["central", "up", "down"]:
                for fl in ["bc", "l"]:
                    print "[btagSF]: Loading calibrator for algo:", algo, ", WP:", wp[1], ", systematic:", syst, ", flavour:", fl
    #                btag_calibrators[algo+wp[1]+"_"+syst+"_"+fl].load(sf_type_map[algo]["file"], fl_enum[fl], sf_type_map[algo][fl])
                    btagSFhandle.btag_calibrators[algo+wp[1]+"_"+syst+"_"+fl] = ROOT.BTagCalibrationReader(wp[0], syst)
                    for i in [0,1,2] :
                        btagSFhandle.btag_calibrators[algo+wp[1]+"_"+syst+"_"+fl].load(btagSFhandle.sf_type_map[algo]["file"], i,btagSFhandle.sf_type_map[algo][fl])


# In[9]:


def get_SF(pt=30., eta=0.0, fl=5, val=0.0, syst="central", algo="DeepCSV", wp="M", shape_corr=False, btagSFhandle=btagSFhandle):
    if not btagSFhandle.initialized :
        initBTagSF()
    # no SF for pT<20 GeV or pt>1000 or abs(eta)>2.4
    if abs(eta)>2.4 or pt>1000. or pt<20.:
        return 1.0

    # the .csv files use the convention: b=0, c=1, l=2. Convert into hadronFlavour convention: b=5, c=4, f=0
    fl_index = min(-fl+5,2)
    # no fl=1 in .csv for CMVAv2 (a bug???)
    #if not shape_corr and "CMVAV2" in algo and fl==4:
    #    fl_index = 0

    if shape_corr:
        if applies(fl,syst):
            sf = btagSFhandle.btag_calibrators[algo+"_iterative_"+syst].eval(fl_index ,eta, pt, val)
            #print sf
            return sf
        else:
            sf = btagSFhandle.btag_calibrators[algo+"_iterative_central"].eval(fl_index ,eta, pt, val)
            #print sf
            return sf 

    # pt ranges for bc SF: needed to avoid out_of_range exceptions
    pt_range_high_bc = (670.-1e-02) if "CSV" in algo else (300.-1e-02)
    pt_range_low_bc = 30.+1e-02

    # b or c jets
    if fl>=4:
        # use end_of_range values for pt in [20,30] or pt in [670,1000], with double error
        out_of_range = False
        if pt>pt_range_high_bc or pt<pt_range_low_bc:
            out_of_range = True        
        pt = min(pt, pt_range_high_bc)
        pt = max(pt, pt_range_low_bc)
        sf = btagSFhandle.btag_calibrators[algo+wp+"_"+syst+"_bc"].eval(fl_index ,eta, pt)
        # double the error for pt out-of-range
        if out_of_range and syst in ["up","down"]:
            sf = max(2*sf - btagSFhandle.btag_calibrators[algo+wp+"_central_bc"].eval(fl_index ,eta, pt), 0.)
        #print sf
        return sf
    # light jets
    else:
        sf = btagSFhandle.btag_calibrators[algo+wp+"_"+syst+"_l"].eval( fl_index ,eta, pt)
        #print sf
        return  sf


# In[10]:


debug_btagSF = True
def get_event_SF(jets=[], syst="central", algo="DeepCSV", btagSFhandle=btagSFhandle):
    weight = 1.0
    for jet in jets:
        weight *= get_SF(pt=jet.pt(), eta=jet.eta(), fl=jet.hadronFlavour(), val=(jet.btag("newpfCombinedInclusiveSecondaryVertexV2BJetTags") if algo=="CSV" else jet.btag('newpfCombinedMVAV2BJetTags')), syst=syst, algo=algo, wp="", shape_corr=False, btagSFhandle=btagSFhandle)
    return weight                             

if debug_btagSF:
    print "POG WP:"
    for algo in ["DeepCSV"]:
        for wp in [ "L","M", "T" ]:
            print algo+wp+":"
            for syst in ["central" 
                         , "up", "down"
                         ]:
                print "\t"+syst+":"
                for pt in [19.,25.,31.,330., 680., 1200.]:
                    print ("\t\tB(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=5, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))
                    print ("\t\tC(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=4, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))
                    print ("\t\tL(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=0, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))


# In[ ]:




