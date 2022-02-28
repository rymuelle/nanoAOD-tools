import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
from root_numpy import tree2array

class bffPreselProducer(Module):
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def ptSel(self, obj, variation, met=0):
        if variation=="jerUp": pt = obj.pt_jerUp
        elif variation=="jerDown": pt = obj.pt_jerDown
        elif variation=="jesTotalUp": pt = obj.pt_jesTotalUp
        elif variation=="jesTotalDown": pt = obj.pt_jesTotalDown
        else: 
            if met: pt = obj.pt
            else:
                if self.isMC:
                    pt = obj.pt_nom
                else: pt = obj.pt
        return pt
    def bjetSel(self, jet, variation):
        btagWP = self.btagWP
        pt = self.ptSel(jet,variation)
        return ((pt > 30) & (self.select_btag(jet)) & (abs(jet.eta) < 2.4))
    def lightjetSel(self, jet, variation):
        btagWP = self.btagWP
        pt = self.ptSel(jet,variation)
        return ((pt > 30) & (not self.select_btag(jet)) & (abs(jet.eta) < 2.4) & (jet.jetId > 3) & ((jet.puId & 1) | (pt > 50)))
    def alljetSel(self, jet, variation):
        btagWP = self.btagWP
        return (self.bjetSel(jet, variation) or self.lightjetSel(jet, variation))
    def __init__(self, btagWP, triggers, btag_type="deepcsv", isMC=False, dr_cut=False):
        self.triggers = triggers
        self.btagWP = btagWP
        #select different btags
        def deepcsv(jet):
            return jet.btagDeepB > self.btagWP
        def deepflavour(jet):
            return jet.btagDeepFlavB > self.btagWP
           #set right filtering function
        if btag_type=="deepcsv":
            self.select_btag = deepcsv
        elif btag_type=="DeepFlavour":
            self.select_btag = deepflavour
        print("btag wp: {} type: {}".format(self.btagWP, btag_type))
        self.muSel = lambda x,pt: ((x.corrected_pt > pt) & (abs(x.eta) < 2.4) & (x.tightId > 0) 
                               & (x.pfRelIso04_all < 0.25))
        self.eleSel = lambda x: (x.cutBased_HEEP > 0)
        self.diLepMass = -1
        self.lep_1 = ROOT.TLorentzVector()
        self.lep_2 = ROOT.TLorentzVector()
        self.isMC = isMC
        self.dr_cut = dr_cut
        self.sysDict = {}
        self.sysDict['nom']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"nom"),
        'bjetSel': lambda sel:self.bjetSel(sel,"nom"),
        'alljetSel': lambda sel:self.alljetSel(sel,"nom"),
        'met': lambda sel:self.ptSel(sel,"nom",met=1)}
        pass

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        '''
        In case of data, the b-tagging scale factors are not produced. 
        Check whether they were produced and if not, drop them from jet selections
        '''
        #cutflow histogram
        # all events, after htlt, after two lepton selection, after one or two jets
        self._cutflow_unweighted = ROOT.TH1D('cutflow_unweighted', 'cutflow_unweighted', 4, .5, 4.5)
        self._cutflow_weighted = ROOT.TH1D('cutflow_weighted', 'cutflow_weighted', 4, .5, 4.5)
        list_of_branches = wrappedOutputTree.tree().GetListOfBranches()
        self._triggers = [trigger for trigger in self.triggers if trigger in list_of_branches]
        print(self._triggers)
        if self.isMC:
            self.sysDict['jerUp']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jerUp"),
            'bjetSel': lambda sel:self.bjetSel(sel,"jerUp"),
            'alljetSel': lambda sel:self.alljetSel(sel,"jerUp"),
            'met': lambda sel:self.ptSel(sel,"jerUp", met=1)}
            self.sysDict['jerDown']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jerDown"),
            'bjetSel': lambda sel:self.bjetSel(sel,"jerDown"),
            'alljetSel': lambda sel:self.alljetSel(sel,"jerDown"),
            'met': lambda sel:self.ptSel(sel,"jerDown",met=1)}
            self.sysDict['jesTotalUp']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jesTotalUp"),
            'bjetSel': lambda sel:self.bjetSel(sel,"jesTotalUp"),
            'alljetSel': lambda sel:self.alljetSel(sel,"jesTotalUp"),
            'met': lambda sel:self.ptSel(sel,"jesTotalUp",met=1)}
            self.sysDict['jesTotalDown']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jesTotalDown"),
            'bjetSel': lambda sel:self.bjetSel(sel,"jesTotalDown"),
            'alljetSel': lambda sel:self.alljetSel(sel,"jesTotalDown"),
            'met': lambda sel:self.ptSel(sel,"jesTotalDown",met=1)}
        self.out = wrappedOutputTree
        self.out.branch("inNregions", "F")
        for key in self.sysDict:
            self.out.branch("nBjets_{}".format(key), "F")
            self.out.branch("nSeljets_{}".format(key), "F")
            self.out.branch("JetSFWeight_{}".format(key), "F")
            self.out.branch("HTLT_{}".format(key), "F")
            self.out.branch("RelMET_{}".format(key), "F")
            self.out.branch("TMB_{}".format(key), "F")
            self.out.branch("TMBMin_{}".format(key), "F")
            self.out.branch("TMBMax_{}".format(key), "F")
            self.out.branch("SR1_{}".format(key), "I")
            self.out.branch("CR10_{}".format(key), "I")
            self.out.branch("CR11_{}".format(key), "I")
            self.out.branch("CR12_{}".format(key), "I")
            self.out.branch("CR13_{}".format(key), "I")
            self.out.branch("CR14_{}".format(key), "I")
            self.out.branch("SR2_{}".format(key), "I")
            self.out.branch("CR20_{}".format(key), "I")
            self.out.branch("CR21_{}".format(key), "I")
            self.out.branch("CR22_{}".format(key), "I")
            self.out.branch("CR23_{}".format(key), "I")
            self.out.branch("CR24_{}".format(key), "I")
        self.out.branch("DiLepMass", "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #write cutflow hist
        prevdir = ROOT.gDirectory
        outputFile.cd()
        self._cutflow_weighted.Write()
        self._cutflow_unweighted.Write()
        prevdir.cd()

    def selectDiMu(self, electrons, muons):
        if len(muons) != 2:
            return False
        if len(electrons) != 0:
            return False
        if (muons[0].charge+muons[1].charge) != 0:
            return False
        self.lep_1 = muons[0].p4()*(muons[0].corrected_pt/muons[0].pt)
        self.lep_2 = muons[1].p4()*(muons[1].corrected_pt/muons[1].pt)
        diLep = self.lep_1 + self.lep_2
        self.diLepMass = diLep.M()
        self.out.fillBranch("DiLepMass", self.diLepMass)
        return True
    def selectDiEle(self, electrons):
        if len(electrons) != 2:
            return False
        if (electrons[0].charge+electrons[1].charge) != 0:
            return False
        self.lep_1 = electrons[0].p4()
        self.lep_2 = electrons[1].p4()
        diLep = self.lep_1 + self.lep_2
        self.diLepMass = diLep.M()
        self.out.fillBranch("DiLepMass", self.diLepMass)
        return True
    def selectEleMu(self, electrons, muons):
        if len(electrons) != 1:
            return False
        if len(muons) != 1:
            return False
        if (electrons[0].charge+muons[0].charge) != 0:
            return False
        self.lep_1 = electrons[0].p4()
        self.lep_2 = muons[0].p4()*(muons[0].corrected_pt/muons[0].pt)
        diLep = self.lep_1 + self.lep_2
        self.diLepMass = diLep.M()
        self.out.fillBranch("DiLepMass", self.diLepMass)
        return True
    
    def get_binary_event_weight(self, event):
        #get +/- gen weight as in countHistogramModule
        if hasattr(event, 'Generator_weight') and event.Generator_weight < 0: return -1
        else: return 1
        
    def fill_cutflow(self, value, binary_gen_weight):
        self._cutflow_unweighted.Fill(value)
        self._cutflow_weighted.Fill(value, binary_gen_weight)
        
    def analyze(self, event):
        #get +/- gen weight as in countHistogramModule
        binary_gen_weight = self.get_binary_event_weight(event)
        #cutflow all events
        self.fill_cutflow(1, binary_gen_weight)
        
        ###
        ##HLT
        ###
        HLT_select = False
        for trigger in self._triggers:
            if event[trigger]: 
                HLT_select = 1
                break
        if not HLT_select: return False
        #cutflow after hlt
        self.fill_cutflow(2, binary_gen_weight)
        
        ###
        ##Objects
        ###        
        electrons = sorted(filter(lambda x: self.eleSel(x), Collection(event, "Electron")), key=lambda x: x.pt)
        bJets = sorted(filter(lambda x: self.bjetSel(x, "nom"), Collection(event, "Jet")), key=lambda x: x.pt)
        
        ###
        ##Electrons
        ###   
        if not self.selectDiEle(electrons): return False
        self.fill_cutflow(3, binary_gen_weight)
        
        ###
        ##Jets
        ###   
        if not len(bJets)==1: return False
        self.fill_cutflow(4, binary_gen_weight) 
        
        return True

