import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
from root_numpy import tree2array

class bffPreselProducer(Module):
    def __init__(self, btagWP, triggers, 
                 isMC=True,
                muon_pt = "corrected_pt", 
                jet_sys = "",
                ele_pt = "pt", 
                metBranchName='MET', 
                metBranchPostFix='', 
                heepBranchName='cutBased_HEEP',
                record_dataframe= False,
                applyHEMfix=False,
                btag_type='',
                ):
        self.btagWP = btagWP
        self.triggers = triggers
        self.isMC = isMC
        #object variation definitions
        self.muon_pt = muon_pt
        self.muon_pt_name = "_muon_{}".format(muon_pt) if len(muon_pt)>0 else ""
        self.jet_sys = jet_sys
        self.jet_sys_name = "_jet_{}".format(jet_sys) if len(jet_sys)>0 else "_jet_nom"
        self.jet_pt_key = 'pt_' + self.jet_sys if len(jet_sys)>0 else "pt"
        self.ele_pt = ele_pt
        self.ele_pt_name = "_ele_{}".format(ele_pt) if len(ele_pt)>0 else ""
        self.applyHEMfix = applyHEMfix
        self.metBranchName = metBranchName
        self.metBranchPostFix = metBranchPostFix
        self.met_pt_name = "pt_{}".format(jet_sys) if len(jet_sys)>0 else "pt"
        self.heepBranchName = heepBranchName
        self.regions = ["SR1","CR10","CR11","CR12","CR13","CR14",
                        "SR2","CR20","CR21","CR22","CR23","CR24"]
        self.key = self.jet_sys_name + self.muon_pt_name + self.ele_pt_name
        if self.key == "": self.key = "_nom"
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        key = self.key
        print(self.key)        
        self.out.branch("GoodMuon{}".format(key), "B", lenVar="nMuon")
        self.out.branch("GoodElectron{}".format(key), "B", lenVar="nElectron")
        self.out.branch("GoodMuonLowPt{}".format(key), "B", lenVar="nMuon")
        self.out.branch("GoodElectronLowPt{}".format(key), "B", lenVar="nElectron")
        self.out.branch("nLep{}".format(key), "I")
        self.out.branch("nLowPtLep{}".format(key), "I")
        self.out.branch("minGoodJetElDR{}".format(key), "F")
        self.out.branch("minGoodJetMuDR{}".format(key), "F")
        self.out.branch("inNregions{}".format(key), "B")
        self.out.branch("GoodJet{}".format(key), "B", lenVar="nJet")
        self.out.branch("GoodLJet{}".format(key), "B", lenVar="nJet")
        self.out.branch("GoodBJet{}".format(key), "B", lenVar="nJet")
        self.out.branch("nBjets{}".format(key), "F")
        self.out.branch("nSeljets{}".format(key), "F")
        self.out.branch("HTLT{}".format(key), "F")
        self.out.branch("RelMET{}".format(key), "F")
        self.out.branch("TMB{}".format(key), "F")
        self.out.branch("TMBMin{}".format(key), "F")
        self.out.branch("TMBMax{}".format(key), "F")
        for reg in self.regions:
            self.out.branch("{}{}".format(reg, key), "I")
        self.out.branch("DiLepMass{}".format(key), "F")
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def analyze(self, event):
        ##
        ## apply first cuts
        ##
        #get event weight
        binary_gen_weight = get_binary_event_weight(event)
        #apply trigger
        event_triggered = apply_trigger(event, self.triggers)
        if not event_triggered: return False
        ##
        ## Get and select objects
        ##
        muons, goodMuons, nMuons = prepare_collection(event, "Muon", muSel, self.muon_pt)
        eles, goodEle, nEles = prepare_collection(event, "Electron", lambda x: eleSel(x, self.heepBranchName), self.ele_pt)
        muonsLowPt, goodMuonsLowPt, nMuonsLowPt = prepare_collection(event, "Muon", muSelLowPt, self.muon_pt)
        elesLowPt, goodElectronsLowPt, nElesLowPt = prepare_collection(event, "Electron", lambda x: eleSelLowPt(x, self.heepBranchName), self.ele_pt)
        leptons = elesLowPt + muonsLowPt
        nLep = nMuons + nEles
        nLepLowPt = nElesLowPt + nMuonsLowPt
        #apply lepton selection
        if (nLep ==2 and nLepLowPt==2) and (leptons[0].charge + leptons[1].charge == 0):
            diLepMass = (leptons[0].p4()+leptons[1].p4()).M()
        else:
            diLepMass = -1 
        #
        #Select Jets, reject if within .04 of a lepton
        #
        #make collection that is not delta r filtered
        jets_no_dr, _, _ = prepare_collection(event, "Jet", lambda x: alljetSel(x, self.btagWP), self.jet_pt_key)
        
        jets, goodJets, nJets = prepare_collection(event, "Jet", lambda x: alljetSel(x, self.btagWP), self.jet_pt_key, 
        dr_collection=leptons)
        ljets, goodLJets, nLJets = prepare_collection(event, "Jet", lambda x: lightjetSel(x, self.btagWP), self.jet_pt_key, 
        dr_collection=leptons)
        bJets, goodBJets, nBJets = prepare_collection(event, "Jet", lambda x: bjetSel(x, self.btagWP), self.jet_pt_key, 
        dr_collection=leptons)
        #use smeared MET for mc
        MET = Object(event, "{}{}".format(self.metBranchName, self.metBranchPostFix))
        MET.pt = MET[self.met_pt_name]
        ###
        ### Get and select objects
        ###
        # lepton veto
        minGoodJetElDR = min([minDR(x, jets_no_dr) for x in eles] + [999])
        minGoodJetMuDR = min([minDR(x, jets_no_dr) for x in muons] + [9999])
        minJetLepDR = min(minGoodJetElDR,minGoodJetMuDR)
        if diLepMass > 0: 
            inRegion, region_dict = calc_regions(nMuons, nEles, nBJets, nLJets, nJets)    
            # select regions
            relmet, htlt, tmb, tmbMin, tmbMax = calculate_bff_variables(jets, leptons, MET.pt, diLepMass)
        else:
            inRegion = 0
            region_dict = {reg:0 for reg in self.regions}
            relmet, htlt, tmb, tmbMin, tmbMax = -1,-1,-1,-1,-1
        ##
        ## fill branches
        ##
        key = self.key
        self.out.fillBranch("GoodMuon{}".format(key), goodMuons)
        self.out.fillBranch("GoodElectron{}".format(key), goodEle)
        self.out.fillBranch("GoodMuonLowPt{}".format(key), goodMuonsLowPt)
        self.out.fillBranch("GoodElectronLowPt{}".format(key), goodElectronsLowPt)
        self.out.fillBranch("nLep{}".format(key), nLep)
        self.out.fillBranch("nLowPtLep{}".format(key), nLepLowPt)
        self.out.fillBranch("minGoodJetElDR{}".format(key), minGoodJetElDR)
        self.out.fillBranch("minGoodJetMuDR{}".format(key), minGoodJetMuDR)
        self.out.fillBranch("inNregions{}".format(key), inRegion)
        self.out.fillBranch("GoodJet{}".format(key), goodJets)
        self.out.fillBranch("GoodLJet{}".format(key), goodLJets)
        self.out.fillBranch("GoodBJet{}".format(key), goodBJets)
        self.out.fillBranch("nBjets{}".format(key), nBJets)
        self.out.fillBranch("nSeljets{}".format(key), nJets)
        self.out.fillBranch("HTLT{}".format(key), htlt)
        self.out.fillBranch("RelMET{}".format(key), relmet)
        self.out.fillBranch("TMB{}".format(key), tmb)
        self.out.fillBranch("TMBMin{}".format(key), tmbMin)
        self.out.fillBranch("TMBMax{}".format(key), tmbMax)
        for reg in self.regions:
            self.out.fillBranch("{}{}".format(reg, key), int(region_dict[reg]))
        self.out.fillBranch("DiLepMass{}".format(key), diLepMass)        
        return True
    
    
def calc_regions(nMuons, nEles, nBJets, nLJets, nJets):
    isDiMu = (nMuons==2) & (nEles==0)
    isEleMu = (nMuons==1) & (nEles==1)
    isDiEle = (nMuons==0) & (nEles==2)
    region_dict = {
    "SR1"  : isDiMu  & (nBJets == 1) & (nLJets == 0),
    "CR10" : isDiMu  & (nBJets == 0) & (nLJets == 1),
    "CR11" : isEleMu & (nBJets == 1) & (nLJets == 0),
    "CR12" : isEleMu & (nBJets == 0) & (nLJets == 1),
    "CR13" : isDiEle & (nBJets == 1) & (nLJets == 0),
    "CR14" : isDiEle & (nBJets == 0) & (nLJets == 1),
    "SR2"  : isDiMu  & (nBJets >= 1) & (nJets == 2),
    "CR20" : isDiMu  & (nBJets == 0) & (nJets == 2),
    "CR21" : isEleMu & (nBJets >= 1) & (nJets == 2),
    "CR22" : isEleMu & (nBJets == 0) & (nJets == 2),
    "CR23" : isDiEle & (nBJets >= 1) & (nJets == 2),
    "CR24" : isDiEle & (nBJets == 0) & (nJets == 2),
    }  
    inRegion = sum([v for k,v in region_dict.items()])
    return inRegion, region_dict
    
    
def calculate_bff_variables(jets, leptons, met_pt, dilepmass):
    relmet = met_pt/dilepmass
    htlt = (sum([j.pt for j in jets]) 
           - sum([lep.pt for lep in leptons]))
    #TMB
    sbm = -1.
    sbmMin = -1.
    sbmMax = -1.
    if len(jets) == 1:
        sbm1 = (leptons[0].p4() + jets[0].p4()).M()
        sbm2 = (leptons[1].p4() + jets[0].p4()).M()
        sbm = min(sbm1, sbm2)
        sbmMin = sbm
        sbmMax = max(sbm1, sbm2)
    elif len(jets) == 2:
        l1j1mass = (leptons[0].p4() + jets[0].p4()).M()
        l1j2mass = (leptons[0].p4() + jets[1].p4()).M()
        l2j1mass = (leptons[1].p4() + jets[0].p4()).M()
        l2j2mass = (leptons[1].p4() + jets[1].p4()).M()
        if abs(l1j1mass - l2j2mass) < abs(l2j1mass - l1j2mass):
            sbm = max(l1j1mass, l2j2mass)
            sbmMin = min(l1j1mass, l2j2mass)
            sbmMax = sbm
        else:
            sbm = max(l1j2mass, l2j1mass)
            sbmMin = min(l1j2mass, l2j1mass)
            sbmMax = sbm
    return relmet, htlt, sbm, sbmMin, sbmMax
        
def prepare_collection(event, collection, filter_func, pt_key, dr_collection=[]):
    collection_obj = Collection(event, collection)
    collection_obj = list(collection_obj)
    collection_obj_new = []
    for i, x in enumerate(collection_obj):
        #remove objects that are within some deltar of some object
        if len(dr_collection)>0:
            if minDR(x,dr_collection) < .4: 
                continue
        x.pt = x[pt_key]
        collection_obj_new.append(x)
    goodObj = [filter_func(x) for x in collection_obj_new]
    filtered_obj = filter_collection(collection_obj_new,goodObj)
    return filtered_obj, goodObj, len(filtered_obj)

def seleDiLep(leptons):
    if len(leptons)!=2: return False
    if (leptons[0].charge+leptons[1].charge) != 0: return False

def apply_trigger(event, triggers):
        HLT_select = False
        for trigger in triggers:
            if event[trigger]: 
                HLT_select = 1
                break
        if not HLT_select: return False
        return True

def get_binary_event_weight(event):
    #get +/- gen weight as in countHistogramModule
    if hasattr(event, 'Generator_weight') and event.Generator_weight < 0: return -1
    else: return 1
        
def minDR(x,ys):
    minDR = 9999
    for y in ys:
        minDR =min(minDR, x.p4().DeltaR(y.p4()))
    return minDR

def filter_collection(coll, mask):
    return [x for x, m in zip(coll, mask) if m]

# jet selections
def bjetSel(jet, btagWP, leptons = []):
    return ((jet.pt > 20) & 
            (jet['btagDeepFlavB'] > btagWP) & 
            (abs(jet.eta) < 2.4) & 
            (jet.jetId > 3) & 
            ((jet.puId & 1) | (jet.pt>50)))
def lightjetSel(jet, btagWP):
    return ((jet.pt > 30) & 
            (not jet['btagDeepFlavB'] > btagWP) 
            & (abs(jet.eta) < 2.4) 
            & (jet.jetId > 3) 
            & ((jet.puId & 1) | (jet.pt > 50)))
def alljetSel(jet, btagWP):
    return (bjetSel(jet, btagWP) or lightjetSel(jet, btagWP))

# lep selection
def muSel(mu):
    pt = 53
    return (mu.pt > pt) & (abs(mu.eta) < 2.4) & (mu.highPtId > 0) & (mu.tkRelIso < .1)
def muSelLowPt(mu):
    pt = 10
    return (mu.pt > pt) & (abs(mu.eta) < 2.4) & (mu.highPtId > 0) & (mu.tkRelIso < .1)
def eleSel(el, heepBranchName):
    pt = 53
    return (el.pt > pt) & (abs(el.eta) < 2.4) & (el[heepBranchName] > 0)
def eleSelLowPt(el, heepBranchName):
    pt = 10
    #make sure to count muons in previous selection
    eleSelection = eleSel(el, heepBranchName)
    return eleSelection or ((el.pt > pt) & (abs(el.eta) < 2.4) & (el.mvaFall17V2Iso_WPL > 0))