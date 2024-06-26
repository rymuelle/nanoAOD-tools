import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
from root_numpy import tree2array

def is_within_deltaR(obj, collection, deltaR=0.4):
    for obj2 in collection:
        if obj2.DeltaR(obj) < deltaR: return True
    return False

def minDR(x,ys):
    minDR = 9999
    for y in ys:
        minDR =min(minDR, x.p4().DeltaR(y.p4()))
    return minDR

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
        elif variation=="jesHEMIssueUp": pt = obj.pt_jesHEMIssueUp
        elif variation=="jesHEMIssueDown": pt = obj.pt_jesHEMIssueDown            
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
        return ((pt > 20) & (self.select_btag(jet)) & (abs(jet.eta) < 2.4) & (jet.jetId > 3) & ((jet.puId & 1) | (pt>50)))
    def lightjetSel(self, jet, variation):
        btagWP = self.btagWP
        pt = self.ptSel(jet,variation)
        return ((pt > 30) & (not self.select_btag(jet)) & (abs(jet.eta) < 2.4) & (jet.jetId > 3) & ((jet.puId & 1) | (pt > 50)))
    def alljetSel(self, jet, variation):
        return (self.bjetSel(jet, variation) or self.lightjetSel(jet, variation))
    def __init__(self, btagWP, triggers, btag_type="deepflavour", isMC=False, dr_cut=False,
                metBranchName='MET', heepBranchName='cutBased_HEEP',
                record_dataframe= False,
                applyHEMfix=False, 
                metBranchPostFix=""):
        self.applyHEMfix = applyHEMfix
        self.record_dataframe = record_dataframe
        self.nselected = 0
        self.metBranchName=metBranchName
        self.triggers = triggers
        self.btagWP = btagWP
        #select different btags
        def deepcsv(jet):
            return jet.btagDeepB > self.btagWP
        def deepflavour(jet):
            return jet.btagDeepFlavB > self.btagWP
           #set right filtering function
        print("btag_type", btag_type)
        if btag_type=="deepcsv":
            self.select_btag = deepcsv
        elif btag_type=="deepflavour":
            self.select_btag = deepflavour
        print(self.select_btag)
        self.btag_type = btag_type
        print("btag wp: {} type: {}".format(self.btagWP, btag_type))
        self.muSel = lambda x,pt: ((x.corrected_pt > pt) & (abs(x.eta) < 2.4) & (x.highPtId > 0) & (x.tkRelIso < .1))
        self.eleSel = lambda x,pt: ((x.pt > pt) & (abs(x.eta) < 2.4) & x[heepBranchName] > 0)
        self.eleSelLowPt = lambda x,pt: ((x.pt > pt) & (abs(x.eta) < 2.4) & x.mvaFall17V2Iso_WPL > 0)
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
        if self.record_dataframe:
            self.list = []
        pass
    def endJob(self):
        if self.record_dataframe:
            import pandas as pd
            df = pd.DataFrame(self.list)
            df.to_csv('event_df.csv')
        print("all selected: {}".format(self.nselected))
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        '''
        In case of data, the b-tagging scale factors are not produced. 
        Check whether they were produced and if not, drop them from jet selections
        '''
        #cutflow histogram
        # all events, after htlt, after two lepton selection, after one or two jets
        self._cutflow_unweighted = ROOT.TH1D('cutflow_unweighted', 'cutflow_unweighted', 4, .5, 4.5)
        self._cutflow_weighted = ROOT.TH1D('cutflow_weighted', 'cutflow_weighted', 4, .5, 4.5)
        self._denis_cutflow_unweighted = ROOT.TH1D('denis_cutflow_unweighted', 'denis_cutflow_unweighted', 10, .5, 10.5)
        self._denis_cutflow_weighted = ROOT.TH1D('denis_cutflow_weighted', 'denis_cutflow_weighted', 10, .5, 10.5)
        list_of_branches = wrappedOutputTree.tree().GetListOfBranches()
        self._triggers = [trigger for trigger in self.triggers if trigger in list_of_branches]
        print(self._triggers)
        if self.isMC:
            if self.applyHEMfix:
                self.sysDict['jesHEMIssueUp']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jesHEMIssueUp"),
                'bjetSel': lambda sel:self.bjetSel(sel,"jesHEMIssueUp"),
                'alljetSel': lambda sel:self.alljetSel(sel,"jesHEMIssueUp"),
                'met': lambda sel:self.ptSel(sel,"jesHEMIssueUp", met=1)}
                self.sysDict['jesHEMIssueDown']= {'lightJetSel': lambda sel:self.lightjetSel(sel,"jesHEMIssueDown"),
                'bjetSel': lambda sel:self.bjetSel(sel,"jesHEMIssueDown"),
                'alljetSel': lambda sel:self.alljetSel(sel,"jesHEMIssueDown"),
                'met': lambda sel:self.ptSel(sel,"jesHEMIssueDown", met=1)}
           
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
        self.out.branch("inNregions", "B")
        self.out.branch("GoodJet", "B", lenVar="nJet")
        self.out.branch("GoodBJet", "B", lenVar="nJet")
        self.out.branch("GoodMuon", "B", lenVar="nMuon")
        self.out.branch("GoodElectron", "B", lenVar="nElectron")
        self.out.branch("GoodMuonLowPt", "B", lenVar="nMuon")
        self.out.branch("GoodElectronLowPt", "B", lenVar="nElectron")
        
        self.out.branch("minGoodJetElDR", "F")
        self.out.branch("minGoodJetMuDR", "F")
        
        self.out.branch("nLep", "I")
        self.out.branch("nLowPtLep", "I")        
        
        self.out.branch("inNregions", "B")
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
        self._denis_cutflow_unweighted.Write()
        self._denis_cutflow_weighted.Write()
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
    def selectDiEle(self, electrons, muons):
        if len(electrons) != 2:
            return False
        if len(muons) != 0:
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
        
    def fill_denis_cutflow(self, value, binary_gen_weight):
        self._denis_cutflow_unweighted.Fill(value)
        self._denis_cutflow_weighted.Fill(value, binary_gen_weight)
        
    def analyze(self, event):
        event_dict = {}
        #get +/- gen weight as in countHistogramModule
        binary_gen_weight = self.get_binary_event_weight(event)
        #cutflow all events
        self.fill_cutflow(1, binary_gen_weight)
        self.fill_denis_cutflow(1, binary_gen_weight)
        HLT_select = False
        for trigger in self._triggers:
            if event[trigger]: 
                HLT_select = 1
                break
        if not HLT_select: return False
        #cutflow after hlt
        self.fill_cutflow(2, binary_gen_weight)
        self.fill_denis_cutflow(2, binary_gen_weight)
        
        jets = Collection(event, "Jet")
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        
        goodJets = [self.alljetSel(jet,"nom") for jet in jets]
        goodBJets = [self.bjetSel(jet,"nom") for jet in jets]
        
        goodMuons = [self.muSel(mu,53) for mu in muons]
        goodMuonsLowPt = [self.muSel(mu,10) for mu in muons]

        goodElectrons = [self.eleSel(el,53) for el in electrons]
        goodElectronsLowPt = [self.eleSelLowPt(el,10) for el in electrons]
        
        jets = sorted(filter(lambda sel: self.alljetSel(sel,"nom"), Collection(event, "Jet")), key=lambda x: self.ptSel(x,'nom'))
        electrons = sorted(filter(lambda x: self.eleSel(x,53), Collection(event, "Electron")), key=lambda x: x.pt)
        muons = sorted(filter(lambda x: self.muSel(x,53), Collection(event, "Muon")), key=lambda x: x.corrected_pt)
        if self.isMC:
            MET = Object(event, "{}_T1Smear".format(self.metBranchName))
        else:
            MET = Object(event, self.metBranchName)
          
        #veto if additional lepton greater than x pt
        lowPtMuonVeto = 10
        electronsLowPt = sorted(filter(lambda x: self.eleSelLowPt(x,lowPtMuonVeto), Collection(event, "Electron")), key=lambda x: x.pt)
        event_dict['nElectronsLowPt'] = len(electronsLowPt)
        #reject low pt electrons that are too close to jets    
        minGoodJetElDR = min([minDR(x, jets) for x in electrons] + [999])
        minGoodJetMuDR = min([minDR(x, jets) for x in muons] + [9999])
        
        electronsLowPt = sorted(filter(lambda x: not is_within_deltaR(x, jets), electronsLowPt))
        event_dict['nElectronsLowPt_post_dr_cut'] = len(electronsLowPt)
        
        muonsLowPt = sorted(filter(lambda x: self.muSel(x,lowPtMuonVeto), Collection(event, "Muon")), key=lambda x: x.corrected_pt)
        nLowPtLep = len(electronsLowPt)+len(muonsLowPt)
        nLep = len(electrons)+len(muons)

        isDiMu = self.selectDiMu(electrons, muons) and nLowPtLep<3
        isDiEle = self.selectDiEle(electrons, muons) and nLowPtLep<3
        isEleMu = self.selectEleMu(electrons, muons) and nLowPtLep<3

        #cutflow after nlep
        self.fill_cutflow(3, binary_gen_weight)
        #Denis specific cuts
        if len(electrons)>=1:
            self.fill_denis_cutflow(3, binary_gen_weight)
            if len(electrons)==2:
                self.fill_denis_cutflow(4, binary_gen_weight)
                if (electrons[0].charge+electrons[1].charge) == 0:
                    self.fill_denis_cutflow(5, binary_gen_weight)
                    if len(muons)==0:
                        self.fill_denis_cutflow(6, binary_gen_weight)
                        if nLowPtLep<3:
                            self.fill_denis_cutflow(7, binary_gen_weight)
        if isDiEle: self.fill_denis_cutflow(8, binary_gen_weight)  
            
        if not (isDiMu or isDiEle or isEleMu):
            return False
      
        eventSelected = False
        for key in self.sysDict:
            lightJetSel = self.sysDict[key]['lightJetSel']
            bjetSel = self.sysDict[key]['bjetSel']
            alljetSel = self.sysDict[key]['alljetSel']
    
            """process event, return True (go to next module) or False (fail, go to next event)"""
            jets = sorted(filter(alljetSel, Collection(event, "Jet")), key=lambda x: self.ptSel(x,key))

            muon_dr = [minDR(jet, muonsLowPt) for jet in jets]
            #questionably useful with deltaR cut above
            electron_dr = [minDR(jet, electronsLowPt) for jet in jets]
            electron_dr += [minDR(jet, electrons) for jet in jets]
            dr_arr = np.array([muon_dr, electron_dr])
            dr_cut_arr =  dr_arr.min(axis=0) > .4
           
            if self.dr_cut:
                jets = np.array(jets)[dr_cut_arr]

            metPt = self.ptSel(MET,key,met=1)

            n_Bjets = len(filter(bjetSel, jets))
            n_lightjets = len(filter(lightJetSel, jets))
            self.out.fillBranch("nBjets_{}".format(key), n_Bjets)
            n_alljets = len(jets)
            self.out.fillBranch("nSeljets_{}".format(key), n_alljets)

            if n_alljets==1 or n_alljets==2:
                eventSelected = True

            jetSFWeight = 1
            if self.isMC:
                for j in jets:
                    jetSFWeight *= j['btagSF_{}_M'.format(self.btag_type)]
            self.out.fillBranch("JetSFWeight_{}".format(key), jetSFWeight)

            htlt = (sum([j.pt for j in jets]) 
                   - sum([ele.pt for ele in electrons]) 
                   - sum([mu.corrected_pt for mu in muons]))
            self.out.fillBranch("HTLT_{}".format(key), htlt)
            
            self.out.fillBranch("RelMET_{}".format(key), metPt/self.diLepMass)
            isSR1  = isDiMu  & (n_Bjets == 1) & (n_lightjets == 0)
            isCR10 = isDiMu  & (n_Bjets == 0) & (n_lightjets == 1)
            isCR11 = isEleMu & (n_Bjets == 1) & (n_lightjets == 0)
            isCR12 = isEleMu & (n_Bjets == 0) & (n_lightjets == 1)
            isCR13 = isDiEle & (n_Bjets == 1) & (n_lightjets == 0)
            isCR14 = isDiEle & (n_Bjets == 0) & (n_lightjets == 1)
            isSR2  = isDiMu  & (n_Bjets >= 1) & (n_alljets == 2)
            isCR20 = isDiMu  & (n_Bjets == 0) & (n_alljets == 2)
            isCR21 = isEleMu & (n_Bjets >= 1) & (n_alljets == 2)
            isCR22 = isEleMu & (n_Bjets == 0) & (n_alljets == 2)
            isCR23 = isDiEle & (n_Bjets >= 1) & (n_alljets == 2)
            isCR24 = isDiEle & (n_Bjets == 0) & (n_alljets == 2)

            self.out.fillBranch("SR1_{}".format(key), isSR1)
            self.out.fillBranch("CR10_{}".format(key), isCR10)
            self.out.fillBranch("CR11_{}".format(key), isCR11)
            self.out.fillBranch("CR12_{}".format(key), isCR12)
            self.out.fillBranch("CR13_{}".format(key), isCR13)
            self.out.fillBranch("CR14_{}".format(key), isCR14)
            self.out.fillBranch("SR2_{}".format(key), isSR2)
            self.out.fillBranch("CR20_{}".format(key), isCR20)
            self.out.fillBranch("CR21_{}".format(key), isCR21)
            self.out.fillBranch("CR22_{}".format(key), isCR22)
            self.out.fillBranch("CR23_{}".format(key), isCR23)
            self.out.fillBranch("CR24_{}".format(key), isCR24)
            #should not be more than one
            if key == "nom": 
                self.out.fillBranch("inNregions",  + isSR1 + isCR10 + isCR11 + isCR12 + isCR13 + isCR14 + isSR2 + isCR20 + isCR21 + isCR22 + isCR23 + isCR24)
                
            sbm = -1.
            sbmMin = -1.
            sbmMax = -1.

            if len(jets) == 1:
                sbm1 = (self.lep_1 + jets[0].p4()).M()
                sbm2 = (self.lep_2 + jets[0].p4()).M()
                sbm = min(sbm1, sbm2)
                sbmMin = sbm
                sbmMax = max(sbm1, sbm2)
            elif len(jets) == 2:
                l1j1mass = (self.lep_1 + jets[0].p4()).M()
                l1j2mass = (self.lep_1 + jets[1].p4()).M()
                l2j1mass = (self.lep_2 + jets[0].p4()).M()
                l2j2mass = (self.lep_2 + jets[1].p4()).M()
                if abs(l1j1mass - l2j2mass) < abs(l2j1mass - l1j2mass):
                    sbm = max(l1j1mass, l2j2mass)
                    sbmMin = min(l1j1mass, l2j2mass)
                    sbmMax = sbm
                else:
                    sbm = max(l1j2mass, l2j1mass)
                    sbmMin = min(l1j2mass, l2j1mass)
                    sbmMax = sbm
            else:
                return False
            if key == "nom": 
                #cutflow after njets
                self.fill_cutflow(4, binary_gen_weight)
                if isDiEle  & (n_Bjets == 1): self.fill_denis_cutflow(9, binary_gen_weight)
                if isCR13: self.fill_denis_cutflow(10, binary_gen_weight)
            self.out.fillBranch("TMB_{}".format(key), sbm)
            self.out.fillBranch("TMBMin_{}".format(key), sbmMin)
            self.out.fillBranch("TMBMax_{}".format(key), sbmMax)
            self.out.fillBranch("GoodJet", goodJets)
            self.out.fillBranch("GoodBJet", goodBJets)
            self.out.fillBranch("GoodMuon", goodMuons)
            self.out.fillBranch("GoodElectron", goodElectrons)
            self.out.fillBranch("GoodMuonLowPt", goodMuonsLowPt)
            self.out.fillBranch("GoodElectronLowPt", goodElectronsLowPt)
            self.out.fillBranch("minGoodJetElDR", minGoodJetElDR)
            self.out.fillBranch("minGoodJetMuDR", minGoodJetMuDR)
            self.out.fillBranch("nLowPtLep", nLowPtLep)
            self.out.fillBranch("nLep", nLep)            
        
        nLowPtLep
        if eventSelected: 
            if self.record_dataframe:
                self.list.append(event_dict)
            self.nselected+=1
            return True
        else: return True