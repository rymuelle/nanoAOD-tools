import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import numpy as np
from array import array

class bffBtagEffProducer(Module):
    def __init__(self, btagWP, triggers, btag_type="deepflavour", **kwargs):
        self.isMC = True
        self.btagWP = btagWP
        def deepcsv(jet):
            return jet.btagDeepB > self.btagWP
        def deepflavour(jet):
            return jet.btagDeepFlavB > self.btagWP
           #set right filtering function
        if btag_type=="deepcsv":
            self.select_btag = deepcsv
        elif btag_type=="deepflavour":
            self.select_btag = deepflavour
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
        return ((pt > 20) & (self.select_btag(jet)) & (abs(jet.eta) < 2.4) & (jet.jetId > 3) & ((jet.puId & 1) | (pt>50)))
    def lightjetSel(self, jet, variation):
        btagWP = self.btagWP
        pt = self.ptSel(jet,variation)
        return ((pt > 30) & (not self.select_btag(jet)) & (abs(jet.eta) < 2.4) & (jet.jetId > 3) & ((jet.puId & 1) | (pt > 50)))
    def alljetSel(self, jet, variation):
        return (self.bjetSel(jet, variation) or self.lightjetSel(jet, variation))

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        '''
        In case of data, the b-tagging scale factors are not produced. 
        Check whether they were produced and if not, drop them from jet selections
        '''
            
        flavorBins = array('d',[-0.5,  0.5,  1.5,  2.5])
        ptBins = array('d',[  20.,   30.,   50.,   70.,  100.,  140.,  200.,  300.,  600.,1000., 7000.])
        self.bTagEffTH2F = ROOT.TH2F('bTagEff', 'bTagEff', 
                              3, flavorBins,
                              10, ptBins
                                 )
        self.TotalTH2F = ROOT.TH2F('total', 'total', 
                              3, flavorBins,
                              10, ptBins
                                 )                                                        
                              
                              
        self.out = wrappedOutputTree
        self.Pass = ROOT.std.vector(ROOT.std.vector('int'))()
        self.Total= ROOT.std.vector(ROOT.std.vector('int'))()
        self.out._tree.Branch("PassBtag", self.Pass)
        self.out._tree.Branch("TotalBtag",self.Total)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        self.bTagEffTH2F.Divide(self.TotalTH2F)
        self.bTagEffTH2F.Write()
        self.TotalTH2F.Write()
        prevdir.cd()
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = sorted(filter(lambda x: self.alljetSel(x,'nom'), Collection(event, "Jet")), key=lambda x: x.pt)
        jetPt = [20.,30.,50.,70.,100.,140.,200.,300.,600.,1000.]
        self.Pass.clear()
        Dummy1 = ROOT.std.vector('int')()
        Dummy1.resize(10,0)
        self.Pass.resize(3,Dummy1)
        self.Total.clear()
        Dummy2 = ROOT.std.vector('int')()
        Dummy2.resize(10,0)
        self.Total.resize(3,Dummy2)

        if self.isMC:
            for j in jets:
                # select flavor
                flavour = 2
                if j.genJetIdx>=0:
                    if j.hadronFlavour==5:
                        flavour=0
                    elif j.hadronFlavour==4:
                        flavour=1  
                if self.select_btag(j): 
                      self.bTagEffTH2F.Fill(flavour, j.pt)
                
                self.TotalTH2F.Fill(flavour, j.pt)        
                        
        if self.isMC:
            for j in jets:
                bin = 9
                for k in reversed(jetPt):
                    if j.pt>k:
                        break
                    bin-=1
                if bin>=0:
                    flavor = 2
                    if j.genJetIdx>=0:
                        if j.hadronFlavour==5:
                            flavor=0
                        elif j.hadronFlavour==4:
                            flavor=1
                    if self.select_btag(j): 
                        self.Pass[flavor][bin] += 1
                    self.Total[flavor][bin]+= 1
	else:
            pass
        self.out._tree.Fill()

        return True
