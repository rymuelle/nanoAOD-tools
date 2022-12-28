import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class triggerFilter(Module):

    def __init__(self, triggers):
        self.triggers = triggers
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        list_of_branches = wrappedOutputTree.tree().GetListOfBranches()
        self._triggers = [trigger for trigger in self.triggers if trigger in list_of_branches]
        print(self._triggers)
        #set up branches for triggers missing in file
        # needed for RDF since it expects all files to have all triggers
        self.missing_triggers = [trigger for trigger in self.triggers if not trigger in list_of_branches]
        self.out = wrappedOutputTree
        for mt in self.missing_triggers:
            print(mt)
            self.out.branch(mt, "B")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        #set missing triggers to 0
        # might not be needed
        for mt in self.missing_triggers:
            self.out.fillBranch(mt, False)
        HLT_select = False
        for trigger in self._triggers:
            if event[trigger]: 
                HLT_select = 1
                break
        if not HLT_select: return False
        else: return True

