#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.trigger.triggerFilter import triggerFilter
from PhysicsTools.NanoAODTools.postprocessing.bff.bffPreselModule import bffPreselProducer 
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsProducer 
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF 
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import lepSF2016
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016 
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import muonScaleRes2016 
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis

isMC = True
dataYear = "2016"
runPeriod = ""

triggers= ['HLT_Mu50','HLT_TkMu50', 'HLT_DoubleEle33_CaloIdL_MW', 'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL'] 

jmeCorrections = createJMECorrector(
    isMC=isMC, 
    dataYear=dataYear, 
    runPeriod=runPeriod,
    jesUncert="Total", 
    applySmearing=True,
    jetType="AK4PFchs",
    noGroom=False)

#keep_and_drop = '${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/scripts/keep_and_drop_bff.txt'
keep_and_drop = 'keep_and_drop_bff.txt'

modules=[
    countHistogramsProducer(),
    triggerFilter(triggers),
    btagSF(dataYear),
    jmeCorrections(),
    puWeight_2016(),
    muonScaleRes2016(),
    lepSF2016(),
    bffPreselProducer(int(dataYear), isMC=True)
    ]

p = PostProcessor(".",
                  inputFiles(),
                  modules=modules,
                  provenance=True,
                  fwkJobReport=True,
                  outputbranchsel=keep_and_drop,
                  )
p.run()

print("DONE")
#note: lepSF has to go before puWeight due to the lepton weight constructor redefining the weight constructor
#-I PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer btagSF$era \