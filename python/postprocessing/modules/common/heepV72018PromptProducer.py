import ROOT
from ROOT.Math import PtEtaPhiMVector
import math

ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def set_bit(bitmap, bit, value):
    '''take bitmap and change bit position to 1 or 0 based on bool of value'''
    if value: return bitmap | (1<<bit)
    else: return bitmap & ~(1<<bit)
    
def Log2(x): 
    return math.log(x)/math.log(2)

def isInteger(x):
    return int(x) == x

def nBit(n):
    log2_value = Log2(n)
    assert isInteger(log2_value), "n is not a power of 2"
    return int(log2_value)

def get_rho(event):
    rhoBranchName = "fixedGridRhoFastjetAll"
    rho = getattr(event, rhoBranchName)
    pv = Object(event, "PV")
    return pv.npvs

def is_endcap_electron(electron):
    '''define if electron is in endcap.
    pg. 6 https://indico.cern.ch/event/787315/contributions/3434898/attachments/1847223/3031172/HEEP_2019_0517_v3.pdf'''
    return abs(electron.eta) > 1.566 and abs(electron.eta) < 2.5
    
# cuts as described on page 9 https://indico.cern.ch/event/831669/contributions/3485543/attachments/1871797/3084930/ApprovalSlides_EE_v3.pdf

def heepV72018Prompt_HOE_cut(electron, energy, rho):
    is_endcap = is_endcap_electron(electron)
    if not is_endcap: 
        return electron.hoe < 1.0 / electron.eCorr +0.05
    else: 
        energy_numerator = (-0.4 + 0.4 * abs(electron.eta)) * rho
        return electron.hoe < energy_numerator / energy + 0.05
        
def em_p_had_iso(electron):
    return electron.dr03EcalRecHitSumEt + electron.dr03HcalDepth1TowerSumEt

def heepV72018Prompt_em_p_had_iso_cut(electron, Et, rho):
    em_p_had_iso_value = em_p_had_iso(electron)
    is_endcap = is_endcap_electron(electron)
    if not is_endcap:
        return em_p_had_iso_value < 2 + 0.03 * Et + 0.28 * rho
    else:
        rho_coef = (0.15 + 0.07 * abs(electron.eta))
        if Et < 50:
            return em_p_had_iso_value < 2.5 +  rho_coef * rho
        else:
            return em_p_had_iso_value < 2.5 + 0.03 * (Et - 50) + rho_coef * rho
        
def passes_heepV72018Prompt_cuts(electron, rho):
    #define electron lorentz vector
    electron_LV = PtEtaPhiMVector(electron.pt, electron.eta, electron.phi, electron.mass)
    #compute new cuts
    hoe_cut_value = heepV72018Prompt_HOE_cut(electron, electron_LV.E(), rho)
    em_p_had_iso_cut_value = heepV72018Prompt_em_p_had_iso_cut(electron, electron_LV.Et(), rho)
    # update bitmap with new cuts
    # GSfEleEMHadD1IsoRhoCut (256, 0x100) and GsfEleHadronicOverEMLinearCut (64, 0x40) 
    # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Applying_Individual_Cuts_of_a_Se
    vidNestedWPBitmapHEEP = electron.vidNestedWPBitmapHEEP
    vidNestedWPBitmapHEEP = set_bit(electron.vidNestedWPBitmapHEEP, nBit(64), hoe_cut_value)
    vidNestedWPBitmapHEEP = set_bit(electron.vidNestedWPBitmapHEEP, nBit(256), em_p_had_iso_cut_value)    
    return vidNestedWPBitmapHEEP==4095

class heepV72018PromptProducer(Module):
    def __init__(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Electron_cutBased_HEEPV7p0_2018Prompt", "B", lenVar="nElectron")
    def beginJob(self):
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        # get rho
        rho = get_rho(event)
        #apply cut
        Electron_cutBased_HEEPV7p0_2018Prompt = list(map(lambda x: passes_heepV72018Prompt_cuts(x, rho), electrons))
        #fill new branch
        self.out.fillBranch("Electron_cutBased_HEEPV7p0_2018Prompt", Electron_cutBased_HEEPV7p0_2018Prompt)
        return True
