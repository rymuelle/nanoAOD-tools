
from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName = '2017_BFF_BFF_350_dbs1p0_BFFZprimeToMuMu_M_350_dbs1p0_deepflavour_bff'
config.General.transferLogs = False
config.General.workArea ='work_areas'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'bash_scripts/bash_script_2017_bff_deepflavour_isMC_True.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['run_processor.py', '../scripts/haddnano.py', 'keep_and_drop_bff.txt']
config.JobType.sendPythonFolder = True
config.section_('Data')
config.Data.inputDataset = '/BFFZprimeToMuMu_M_350_dbs1p0/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/NANOAODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
if '0'!='0':
    config.Data.lumiMask = '0'

config.Data.outLFNDirBase = '/store/group/phys_exotica/bffZprime/nanoAODskimmed/crab/2017'
config.Data.publication = False
config.Data.outputDatasetTag = 'BFF_350_dbs1p0_deepflavour_bff'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
# config.section_('User')
#config.User.voGroup = 'dcms'
config.JobType.allowUndistributedCMSSW = True