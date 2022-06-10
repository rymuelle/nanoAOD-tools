
from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName = '2018_data_SingleMuonA_SingleMuon_deepflavour_bff'
config.General.transferLogs = False
config.General.workArea ='work_areas'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'bash_scripts/bash_script_2018_bff_deepflavour_isMC_False.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['run_processor.py', '../scripts/haddnano.py', 'keep_and_drop_bff.txt']
config.JobType.sendPythonFolder = True
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2018A-02Apr2020-v1/NANOAOD'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
if 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'!='0':
    config.Data.lumiMask = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'

config.Data.outLFNDirBase = '/store/group/phys_exotica/bffZprime/nanoAODskimmed/crab/2018'
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuonA_deepflavour_bff'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
# config.section_('User')
#config.User.voGroup = 'dcms'
config.JobType.allowUndistributedCMSSW = True