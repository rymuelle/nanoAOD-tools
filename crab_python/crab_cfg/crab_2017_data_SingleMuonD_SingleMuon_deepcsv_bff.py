
from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName = '2017_data_SingleMuonD_SingleMuon_deepcsv_bff'
config.General.transferLogs = False
config.General.workArea ='work_areas'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'bash_scripts/bash_script_2017_bff_deepcsv_isMC_False.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['run_processor.py', '../scripts/haddnano.py', 'keep_and_drop_bff.txt']
config.JobType.sendPythonFolder = True
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2017D-02Apr2020-v1/NANOAOD'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
if 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'!='0':
    config.Data.lumiMask = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'

config.Data.outLFNDirBase = '/store/group/phys_exotica/bffZprime/nanoAODskimmed/crab/2017'
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuonD_deepcsv_bff'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
# config.section_('User')
#config.User.voGroup = 'dcms'
config.JobType.allowUndistributedCMSSW = True