from WMCore.Configuration import Configuration
config = Configuration()
#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'BuToKMuMu_Data_2012C_8TeV_v8'
config.General.transferLogs = False

config.section_('JobType')
config.JobType.psetName = './btokmumu_Data_8TeV.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.outputFiles = ['BuToKMuMu_Data_2012C_8TeV.root']
config.JobType.outputFiles = ['BToKMuMu.root']

config.section_('Data')
config.Data.inputDataset = '/MuOniaParked/Run2012C-22Jan2013-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.lumiMask = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/gechen/crab3_run/'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'BuToKMuMu_Data_2012C_8TeV_v8'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
