from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'BuToJpsiK_MC_GENOnly_8TeV_Ntuples_v5'
config.section_('JobType')
config.JobType.psetName = './BuToJpsiK_MC_GENOnly.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'privateMC'
#config.JobType.outputFiles = ['BToJpsiK_GENOnly_8TeV_Ntuple.root']
config.section_('Data')
config.Data.inputDataset = '/PYTHIA6_BuToJpsiK_GENOnly_8TeV/gechen-crab_BuToJpsiKMuMu_MC_GENOnly_8TeV-387bf2b3df13ffa8b4f3dd9f3950e077/USER'
#config.Data.outputPrimaryDataset = 'PYTHIA6_BuToJpsiK_GENOnly_8TeV_Ntuple_v4'
config.Data.outputDatasetTag = 'PYTHIA6_BuToJpsiK_GENOnly_8TeV_Ntuple_v5'
config.Data.unitsPerJob = 2
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = True
config.Data.publishDBS = 'phys03'
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/gechen/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ["T2_CH*"]




