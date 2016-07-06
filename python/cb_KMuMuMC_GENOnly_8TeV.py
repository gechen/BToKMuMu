from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'BuToKMuMu_MC_GENOnly_8TeV_v4'
config.section_('JobType')
config.JobType.psetName = './PYTHIA6_Bu2MuMuK_TuneZ2star_8TeV_GENOnly.py'
config.JobType.pluginName = 'privateMC'
#config.JobType.outputFiles = ['BToJpsiKMuMu_MC_OnlyGEN_8TeV.root']
config.JobType.generator = 'pythia'
config.section_('Data')
#config.Data.outputDatasetTag = 'PYTHIA6_BuToJpsiK_TuneZ2star_GENOnly_8TeV'
config.Data.outputPrimaryDataset = 'PYTHIA6_Bu2MuMuK_GENOnly_8TeV'
config.Data.publishDBS = 'phys03'
config.Data.publication = True
config.Data.unitsPerJob = 1000000
config.Data.splitting = 'EventBased'
config.Data.totalUnits = 200000000
config.section_('User')
config.section_('Site')
config.Site.whitelist = ['T2_US_FNAL']
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T2_US_FNAL'
