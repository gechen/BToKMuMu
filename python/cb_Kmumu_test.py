from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'BuToKMuMu_MC_8TeV_LooseFilter_v4'

config.section_('JobType')
config.JobType.psetName = './btokmumu_MC_8TeV.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['BToKMuMu.root']

config.section_('Data')
#config.Data.inputDataset = '/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/Summer12_DR53X-PU_RD2_START53_V19F-v1/AODSIM'
config.Data.inputDataset = '/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/Summer12DR53X-PU_RD2_START53_V19F-v1/AODSIM'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '/store/user/gechen/'
#config.Data.publishDBS = 'phys03'
#config.Data.publication = True
config.Data.splitting = 'FileBased'
config.Data.outputDatasetTag = 'Bu2MuMuK_LooseFilter'

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


