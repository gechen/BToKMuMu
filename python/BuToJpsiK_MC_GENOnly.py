import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'root://xrootd.unl.edu//store/user/gechen/PYTHIA6_BuToJpsiK_GENOnly_8TeV/crab_BuToJpsiKMuMu_MC_GENOnly_8TeV/160529_220333/0000/BToJpsiK_MC_OnlyGEN_8TeV_1.root',
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')r
#process.GlobalTag.globaltag = cms.string('START53_V23::All')
#process.GlobalTag.globaltag = cms.string('START53_V7G::All')
process.GlobalTag.globaltag = cms.string('START53_V19F::All')

process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(True)
process.p = cms.Path(process.ntuple)
