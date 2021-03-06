import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
  #'file:/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/RootFiles/MC_GENOnly/BToKMuMu_MC_OnlyGEN_8TeV.root',
  'file:/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/aachen/work/CMSSW_5_3_20/src/ana/BToKMuMu/python/BToKMuMu_MC_OnlyGEN_8TeV.root',
    
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')r
#process.GlobalTag.globaltag = cms.string('START53_V23::All')
#process.GlobalTag.globaltag = cms.string('START53_V7G::All')
process.GlobalTag.globaltag = cms.string('START53_V19F::All')

process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(True)
process.p = cms.Path(process.ntuple)
