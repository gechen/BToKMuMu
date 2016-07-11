import FWCore.ParameterSet.Config as cms
from btokmumu_2012_cfi import process

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'root://xrootd.unl.edu//store/mc/Summer12DR53X/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/60C7051C-68A9-E411-9FFA-0025907FD428.root',
#'root://xrootd.unl.edu//store/mc/Summer12DR53X/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/84D2C859-65A9-E411-B1C9-002590DB923C.root',
#'root://xrootd.unl.edu//store/mc/Summer12DR53X/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/C05863CD-71A9-E411-B86F-002590DB9262.root',
#'root://xrootd.unl.edu//store/mc/Summer12DR53X/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/2C4C5D3F-65A9-E411-A977-002590DB924E.root',
'root://xrootd.unl.edu//store/mc/Summer12_DR53X/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/000771A0-9408-E411-94C0-002590DB92A8.root',
#'/store/mc/Summer12DR53X/BuToMuMuK_MuFilterLoose_TuneZ2star_8TeV-pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/2C4C5D3F-65A9-E411-A977-002590DB924E.root',
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/0021569C-5709-E411-912B-003048D3739A.root'
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/000771A0-9408-E411-94C0-002590DB92A8.root',
#'root://eoscms.cern.ch//store/mc/Summer12_DR53X/Bu2MuMuK_TuneZ2star_8TeV_Pythia6/AODSIM/PU_RD2_START53_V19F-v1/00000/000771A0-9408-E411-94C0-002590DB92A8.root',
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BdToKstarMuMu_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/00100578-164A-E311-80CC-485B39800C31.root',
#'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BpToPsiMuMu_2MuPtEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v1/0000/00388EB8-FDF9-E111-83E3-0030487F1BD7.root'
#'file:/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/aachen/work/test/CMSSW_5_3_20/src/PYTHIA6_Bu2MuMuK_TuneZ2star_8TeV_cff_py_RAW2DIGI_L1Reco_RECO.root',
#    '/store/mc/Summer12_DR53X/BuToMuMuK_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/20000/9C4AFE3E-BA65-E211-821F-AC162DACC3E8.root',
    
                            )
                        )
#process.GlobalTag.globaltag = cms.string('GR_R_42_V25::All')
#process.GlobalTag.globaltag = cms.string('START53_V23::All')
#process.GlobalTag.globaltag = cms.string('START53_V7G::All')
process.GlobalTag.globaltag = cms.string('START53_V19F::All')

# do trigger matching for muons
triggerProcessName = 'HLT'


process.cleanMuonTriggerMatchHLT0 = cms.EDProducer(
            # match by DeltaR only (best match by DeltaR)
            'PATTriggerMatcherDRLessByR',
                        src                   = cms.InputTag('cleanPatMuons'),
                        # default producer label as defined in
                        # PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
                        matched               = cms.InputTag('patTrigger'),
                        matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced*",0,0)'),
                        maxDeltaR             = cms.double(0.1),
                        # only one match per trigger object
                        resolveAmbiguities    = cms.bool(True),
                        # take best match found per reco object (by DeltaR here, see above)
                        resolveByMatchQuality = cms.bool(False))


from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT0'],
                               hltProcess = triggerProcessName, outputModule = '')

g_TriggerNames_LastFilterNames = [
            ('HLT_DoubleMu3p5_LowMass_Displaced',  'hltDisplacedmumuFilterDoubleMu3p5LowMass')
            ]

g_TriggerNames = [i[0] for i in g_TriggerNames_LastFilterNames]
g_LastFilterNames = [i[1] for i in g_TriggerNames_LastFilterNames]


process.ntuple.TriggerNames = cms.vstring(g_TriggerNames)
process.ntuple.LastFilterNames = cms.vstring(g_LastFilterNames)
process.ntuple.IsMonteCarlo = cms.untracked.bool(True)
process.ntuple.KeepGENOnly  = cms.untracked.bool(False)
