import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'MC_37Y_V5::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/scratch/appelte1/icali_NZS/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_1002_1_cRm.root'
    )
)

from RecoLocalTracker.SiStripZeroSuppression.DefaultAlgorithms_cff import *

process.cmn = cms.EDAnalyzer('SiStripCMNAnalyzer',
    Algorithms = DefaultAlgorithms,
    RawDigiProducersList = cms.InputTag('simSiStripDigis','VirginRaw')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cmn.root')
)



process.p = cms.Path(process.cmn)
