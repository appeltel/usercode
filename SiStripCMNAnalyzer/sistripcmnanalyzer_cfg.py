import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.cmn = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
)

process.GlobalTag.globaltag = 'MC_37Y_V5::All'

process.poolDBESSource = cms.ESSource("PoolDBESSource",
   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
   DBParameters = cms.PSet(
       messageLevel = cms.untracked.int32(2),
       authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
   ),
   timetype = cms.untracked.string('runnumber'),
   connect = cms.string('frontier://FrontierProd/CMS_COND_31X_STRIP'),
   toGet = cms.VPSet(
       cms.PSet(
           record = cms.string('SiStripNoisesRcd'),
           tag = cms.string('SiStripNoise_CRAFT09_DecMode_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripPedestalsRcd'),
           tag = cms.string('SiStripPedestal_CRAFT09_DecMode_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripFedCablingRcd'),
           tag = cms.string('SiStripFedCabling_CRAFT09_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripBadChannelRcd'),
           tag = cms.string('SiStripBadChannelsFromO2O_CRAFT09_DecMode_ForTrackerSim')
       )
   )
)

process.es_prefer_my =cms.ESPrefer("PoolDBESSource","poolDBESSource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/net/hibat0007/d00/scratch/appelte1/DiJet_B0_2760GeV_MC_37Y_V4_RECO_TrkZS.root'
    )
)

from RecoLocalTracker.SiStripZeroSuppression.DefaultAlgorithms_cff import *

process.cmn = cms.EDAnalyzer('SiStripCMNAnalyzer',
    Algorithms = DefaultAlgorithms,
    RawDigiProducersList = cms.InputTag('simSiStripDigis','VirginRaw'),
    RawDigiProducersListNoise = cms.InputTag('simSiStripDigis','VirginRawNoise'),
    RawDigiProducersListSignal = cms.InputTag('simSiStripDigis','VirginRawSignal'),
    doNoiseAndSignal = cms.bool(True)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cmn_1sig.root')
)



process.p = cms.Path(process.cmn)
