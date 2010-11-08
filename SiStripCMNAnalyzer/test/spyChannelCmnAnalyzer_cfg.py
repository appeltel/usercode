import FWCore.ParameterSet.Config as cms

process = cms.Process("SPY")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Services_cff")

## CMN analyzer takes a random seed
process.RandomNumberGeneratorService.cmn = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V9::All'

process.Timing = cms.Service("Timing")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    ## Eric's unpacked SpyChannel files with VR digis
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_1.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_10.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_11.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_12.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_13.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_14.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_15.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_16.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_17.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_18.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_19.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_2.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_20.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_21.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_22.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_23.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_24.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_25.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_26.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_27.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_28.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_29.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_30.root',
    'rfio:/castor/cern.ch/user/e/ebutz/SpyDigis/142422/SiStripSpyDigis_142422_31.root',
    )
)

from RecoLocalTracker.SiStripZeroSuppression.DefaultAlgorithms_cff import *
from RecoLocalTracker.SiStripClusterizer.DefaultClusterizer_cff import *

process.cmn = cms.EDAnalyzer('SiStripCMNAnalyzer',
    Clusterizer = DefaultClusterizer,
    Algorithms = DefaultAlgorithms,
    RawDigiProducersList = cms.InputTag('SiStripSpyDigiConverter','VirginRaw'), ## SpyDigi
    RawDigiProducersListNoise = cms.InputTag('SiStripSpyDigiConverter','VirginRawSignal'),
    RawDigiProducersListSignal = cms.InputTag('SiStripSpyDigiConverter','VirginRawNoise'),
    doNoiseAndSignal = cms.bool(False),
    galleryClusterMin = cms.int32(1),
    galleryClusterMin = cms.int32( 2 ),
    nSigmaNoiseDerTh =  cms.uint32(4),
    consecThreshold = cms.uint32(5),
    hitStripThreshold = cms.uint32(40),    
    nSmooth = cms.uint32(9),      
    minStripsToFit = cms.uint32(4)                         
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cmn.root')
)


process.p = cms.Path(process.cmn)
