import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackAnalysis.root')
)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
#  pPbRunFlip = cms.uint32(99999999)
  )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/express/HIRun2013/ExpressPhysics/FEVT/Express-v1/000/210/412/00000/3A1A368E-C461-E211-B719-001D09F2B2CF.root'
  )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]
#process.hltSingleTrigger.HLTPaths = ["HLT_PAMinBiasHfOrBHC_v1"]

process.options = cms.untracked.PSet(
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.GlobalTag.globaltag = 'GR_E_V33::All'

process.p = cms.Path( process.hltSingleTrigger * 
                      process.PAcollisionEventSelection *
#                      process.siPixelRecHits *
                      process.pACentrality * 
                      process.trkAnaMinBias
)
