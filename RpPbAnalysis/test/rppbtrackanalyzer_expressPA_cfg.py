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

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/express/HIRun2013A/ExpressPhysics/FEVT/Express-v1/000/209/889/E6E79F5D-8F5C-E211-A086-BCAEC518FF80.root'
  )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'GR_E_V33::All'

process.p = cms.Path( process.hltSingleTrigger * 
                      process.PAcollisionEventSelection * 
                      process.trkAnaMinBias
)
