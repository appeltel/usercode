import FWCore.ParameterSet.Config as cms

process = cms.Process('SIMPLETRACKANA')


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('simpleTracks.root')
)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    #nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2A14D410-44EE-DF11-949B-003048F118DE.root'
  )
)


process.GlobalTag.globaltag = 'GR10_P_V12::All'

process.simpleTrack = cms.EDAnalyzer('SimpleTrackAnalyzer',
   trackSrc = cms.InputTag("hiSelectedTracks"),
   vertexSrc = cms.InputTag("hiSelectedVertex")
)

process.p = cms.Path( process.simpleTrack )
