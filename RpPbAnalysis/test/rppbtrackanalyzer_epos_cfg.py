import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackAnalysis.root')
)

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
# private reco at MIT
'/store/user/kjung/VZhukova/EposBatch13/pPbRecoShift_f1521_e.root'
)

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'STARTHI53_V17::All'

process.trkAna_etaFull = cms.EDAnalyzer('RpPbTrackingAnalyzer',
   trackSrc = cms.InputTag("generalTracks"),
   vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
   etaMin = cms.double(-2.5),
   etaMax = cms.double(2.5)
)

process.trkAna_etaN2 = process.trkAna_etaFull.clone(
   etaMin = cms.double(-2.5),
   etaMax = cms.double(-1.5)
)

process.trkAna_etaN1 = process.trkAna_etaFull.clone(
   etaMin = cms.double(-1.5),
   etaMax = cms.double(-0.5)
)

process.trkAna_eta0 = process.trkAna_etaFull.clone(
   etaMin = cms.double(-0.5),
   etaMax = cms.double(0.5)
)

process.trkAna_etaP1 = process.trkAna_etaFull.clone(
   etaMin = cms.double(0.5),
   etaMax = cms.double(1.5)
)

process.trkAna_etaP2 = process.trkAna_etaFull.clone(
   etaMin = cms.double(1.5),
   etaMax = cms.double(2.5)
)



process.p = cms.Path( process.PAcollisionEventSelection * 
                      process.trkAna_etaFull * 
                      process.trkAna_etaN2 * 
                      process.trkAna_etaN1 * 
                      process.trkAna_eta0 * 
                      process.trkAna_etaP1 * 
                      process.trkAna_etaP2
)

