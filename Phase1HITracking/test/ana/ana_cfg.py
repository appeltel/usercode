import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbTrackingCorrections_cfi')

# input variable parsing

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.outputFile = 'ana_test.root'

ivars.inputFiles = 'reco_test.root'

ivars.parseArguments()


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile)
)



process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff")

process.tpRecoAssocGeneralTracks = process.trackingParticleRecoTrackAsssociation.clone()
#process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("hiSelectedTracks")
#process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("generalTracks")
#process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("hiQuadPixelTracks")
#process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("hiGeneralTracks")
process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("hiConformalPixelTracks")
#process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("hiPhase1PixelTracks")

#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
#process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')



# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
#'file:step3P1_0ppTrk.root'
    ivars.inputFiles
    )
)

process.pPbTrkCorr.occByCentrality = False
process.pPbTrkCorr.occByNull = True
process.pPbTrkCorr.applyTrackCuts = False
#process.pPbTrkCorr.applyTrackCuts = True
process.pPbTrkCorr.dxyErrMax = 999.0
process.pPbTrkCorr.dzErrMax = 999.0
process.pPbTrkCorr.ptErrMax = 999.0
#process.pPbTrkCorr.trackSrc = "hiGeneralTracks"
#process.pPbTrkCorr.trackSrc = "hiSelectedTracks"
process.pPbTrkCorr.vertexSrc = "hiSelectedVertex"
#process.pPbTrkCorr.trackSrc = "generalTracks"
process.pPbTrkCorr.trackSrc = "hiConformalPixelTracks"
#process.pPbTrkCorr.trackSrc = "hiPhase1PixelTracks"
#process.pPbTrkCorr.vertexSrc = "offlinePrimaryVertices"
process.pPbTrkCorr.fillTrkPerfHistos = True

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.GlobalTag.globaltag = 'STAR17_61_V1A::All'

process.p = cms.Path( 
                      process.tpRecoAssocGeneralTracks *
                      process.pPbTrkCorr 
)

