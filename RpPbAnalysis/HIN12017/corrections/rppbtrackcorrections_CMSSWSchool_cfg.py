import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('Appeltel.RpPbAnalysis.RpPbSimpleTrackingCorrections_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackCorrections.root')
)


process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpRecoAssocGeneralTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("generalTracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
'/store/user/vzhukova/HIJING_GEN-SIM_YUE-SHI_Minbias_2_v1/HIJING_RECO_YUE-SHI_Minbias__2_v1/b7d33bba7673cdb1ee6f4983c0800c79/hijing_reco_fix_100_1_Y3O.root' 
   )
)


process.GlobalTag.globaltag = 'STARTHI53_V25::All'

process.p = cms.Path( 
    process.tpRecoAssocGeneralTracks *
    process.trkCorr_HIN12017 
)

