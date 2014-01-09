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
process.load('Appeltel.RpPbAnalysis.RpPbTrackingCorrections_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

#process.load('Appeltel.RpPbAnalysis.HLT_PIon_partial_cff')
#process.load('Appeltel.RpPbAnalysis.HLT_PIon_partial_cff')
process.load('HLTrigger.Configuration.HLT_PIon_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trackCorrections.root')
)

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string("Hijing"),
  centralitySrc = cms.InputTag("pACentrality"),
  pPbRunFlip = cms.untracked.uint32(99999999)
  )

process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpRecoAssocGeneralTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocGeneralTracks.label_tr = cms.InputTag("generalTracks")

process.tpRecoAssocHLTTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocHLTTracks.label_tr = cms.InputTag("hltPAGoodFullTracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames =  cms.untracked.vstring(
'file:/scratch/appelte1/HIJING_mix_Dijet120_1000ev.root'
    )
)

# JETS

process.load('CmsHi.JetAnalysis.ExtraPfReco_cff')
process.load('CmsHi.JetAnalysis.ExtraJetReco_cff')
process.load('CmsHi.JetAnalysis.PatAna_cff')
process.PFTowers.src = cms.InputTag("particleFlow")

#process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
#process.hltSingleTrigger = process.hltHighLevel.clone()
#process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

process.trkCorr_HIN12017.doMomRes = cms.bool(False)

process.trkCorr_HIN12017.doHLT = cms.bool(True)

#process.trkCorr_HIN12017.trackSrc = 'hltPAGoodFullTracks'
#process.trkCorr_HIN12017.associatorMap = 'tpRecoAssocHLTTracks'

process.trkCorr_HIN12017.fillTrkPerfHistos = cms.bool(True)
process.trkCorr_HIN12017.fillTrkPerfHistosRF = cms.bool(True)

#process.trkCorr_HIN12017.applyTrackCuts = cms.bool(False)

process.trkCorr_HIN12017_dzloose = process.trkCorr_HIN12017.clone(
    reweightDz = cms.bool(True),
    reweightDzFunc = cms.string("1.0-0.4*exp(-(x/3.0)**2)")
)

process.trkCorr_HIN12017_dztight = process.trkCorr_HIN12017.clone(
    reweightDz = cms.bool(True),
    reweightDzFunc = cms.string("1.0+0.4*exp(-(x/3.0)**2)")
)

process.trkCorr_HIN12017_dxyloose = process.trkCorr_HIN12017.clone(
    reweightDxy = cms.bool(True),
    reweightDxyFunc = cms.string("1.0-0.4*exp(-(x/3.0)**2)")
)

process.trkCorr_HIN12017_nhitloose = process.trkCorr_HIN12017.clone(
    reweightNhit = cms.bool(True),
    reweightNhitFunc = cms.string("1.2-0.4*x/30.")
)

process.trkCorr_HIN12017_nhittight = process.trkCorr_HIN12017.clone(
    reweightNhit = cms.bool(True),
    reweightNhitFunc = cms.string("1.0+0.4*x/30.")
)

process.trkCorr_HIN12017_dxytight = process.trkCorr_HIN12017.clone(
    reweightDxy = cms.bool(True),
    reweightDxyFunc = cms.string("1.0+0.4*exp(-(x/3.0)**2)")
)

process.trkCorr_HIN12017_NC = process.trkCorr_HIN12017.clone(
    applyTrackCuts = cms.bool(False)
)

process.trkCorr_HIN12017_pos = process.trkCorr_HIN12017.clone(
    chargeCut = cms.int32(1)
)

process.trkCorr_HIN12017_neg = process.trkCorr_HIN12017.clone(
    chargeCut = cms.int32(-1)
)


process.trkCorr_HIN12017_J2 = process.trkCorr_HIN12017.clone(
    occBins = cms.vdouble( 
        0., 20., 40., 60., 80., 100., 120., 140., 160., 
        180., 200., 220., 240., 260., 280., 300.,
        350., 400., 500.
    )
)

process.trkCorr_HIN12017_MatchForest = process.trkCorr_HIN12017.clone(
    ptBins = cms.vdouble(
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.05, 1.1, 1.15, 1.2,
        1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0,
        20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 
        180.0, 300.0, 500.0
    ),
    etaBins = cms.vdouble( 
        -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
        0.4, 0.8, 1.2, 1.6, 2.0, 2.4
    ),
    occBins = cms.vdouble(
        0.,40., 60., 80.,100., 120., 160., 200., 300., 500., 1000.
    ),
    jetEtaMax = cms.double(2.0),
    applyVertexZCut = cms.bool(True),
    vertexZMax = cms.double(15.),
    jetSrc = cms.InputTag('ak3PFpatJets'),
    doVtxReweighting = cms.bool(False)
)

process.trkCorr_HIN12017_MatchForestPU = process.trkCorr_HIN12017_MatchForest.clone(
    jetSrc = cms.InputTag('akPu3PFpatJets')
)

process.trkCorr_HIN12017_MatchForestRW = process.trkCorr_HIN12017_MatchForest.clone(
    doVtxReweighting = cms.bool(True)
)

process.GlobalTag.globaltag = 'STARTHI53_V27::All'

process.jetReco= cms.Sequence(
    process.PFTowers *
    process.akPu3PFJets*
    process.akPu3PFcorr*
    process.akPu3PFpatJets*
    process.ak3PFJets*
    process.ak3PFcorr*
    process.ak3PFpatJets*
    process.akPu3CaloJets*
    process.akPu3CaloJetID*
    process.akPu3Calocorr*
    process.akPu3CalopatJets
)


# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

#process.hltTrkReco_v1 = cms.Sequence(
#    process.HLTDoLocalPixelSequence *
#    process.HLTRecopixelvertexingForHighMultPASequence *
#    process.HLTDoLocalStripSequence *
#    process.HLTPAIterativeTracking *
#    process.hltPAGoodFullTracks
#)

process.hltTrkReco_v1 = cms.Sequence( 
    process.hltScalersRawToDigi +
    process.hltOnlineBeamSpot +
    process.HLTRecoJetSequencePrePF + 
    process.HLTDoLocalPixelSequence + 
    process.HLTRecopixelvertexingForHighMultPASequence + 
    process.HLTDoLocalStripSequence + 
    process.HLTIterativeTrackingForPA + 
    process.hltPAGoodFullTracks 
)


#process.hltTrkReco_v2 = cms.Sequence(
#    process.HLTDoLocalHfSequence +
#    process.hltHcalPM2Tower3GeVFilter +
#    process.HLTRecoJetSequencePrePF +
#    process.HLTDoLocalPixelSequence +
#    process.HLTRecopixelvertexingForHighMultPASequence +
#    process.HLTDoLocalStripSequence +
#    process.HLTIterativeTrackingForPA +
#    process.hltPAGoodFullTracks
#)

process.p = cms.Path( 
    process.siPixelRecHits *
    process.PAcollisionEventSelection * 
    process.pACentrality *
    process.jetReco *
    process.hltTrkReco_v1 *
    process.tpRecoAssocHLTTracks *
    process.tpRecoAssocGeneralTracks *
    process.trkCorr_HIN12017 
#    process.trkCorr_HIN12017_dxyloose * 
#    process.trkCorr_HIN12017_dxytight * 
#    process.trkCorr_HIN12017_dzloose * 
#    process.trkCorr_HIN12017_dztight * 
#    process.trkCorr_HIN12017_nhitloose * 
#    process.trkCorr_HIN12017_nhittight  
#    process.trkCorr_HIN12017_J2 *
#    process.trkCorr_HIN12017_MatchForest *
#    process.trkCorr_HIN12017_MatchForestPU *
#    process.trkCorr_HIN12017_MatchForestRW
#    process.trkCorr_HIN12017_loose *
#    process.trkCorr_HIN12017_loose2 *
#    process.trkCorr_HIN12017_loose3 *
#    process.trkCorr_HIN12017_tight *
#    process.trkCorr_HIN12017_tight2 *
#    process.trkCorr_HIN12017_tight3 *
#    process.trkCorr_HIN12017_Z25 *
#    process.trkCorr_HIN12017_NoRW *
#    process.trkCorr_HIN12017_wideJetBin *
#    process.trkCorr_HIN12017_wideEtaBin *
#    process.trkCorr_HIN12017_Pbp


)

process.schedule = cms.Schedule( process.p )
