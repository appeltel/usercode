import FWCore.ParameterSet.Config as cms

process = cms.Process('CORRFLOWSKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.load('Appeltel.PixelTracksRun2010.HILowPtPixelTracks_cfi')
process.load('Appeltel.PixelTracksRun2010.HILowPtWideCutPixelTracks_cfi')
process.load('Appeltel.PixelTracksRun2010.HIInvertedLowPtPixelTracks_cfi')
process.load('Appeltel.PixelTracksRun2010.HICorrFlowSkimEventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'file:/scratch/appelte1/150431-test/948A312C-3EEB-DF11-8811-0030487CAF0E.root'
)


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.HICorrFlowEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('filter_step')
                ),
    fileName = cms.untracked.string( 'file:/scratch/appelte1/skim_data_test.root')
)


#Trigger Selection

import HLTrigger.HLTfilters.hltHighLevel_cfi
hltMinBiasHF = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltMinBiasHF.HLTPaths = ["HLT_HIMinBiasHF"]


# Charged Candidates

process.allPixelTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiLowPtPixelTracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allInvertedLowPtPixelTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiInvertedLowPtPixelTracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allSelectedTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiSelectedTracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allMergedNoPtSplitTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedNoPtSplitTracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allMergedPtSplit09Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit09Tracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allMergedPtSplit12Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit12Tracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allMergedPtSplit15Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit15Tracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)

process.allMergedPtSplit18Tracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiMergedPtSplit18Tracks"), # or whatever you call the more tightly selected collection
    particleType = cms.string('pi+')
)


import Appeltel.PixelTracksRun2010.hiTrackListMerger_cfi

process.hiMergedNoPtSplitTracks = Appeltel.PixelTracksRun2010.hiTrackListMerger_cfi.hiTrackListMerger.clone(
    TrackProducer1 = 'hiSelectedTracks',
    TrackProducer2 = 'hiLowPtPixelTracks',
    promoteTrackQuality = False,  
    allowFirstHitShare = False   
)

process.hiMergedPtSplit09Tracks = process.hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 0.9,
    trackCollection1MinPt = 0.9
)

process.hiMergedPtSplit12Tracks = process.hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.2,
    trackCollection1MinPt = 1.2
)

process.hiMergedPtSplit15Tracks = process.hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.5,
    trackCollection1MinPt = 1.5
)

process.hiMergedPtSplit18Tracks = process.hiMergedNoPtSplitTracks.clone(
    trackCollection2MaxPt = 1.8,
    trackCollection1MinPt = 1.8
)


# Other statements
process.GlobalTag.globaltag = 'GR10_P_V12::All'



# Path and EndPath definitions

process.filter_step = cms.Path( process.hltMinBiasHF )

process.lowptpixel_step = cms.Path(
    process.siPixelRecHits * 
    process.hiLowPtPixelTracks *
    process.hiLowPtWideCutPixelTracks *
    process.hiInvertedLowPtPixelTracks 
)
process.merge_step = cms.Path(
    process.hiMergedNoPtSplitTracks *
    process.hiMergedPtSplit09Tracks *
    process.hiMergedPtSplit12Tracks *
    process.hiMergedPtSplit15Tracks *
    process.hiMergedPtSplit18Tracks
)
process.candidates_step = cms.Path(
    process.allPixelTracks *
    process.allInvertedLowPtPixelTracks *
    process.allSelectedTracks *
    process.allMergedNoPtSplitTracks *
    process.allMergedPtSplit09Tracks *
    process.allMergedPtSplit12Tracks *
    process.allMergedPtSplit15Tracks *
    process.allMergedPtSplit18Tracks 
)

process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(
    process.filter_step,
    process.lowptpixel_step,
    process.merge_step,
    process.candidates_step,
    process.endjob_step,
    process.out_step
)

