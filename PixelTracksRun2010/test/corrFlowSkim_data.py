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

process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.load('Appeltel.PixelTracksRun2010.HiTrackCandidates_cff')
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
process.hltMinBiasHF = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMinBiasHF.HLTPaths = ["HLT_HIMinBiasHF"]

# Other statements
process.GlobalTag.globaltag = 'GR10_P_V12::All'

# Path and EndPath definitions

process.filter_step = cms.Path( process.hltMinBiasHF )

process.lowptpixel_step = cms.Path(process.lowPtPixelTrackReco)
process.merge_step = cms.Path(process.multipleMergedTracks)
process.candidates_step = cms.Path( process.hiAllTrackCandidates)
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

