# Auto generated configuration file
# using: 
# Revision: 1.182 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: hiRecoJEX -n 2 -s RAW2DIGI,RECO --scenario HeavyIons --geometry DB --conditions auto:mc --datatier RECO --eventcontent=FEVTDEBUG --filein=file.root --no_exec
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

#process.Timing = cms.Service("Timing")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
)


################################################################################
# Input Variable Parsing
#

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')


ivars.inputFiles = [ 'file:/scratch/appelte1/store/results/heavy-ions/MinimumBiasHI/StoreResults-Spring10-JulyAnalysisExercise_MC_37Y_V4-RECO-prelim-v1-599a3a4efd1a163b86e329cbf44e9f6a/MinimumBiasHI/USER/StoreResults-Spring10-JulyAnalysisExercise_MC_37Y_V4-RECO-prelim-v1-599a3a4efd1a163b86e329cbf44e9f6a/0000/3051F0E6-448B-DF11-A753-0026B939DCF3.root' ]
ivars.outputFile = 'reco_test.root'

ivars.parseArguments()


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( ivars.inputFiles )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring('drop *'),
    fileName = cms.untracked.string( ivars.outputFile )
)

process.output.outputCommands += [ 'keep recoTracks_hiSelectedTracks_*_*' ]
process.output.outputCommands += [ 'keep recoTracks_hiLowPtPixelTracks_*_*' ]
process.output.outputCommands += [ 'keep *_hiSelectedVertex_*_*' ]
process.output.outputCommands += [ 'keep *_hiCentrality_*_*' ]
process.output.outputCommands += [ 'keep *_hiEvtPlane_*_*' ]
process.output.outputCommands += [ 'keep recoCaloJets_*_*_*' ]
process.output.outputCommands += [ 'keep *_allPixelTracks_*_*' ]
process.output.outputCommands += [ 'keep *_allSelectedTracks_*_*' ]
process.output.outputCommands += [ 'keep *_allMergedNoPtSplitTracks_*_*' ]
process.output.outputCommands += [ 'keep *_allMergedPtSplit09Tracks_*_*' ]
process.output.outputCommands += [ 'keep *_allMergedPtSplit12Tracks_*_*' ]
process.output.outputCommands += [ 'keep *_allMergedPtSplit15Tracks_*_*' ]
process.output.outputCommands += [ 'keep *_allMergedPtSplit18Tracks_*_*' ]

process.hiLowPtPixelTracks = process.hiPixel3PrimTracks.clone()
process.hiLowPtPixelTracks.RegionFactoryPSet.RegionPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.ptMin = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaTipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.tipMax = cms.double(0.2)
process.hiLowPtPixelTracks.FilterPSet.nSigmaLipMaxTolerance = cms.double(4.0)
process.hiLowPtPixelTracks.FilterPSet.lipmax = cms.double(0.2)

process.hiLowPtPixelTracks.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = 5000000

process.allPixelTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hiLowPtPixelTracks"), # or whatever you call the more tightly selected collection
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


import Appeltel.JEX.hiTrackListMerger_cfi

process.hiMergedNoPtSplitTracks = Appeltel.JEX.hiTrackListMerger_cfi.hiTrackListMerger.clone(
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

# Additional output definition
#process.load("edwenger.JulyExercise.BeamSpotJEX_cff")

#process.BeamSpotDBSource.connect = cms.string('sqlite_file:PreliminaryBeamSpot.db')

# Other statements
process.GlobalTag.globaltag = 'MC_37Y_V4::All'


########################### tunes for Photon ##############
process.photonCore.minSCEt = cms.double(5.0)
process.photons.minSCEtBarrel = cms.double(5.0)
process.photons.minSCEtEndcap = cms.double(5.0)
process.photons.minR9Barrel = cms.double(0.01)  #0.94
process.photons.minR9Endcap = cms.double(0.01)   #0.95
process.photons.maxHoverEEndcap = cms.double(0.99)  #0.5
process.photons.maxHoverEBarrel = cms.double(0.99)  #0.53
###########################################################      

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.lowptpixel_step = cms.Path(process.siPixelRecHits * process.hiLowPtPixelTracks)
process.merge_step = cms.Path(
    process.hiMergedNoPtSplitTracks *
    process.hiMergedPtSplit09Tracks *
    process.hiMergedPtSplit12Tracks *
    process.hiMergedPtSplit15Tracks *
    process.hiMergedPtSplit18Tracks
)
process.candidates_step = cms.Path(
    process.allPixelTracks *
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
#    process.raw2digi_step,
#    process.reconstruction_step,
    process.lowptpixel_step,
    process.merge_step,
    process.candidates_step,
    process.endjob_step,
    process.out_step
)

