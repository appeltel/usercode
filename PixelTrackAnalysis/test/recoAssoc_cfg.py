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

process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.hiPostGlobalPrimTracks = process.hiGoodTracks.clone(src = cms.InputTag("hiGlobalPrimTracks"),
                                            qualityBit = cms.string(''),
                                            min_relpterr = cms.double(0.06),
                                            min_nhits = cms.uint32(12),
                                            chi2n_par = cms.double(0.4),
                                            d0_par2 =cms.vdouble(999.0, 0.0),
                                            dz_par2 = cms.vdouble(999.0, 0.0))

process.hiGoodTightTracks = process.hiGoodTracks.clone(src = cms.InputTag("hiPostGlobalPrimTracks"),
                                       qualityBit = cms.string('highPurity'),
                                       min_nhits = cms.uint32(13),
                                       chi2n_par = cms.double(0.15))

process.hiGoodTightMergedTracks = process.hiGoodMergedTracks.clone(
        TrackProducer1 = 'hiGoodTightTracks'
)

################################################################################
# Input Variable Parsing
#

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.inputFiles = '/store/himc/Fall10/AMPT_Organ_MinBias_2760GeV/GEN-SIM-RAWDEBUG/START39_V5HI-v1/0001/F45C0034-BCF5-DF11-A36B-00163691D99E.root'

ivars.outputFile = 'pixel_assoc_test.root'

ivars.parseArguments()


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( ivars.inputFiles )
)


# Other statements
process.GlobalTag.globaltag = 'START39_V5HI::All'


process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ivars.outputFile)
)



# Changes to full tracking
#from RecoPixelVertexing.PixelTrackFitting.PixelFitterByConformalMappingAndLine_cfi import *
#process.hiPixel3PrimTracks.FitterPSet = cms.PSet( PixelFitterByConformalMappingAndLine )
#process.hiPixel3PrimTracks.FitterPSet.fixImpactParameter = cms.double(0.0)
#process.hiPixel3PrimTracks.FilterPSet.nSigmaLipMaxTolerance = 40.0
#process.hiPixel3PrimTracks.FilterPSet.chi2 = 500.0
#process.CkfTrajectoryBuilder.maxCand = 20
#process.hiPixel3PrimTracks.FitterPSet.TTRHBuilder = 'TTRHBuilderWithoutAngle4PixelTriplets'  

#process.Timing = cms.Service("Timing")

#####################################################
#
# TrackAssociation
#
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpRecoAssocConformalPixelTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocConformalPixelTracks.label_tr = cms.InputTag("hiConformalPixelTracks")

process.tpRecoAssocGoodMergedTracks = process.trackingParticleRecoTrackAsssociation.clone()
process.tpRecoAssocGoodMergedTracks.label_tr = cms.InputTag("hiGoodTightMergedTracks")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.hitrkPixelEffAnalyzer = cms.EDAnalyzer('HiPixelTrkEffAnalyzer',
    tracks = cms.untracked.InputTag('hiGoodTightMergedTracks'),
    label_tp_effic = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
    label_tp_fake = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),                 
    associatormap = cms.untracked.InputTag('tpRecoAssocGoodMergedTracks'),
    vertices = cms.untracked.InputTag('hiSelectedVertex',''),
    beamspot = cms.untracked.InputTag('offlineBeamSpot'),                         
    fillHistograms = cms.bool(True),
    fillNtuples = cms.bool(False),
    constPtBins = cms.bool(False),
    lowPtMode = cms.bool(False),
    flow2010ShengquanMode = cms.bool(False),
    flow2010Mode = cms.bool(True),
    doAssociation = cms.untracked.bool(False),
    doFlow = cms.bool(True),
    generators = cms.vstring("generator"),
    hasSimInfo = cms.untracked.bool(True),
    doRP = cms.bool(True),
    genParticles = cms.untracked.InputTag('hiGenParticles')
)

process.prePixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    src = cms.InputTag("siPixelClusters"),
    # untracked string ClusterCollLabel   = "siPixelClusters"
    CPE = cms.string('PixelCPEGeneric'),
    VerboseLevel = cms.untracked.int32(0),
    #TanLorentzAnglePerTesla = cms.double(0.106),
    #Alpha2Order = cms.bool(True),
    #speed = cms.int32(0)
)


process.siPixelRecHits = cms.EDProducer('PixelFiducialRemover',
    src = cms.InputTag("prePixelRecHits"),
    sides = cms.vint32( 4 ),
    disks = cms.vint32( 2 ),
    blades = cms.vint32( 4,5,6 )
)

# vertex cut

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("hiSelectedVertex"),
    cut = cms.string("abs(z) <= 10"), 
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#



# kinematic  track cuts

#process.hiGoodTracks.chi2n_par = 20
#process.hiGoodTracks.dz_par2 = cms.vdouble( 2., 0.0)
#process.hiGoodTracks.d0_par2 = cms.vdouble(2., 0.0)
#process.hiGoodTracks.min_relpterr = cms.double(0.03)
#process.hiGoodTracks.min_nhits = cms.uint32(14)

process.hiConformalPixelTracks.FilterPSet.chi2 = 120.
process.hiConformalPixelTracks.FilterPSet.nSigmaLipMaxTolerance = 18.0


# process.hiConformalPixelTracks.FilterPSet.useClusterShape = cms.bool( False )

# Path and EndPath definitions

process.raw2digi_step = cms.Path(process.RawToDigi)
#process.reconstruction_step = cms.Path(process.offlineBeamSpot*process.siPixelClusters*process.prePixelRecHits*process.muonReco*process.heavyIonTracking)
process.reconstruction_step = cms.Path(process.siPixelClusters*process.prePixelRecHits*process.reconstructionHeavyIons)
process.rechit_step = cms.Path(process.primaryVertexFilter*process.prePixelRecHits*process.siPixelRecHits)
process.fulltrack_step = cms.Path( process.primaryVertexFilter*process.heavyIonTracking )
process.trksel_step  = cms.Path(process.primaryVertexFilter*process.hiPostGlobalPrimTracks*process.hiGoodTightTracks)
process.lowptpixel_step = cms.Path(process.primaryVertexFilter* process.hiConformalPixelTracks)
process.merge_step = cms.Path(process.primaryVertexFilter* process.hiGoodTightMergedTracks)
#process.candidates_step = cms.Path( process.hiAllTrackCandidates)
process.association_step = cms.Path(process.primaryVertexFilter*process.tpRecoAssocGoodMergedTracks ) 
process.val_step = cms.Path(process.primaryVertexFilter* process.hitrkPixelEffAnalyzer )
#process.endjob_step = cms.Path(process.endOfProcess)
#process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.reconstruction_step,
    process.rechit_step,
    process.fulltrack_step,
    process.trksel_step,
    process.lowptpixel_step,
    process.merge_step,
#    process.candidates_step,
    process.association_step,
    process.val_step
#    process.endjob_step,
#    process.out_step
)

