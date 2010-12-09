import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGEDTRACKANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.load("edwenger.HiTrkEffAnalyzer.TrackSelections_cff")

process.load('Appeltel.PixelTracksRun2010.HiLowPtPixelTracksFromReco_cff')
process.load('Appeltel.PixelTracksRun2010.HiMultipleMergedTracks_cff')
process.load('Appeltel.PixelTracksRun2010.HiTrackCandidates_cff')

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.GlobalTag.globaltag = 'START39_V4HI::All'

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #   '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_253_1_kCA.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_252_1_m3O.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_250_1_2Aw.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_24_1_Qt9.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_249_1_htn.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_248_1_z4G.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_247_1_QSk.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_246_1_TJm.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_245_1_BDJ.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_244_1_vr4.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_243_1_cBA.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_242_1_HFQ.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_241_1_nGz.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_240_2_YVg.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_23_1_1xS.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_239_1_Zjo.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_238_1_Zm7.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_237_1_wwJ.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_236_1_L4f.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_235_1_l8z.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_233_1_2As.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_231_2_11l.root',
       '/store/user/davidlw/Hydjet_Quenched_MinBias_2760GeV/Fall10-START39_V4HI-391redigi-v5-GEN-SIM-RECO/597653989fad95507f68df8cea7f0ee4/hiReco_RAW2DIGI_RECO_230_1_kM3.root'
    )
)


process.mergedTrackAna = cms.EDAnalyzer('MergedTrackAnalyzer',
    trackSrc = cms.InputTag('hiGoodMergedTracks'),
    vertexSrc = cms.InputTag('hiSelectedVertex'),
    genSrc = cms.InputTag('hiGenParticles'),
    doGenerator = cms.bool(True),
    etaCut = cms.double(2.4),
    ptMin = cms.double(0.2)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('mergedTracks_eta24.root')
)


process.p = cms.Path(
               process.hiGoodTracksSelection *
               process.conformalPixelTrackReco *
               process.hiGoodMergedTracks * 
               process.mergedTrackAna)

