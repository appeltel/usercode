import FWCore.ParameterSet.Config as cms

import Appeltel.RpPbAnalysis.RpPbTrackingAnalyzer_cfi


trkAna_quality = trkAna.clone(
        applyCuts = cms.bool(True)
)

trkAna_pixel = trkAna.clone(
   trackSrc = cms.InputTag("pixelTracks")
)

trkAna_noBS = trkAna.clone(
   vertexSrc = cms.InputTag("offlinePrimaryVertices")
)

trkAnaMinBias = cms.Sequence( trkAna 
                              * trkAna_quality
                              * trkAna_pixel
                              * trkAna_noBS
)


trkAna_HIN12017 = trkAna.clone(
   applyCuts = cms.bool(True),
   qualityString = cms.string('highPurity'),
   dxyErrMax = cms.double(3.0),
   dzErrMax = cms.double(3.0),
   ptErrMax = cms.double(0.1),
   vertexSrc = cms.InputTag("offlinePrimaryVertices"),
   ptBins = cms.vdouble(
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4,
        7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6,
        130.0
   ),
   etaBins = cms.vdouble( 
        -2.46, -2.36, -2.26, -2.16, -2.06, -1.96, -1.86, -1.76, -1.66, -1.56,
        -1.46, -1.36, -1.26, -1.16, -1.06, -0.96, -0.86, -0.76, -0.66, -0.56,
        -0.46, -0.36, -0.26, -0.16, -0.06, 0.04, 0.14, 0.24, 0.34, 0.44, 0.54, 
        0.64, 0.74, 0.94, 0.94, 1.04, 1.14, 1.24, 1.34, 1.44, 1.54, 1.64, 1.74,
        1.84, 1.94, 2.04, 2.14, 2.24, 2.34, 2.44
   ),
   occByLeadingJetEt = cms.bool(True),
   jetEtaMax = cms.double(2.5),
   occBins = cms.vdouble( 
        0., 20., 40., 60., 80., 100., 120., 140., 160., 
        180., 200., 250., 300., 500.
   ),
   vertexZMax = cms.double(15.),
   doTrigEffCorrection = cms.bool(True),
   trigEffByMult = cms.vdouble( 
       0.0, 0.056, 0.483, 0.809, 0.932, 0.966, 0.973, 
       0.981, 0.985, 0.991, 0.993, 0.994, 0.995 
   ),
   zeroMultFraction = cms.double( 0.00083 )                        
)

trkAna_HIN12017_MC = trkAna_HIN12017.clone(
    doMC = cms.bool(True),
    genSrc = cms.InputTag("genParticles")
)

trkAna_HIN12017_TP = trkAna_HIN12017.clone(
    doMC = cms.bool(True),
    genSrc = cms.InputTag("genParticles"),
    doMCbyTP = cms.bool(True)
)
