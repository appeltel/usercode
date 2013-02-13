import FWCore.ParameterSet.Config as cms

pPbTrkCorr  = cms.EDAnalyzer('RpPbTrackingCorrections',
    trackSrc = cms.InputTag("generalTracks"),
    vertexSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
    tpEffSrc = cms.InputTag('mergedtruth','MergedTrackTruth'),
    tpFakSrc = cms.InputTag('mergedtruth','MergedTrackTruth'),
    associatorMap = cms.InputTag('tpRecoAssocGeneralTracks'),
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
    occBins = cms.vdouble( 0., 100.),
    occByCentrality = cms.bool(True),
    occByNPixelTrk = cms.bool(False),
    applyCuts = cms.bool(True),
    qualityString = cms.string("highPurity"),
    dxyErrMax = cms.double(3.0),
    dzErrMax = cms.double(3.0),
    ptErrMax = cms.double(0.1)
)
