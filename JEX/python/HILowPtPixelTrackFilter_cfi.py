import FWCore.ParameterSet.Config as cms

HiInvertedLowPtFilterBlock = cms.PSet(
    ComponentName = cms.string( "HIInvertedPixelTrackFilter" ),
    ptMin = cms.double( 0.2 ),
    chi2 = cms.double( 0.0 ),
    useClusterShape = cms.bool( False ),
    VertexCollection = cms.InputTag("hiSelectedVertex"),
    nSigmaTipMinTolerance = cms.double( 5.0 ),
    tipMin = cms.double( 0.25 ),
    nSigmaLipMinTolerance = cms.double( 5.0 ),
    lipMin = cms.double( 0.25 )
    )

HiLowPtFilterBlock = cms.PSet(
    ComponentName = cms.string( "HIPixelTrackFilter" ),
    ptMin = cms.double( 0.2 ),
    chi2 = cms.double( 1000.0 ),
    useClusterShape = cms.bool( False ),
    VertexCollection = cms.InputTag("hiSelectedVertex"),
    nSigmaTipMaxTolerance = cms.double( 4.0 ),
    tipMax = cms.double( 0.2 ),
    nSigmaLipMaxTolerance = cms.double( 4.0 ),
    lipMax = cms.double( 0.2 )
    )

HiLowPtWideCutFilterBlock = cms.PSet(
    ComponentName = cms.string( "HIPixelTrackFilter" ),
    ptMin = cms.double( 0.2 ),
    chi2 = cms.double( 1000.0 ),
    useClusterShape = cms.bool( False ),
    VertexCollection = cms.InputTag("hiSelectedVertex"),
    nSigmaTipMaxTolerance = cms.double( 6.0 ),
    tipMax = cms.double( 0.3 ),
    nSigmaLipMaxTolerance = cms.double( 6.0 ),
    lipMax = cms.double( 0.3 )
    )
