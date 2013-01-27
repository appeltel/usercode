import FWCore.ParameterSet.Config as cms

pileupVertexFilter  = cms.EDFilter('PAPileUpVertexFiliter',
    vtxSrc = cms.InputTag("offlinePrimaryVertices"),
    dxyCut = cms.double(0.02),
    trkCut = cms.int32(20),
    dzCutByNtrk = cms.vdouble(
        999., 999., 2.5,  2.0,
        1.5,  1.2,  1.0,  0.8,
        0.5,  0.4,  0.35, 0.3,
        0.25, 0.25, 0.2,  0.2,
        0.2,  0.15, 0.15, 0.1
    )
)
