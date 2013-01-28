import FWCore.ParameterSet.Config as cms

import Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cfi

pileupVertexFilterCutG = Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cfi.pileupVertexFilter.clone()

pileupVertexFilterCutGloose = pileupVertexFilterCutG.clone(
    dzCutByNtrk = cms.vdouble(
        999., 4.5, 3.2, 3.0, 1.8, 1.8, 1.35, 0.9
    )
)

pileupVertexFilterCutGtight = pileupVertexFilterCutG.clone(
    dzCutByNtrk = cms.vdouble(
        999., 2.0, 1.6, 1.333, 0.8, 0.8, 0.6, 0.4
    )
)

pileupVertexFilterCutE = pileupVertexFilterCutG.clone(
    doDzNtrkCut = cms.bool(False),
    doDxyDzCut = cms.bool(True)
)

pileupVertexFilterCutEandG = pileupVertexFilterCutG.clone(
    doDzNtrkCut = cms.bool(True),
    doDxyDzCut = cms.bool(True)
)

pileupVertexFilterDevelBits = cms.Sequence( 
    pileupVertexFilterCutG +
    pileupVertexFilterCutGloose +
    pileupVertexFilterCutGtight +
    pileupVertexFilterCutE + 
    pileupVertexFilterCutEandG
)
