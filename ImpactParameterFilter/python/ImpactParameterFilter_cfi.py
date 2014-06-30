import FWCore.ParameterSet.Config as cms

impactParameterFilter  = cms.EDFilter('ImpactParameterFilter',
    generators = cms.vstring("generator"),
    bMin = cms.double(-1.0),
    bMax = cms.double(9999.0)
)
