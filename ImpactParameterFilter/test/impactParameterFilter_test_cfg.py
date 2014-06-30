import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Appeltel.ImpactParameterFilter.ImpactParameterFilter_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
'/store/himc/HiFall13/Hydjet1p8_TuneDrum_Quenched_MinBias_2760GeV/GEN-SIM/STARTHI53_V28-v2/00000/02B3F103-995E-E311-BD92-003048C9AC48.root'
  )
)

process.ipFilter = process.impactParameterFilter.clone(
  bMin = -1.0,
  bMax = 4.0
)


#process.GlobalTag.globaltag = 'GR_R_53_LV6::All'


process.p = cms.Path( process.ipFilter)
