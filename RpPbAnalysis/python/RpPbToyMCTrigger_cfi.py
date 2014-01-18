import FWCore.ParameterSet.Config as cms

trkAna  = cms.EDAnalyzer('RpPbToyMCTrigger',
   ptBins = cms.vdouble(
        0.0, 0.1, 0.2, 0.3, 0.4, 
        0.5, 0.6, 0.7, 0.8, 0.9, 
        1.1, 1.2, 1.4, 1.6, 1.8, 
        2.0, 2.2, 2.4, 3.2, 4.0, 
        4.8, 5.6, 6.4, 7.2, 9.6, 
        12.0, 14.4, 19.2, 24.0, 28.8, 
        35.2, 41.6, 48.0, 60.8, 73.6, 
        86.4, 103.6, 130.0
   ),
   etaMin = cms.double(-1.465),
   etaMax = cms.double(0.535),
   etaShift = cms.double(-0.465),
   genSrc = cms.InputTag("genParticles"),
   mbPrescale = cms.double(1000.0),
   trk12Prescale = cms.double(20.0),
   trk20Prescale = cms.double(5.0),
   correctionFile = cms.string("trackCorrections_HIN12017HLTv2_XSecWeighted.root"),
   correctionDir = cms.string("trkCorr_HIN12017")
)
