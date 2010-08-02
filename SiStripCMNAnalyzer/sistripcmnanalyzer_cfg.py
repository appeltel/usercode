import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.RandomNumberGeneratorService.cmn = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
)

process.GlobalTag.globaltag = 'MC_37Y_V5::All'

process.poolDBESSource = cms.ESSource("PoolDBESSource",
   BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
   DBParameters = cms.PSet(
       messageLevel = cms.untracked.int32(2),
       authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
   ),
   timetype = cms.untracked.string('runnumber'),
   connect = cms.string('frontier://FrontierProd/CMS_COND_31X_STRIP'),
   toGet = cms.VPSet(
       cms.PSet(
           record = cms.string('SiStripNoisesRcd'),
           tag = cms.string('SiStripNoise_CRAFT09_DecMode_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripPedestalsRcd'),
           tag = cms.string('SiStripPedestal_CRAFT09_DecMode_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripFedCablingRcd'),
           tag = cms.string('SiStripFedCabling_CRAFT09_ForTrackerSim')
       ),
       cms.PSet(
           record = cms.string('SiStripBadChannelRcd'),
           tag = cms.string('SiStripBadChannelsFromO2O_CRAFT09_DecMode_ForTrackerSim')
       )
   )
)

process.es_prefer_my =cms.ESPrefer("PoolDBESSource","poolDBESSource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_9_1_Uea.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_99_1_iVy.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_98_1_IW6.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_97_1_jsZ.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_96_1_DkZ.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_95_1_VGg.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_94_1_ZrQ.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_93_1_wDk.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_723_1_3ls.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_722_1_Uvb.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_695_1_kxa.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_694_1_dDK.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_693_1_Gxh.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_692_1_P1f.root',
       '/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_OnlyBS_FiveTimesBigger/35475ada32de29e04cd11322eaedaaea/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_OnlyBS_691_1_mAg.root'
    )
)

from RecoLocalTracker.SiStripZeroSuppression.DefaultAlgorithms_cff import *
from RecoLocalTracker.SiStripClusterizer.DefaultClusterizer_cff import *

process.cmn = cms.EDAnalyzer('SiStripCMNAnalyzer',
    Clusterizer = DefaultClusterizer,
    Algorithms = DefaultAlgorithms,
    RawDigiProducersList = cms.InputTag('simSiStripDigis','VirginRaw'),
    RawDigiProducersListNoise = cms.InputTag('simSiStripDigis','VirginRawSignal'),
    RawDigiProducersListSignal = cms.InputTag('simSiStripDigis','VirginRawNoise'),
    doNoiseAndSignal = cms.bool(False)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('cmn.root')
)



process.p = cms.Path(process.cmn)
