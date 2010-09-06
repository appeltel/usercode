import FWCore.ParameterSet.Config as cms

process = cms.Process('NEWCMN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('SimGeneral.MixingModule.himixDIGI_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.ValidationHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(1)
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(

)

from RecoLocalTracker.SiStripZeroSuppression.customiseCMN import *

customiseAPVRestore(process)

customiseMedian(process)
#customisePercentile(process)
#customiseIteratedMedian(process)


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


process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
oncePerEventMode = cms.untracked.bool(False),
ignoreTotal = cms.untracked.int32(0)
)


# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#'dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/vzhukova/Hydjet_Quenched_MinBias_2760GeV/Summer10_b0_DiJet_2760GeV_MC_37Y_V4-prelim/a8253289b6f469855700395662aa23dd/DiJet_2760GeV_MC_37Y_V4_9_1_eYj.root'
#'/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_OnlySSTRedigi_VR_FullStatistic/585635c9d4b54265067465e6a8a60b81/DiJet_B0_2760GeV_MC_37Y_V5_SSTREDIGI_TrkVR_99_1_cls.root'
'/store/user/icali/Hydjet_Quenched_MinBias_2760GeV/Hydjet_Quenched_b0_DiJet120_2760GeV_MC_37Y_V5_TrkVR_Full/8edde60506a6e07603d62be0000b77d0/DiJet_B0_2760GeV_MC_37Y_V5_RECO_TrkVR_984_1_wd6.root'
                                   ),
   inputCommands = cms.untracked.vstring(
       'keep *',
       'drop *_hiGlobalPrimTracks_*_*',
       'drop *_hiSelectedTracks_*_*',
       'drop *_siPixelClusters_*_*',
       'drop *_hiPixel3PrimTracks_*_*',
       'drop *_hiPixel3ProtoTracks_*_*',
       'drop *_hiSelectedProtoTracks_*_*'
   )

)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
   splitLevel = cms.untracked.int32(0),
   #outputCommands = process.FEVTDEBUGEventContent.outputCommands,
   outputCommands = cms.untracked.vstring( 'drop *_*_*_*', 'keep *_MEtoEDMConverter_*_*' ),
   fileName = cms.untracked.string('file:validation.root'
                             ),
   dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO'),
        filterName = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'MC_37Y_V5::All'


process.validationHeavyIons.remove(process.hiJetValidation)


# Path and EndPath definitions
#process.mix_step = cms.Path(process.mix)
#process.localreco_step = cms.Path( process.newlocalreco)
process.raw2digi_step = cms.Path( process.RawToDigi )
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)
process.validation_step = cms.Path(process.hiTrackValidation)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)




# Schedule definition
process.schedule = cms.Schedule(
   process.raw2digi_step,
   process.reconstruction_step,
   process.validation_step,
   process.endjob_step,
   process.out_step
)
