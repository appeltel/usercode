# Auto generated configuration file
# using: 
# Revision: 1.232.2.6 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: DMPreProcess -n 10 -s RAW2DIGI --eventcontent DATAMIXER --conditions FrontierConditions_GlobalTag,GR10_P_V12::All --data --dbsquery find file,file.parent where run=151353 and lumi=100 and dataset like *HICorePhysics*/RECO --customise SimGeneral.DataMixingModule.customiseMatchRecVertex --cust_function customiseDataMixerOutput --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RAW2DIGI')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.232.2.6 $'),
    annotation = cms.untracked.string('DMPreProcess nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

################################################################################
# Input Variable Parsing
#

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.inputFiles = 'file:///scratch/appelte1/HIRAWSelected/selectedRAW_1.root'

ivars.outputFile = 'preMix.root'

ivars.parseArguments()


# Input source
process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(ivars.inputFiles)
)

process.options = cms.untracked.PSet(

)

# Output definition

process.DATAMIXERoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.DATAMIXEREventContent.outputCommands,
    fileName = cms.untracked.string(ivars.outputFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'GR10_P_V12::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)

process.endjob_step = cms.Path(process.endOfProcess)

process.DATAMIXERoutput_step = cms.EndPath(process.DATAMIXERoutput)


# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.DATAMIXERoutput_step)

# customisation of the process


# Automatic addition of the customisation function from SimGeneral.DataMixingModule.customiseMatchRecVertex

##############################################################################

def customiseDataMixerOutput(process):
    process.DATAMIXERoutput.outputCommands.extend(
        ['keep *_hiSelectedVertex_*_*',
         'keep *_siStripDigis_VirginRaw_*'])
    return process

##############################################################################

def customiseCloneGenerator(process):
    process.generator = process.hiSignal.clone(embeddingMode=cms.bool(False))
    process.ProductionFilterSequence.replace(process.hiSignal,process.generator)
    return process

def customiseMatchRecVertex(process):
    process.load("SimGeneral.MixingModule.MatchRecVtx_cfi")
    process.matchRecVtx.heavyIonLabel = "hiSelectedVertex"
    process.pgen.replace(process.VertexSmearing,process.matchRecVtx)
    process.mixData.input.sequential = True
    process.mixData.input.fileNames = process.source.fileNames
    return process

def customiseSiStripRawDigi(process):
    process.mixData.SiStripRawDigiSource = "PILEUP" # use raw digis from "PILEUP" or "SIGNAL"
    process.SiStripDigiToRaw.FedReadoutMode = cms.string('VIRGIN_RAW') # pack rawdigis as VR

    process.stripConditions = cms.ESSource("PoolDBESSource",
        process.CondDBSetup,
        timetype = cms.untracked.string('runnumber'),
        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_STRIP'),
        toGet = cms.VPSet(cms.PSet(
          record = cms.string('SiStripFedCablingRcd'),
          tag = cms.string('SiStripFedCabling_GR10_v1_hlt')))
    )
    process.es_prefer_strips = cms.ESPrefer("PoolDBESSource","stripConditions")
    return process

def customiseCloneMatchRaw(process):
    customiseCloneGenerator(process)
    customiseMatchRecVertex(process)
    customiseSiStripRawDigi(process)
    return process

def customiseSiStripRawSignal(process):
    customiseSiStripRawDigi(process)
    process.mixData.SiStripRawDigiSource = "SIGNAL"
    process.mixData.ZDCdigiCollectionSig = ""
    process.mixData.ZDCPileInputTag = ""
    return process

##############################################################################

def customiseBeamSpot(process):
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BeamSpotObjectsRcd"),
                 tag = cms.string("Realistic2.76ATeVCollisions_STARTUP_v0_mc"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_BEAMSPOT")
                 )
        )
    return process

def customiseSiStripConditions(process):
    process.stripConditions = cms.ESSource("PoolDBESSource",
    process.CondDBSetup,
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('frontier://FrontierProd/CMS_COND_31X_STRIP'),
    toGet = cms.VPSet(
        cms.PSet(
        record = cms.string('SiStripNoisesRcd'),
        tag = cms.string('SiStripNoise_GR10_v1_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripPedestalsRcd'),
        tag = cms.string('SiStripPedestals_GR10_v1_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripFedCablingRcd'),
        tag = cms.string('SiStripFedCabling_GR10_v1_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripBadChannelRcd'),
        tag = cms.string('SiStripBadChannel_FromOnline_GR10_v1_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripLatencyRcd'),
        tag = cms.string('SiStripLatency_GR10_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripThresholdRcd'),
        tag = cms.string('SiStripThreshold_GR10_v1_hlt')
        ),
        cms.PSet(
        record = cms.string('SiStripBadFiberRcd'),
        tag = cms.string('SiStripBadChannel_FromOfflineCalibration_GR10_v2_offline')
        ),
        cms.PSet(
        record = cms.string('SiStripBadModuleRcd'),
        tag = cms.string('SiStripBadChannel_FromEfficiencyAnalysis_GR10_v1_offline')
        ),
        cms.PSet(
        record = cms.string('SiStripConfObjectRcd'),
        tag = cms.string('SiStripShiftAndCrosstalk_GR10_v1_offline')
        ),
        cms.PSet(
        record = cms.string('SiStripLorentzAngleRcd'),
        tag = cms.string('SiStripLorentzAngle_GR10_v1_offline')
        ),
        cms.PSet(
        record = cms.string('SiStripApvGain2Rcd'),
        tag = cms.string('SiStripApvGain_FromParticles_GR10_v2_offline')
        ),
        cms.PSet(
        record = cms.string('SiStripApvGainRcd'),
        tag = cms.string('SiStripApvGain_GR10_v1_hlt')
        )
        )
    )
    
    process.es_prefer_strips = cms.ESPrefer("PoolDBESSource","stripConditions")
    return process

def customiseDummyVtx(process):
    # Dummy GEN -> RECO vertex
    process.hiSelectedVertex = cms.EDProducer("GenToRecoVtxProducer",
        signalLabel=cms.InputTag("generator"))
    return process

def customiseVtxPP(process):
    # Run RECO around p+p vertex from embedded RAW-RECO event
    #process.hiCentrality.srcVertex = 'offlinePrimaryVertices'
    process.photons.primaryVertexProducer = 'offlinePrimaryVertices'
    process.hiPixel3PrimTracks.RegionFactoryPSet.RegionPSet.VertexCollection = 'offlinePrimaryVertices'
    process.hiPixel3PrimTracks.FilterPSet.VertexCollection = 'offlinePrimaryVertices'
    process.hiTracksWithLooseQuality.vertices = 'offlinePrimaryVertices'
    process.hiTracksWithTightQuality.vertices = 'offlinePrimaryVertices'
    process.hiSelectedTracks.vertices = 'offlinePrimaryVertices'
    return process

def customiseBeamRaw(process):
    customiseBeamSpot(process)
    customiseSiStripConditions(process)
    customiseDummyVtx(process)
    return process


process = customiseDataMixerOutput(process)


# End of customisation functions

