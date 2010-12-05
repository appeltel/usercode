import FWCore.ParameterSet.Config as cms

process = cms.Process('RAWMIX')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedRealistic7TeVCollision_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.DigiDM_cff')
process.load('Configuration.StandardSequences.DataMixerDataOnSim_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
process.load('Configuration.StandardSequences.DigiToRawDM_cff')
process.load('HLTrigger.Configuration.HLT_HIon_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

################################################################################
# Input Variable Parsing
#

import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

#ivars.register ('randomNumber',
#                 1,
#                ivars.multiplicity.singleton,
#                ivars.varType.int,
#                "Random Seed")



#ivars.randomNumber = 1

ivars.inputFiles = 'file:/scratch/appelte1/HIRAWSelected/selectedRAW_14.root'

ivars.outputFile = 'pgunMix_test.root'

ivars.parseArguments()


#################################################################################
# Random seed
#

process.load('Configuration/StandardSequences/SimulationRandomNumberGeneratorSeeds_cff')
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

########################################################################################
# Timing and memory services
#

process.Timing = cms.Service("Timing")


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(ivars.inputFiles)
)


process.options = cms.untracked.PSet(

)

# Output definition

process.RAWDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RAWDEBUGEventContent.outputCommands,
    fileName = cms.untracked.string(ivars.outputFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAWDEBUG')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START39_V7HI::All'

process.hiSignal = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(5.),
        MinPt = cms.double(0.2),
        PartID = cms.vint32(
          211, 211, 211, 211, 211,
          211, 211, 211, 211, 211,
          -211, -211, -211, -211, -211,
          -211, -211, -211, -211, -211
        ),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359) ## in radians

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('20 pions pt 0.2 to 5'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)

process.ProductionFilterSequence = cms.Sequence(process.hiSignal)



# Path and EndPath definitions

process.generation_step = cms.Path(process.pgen)

process.simulation_step = cms.Path(process.psim)

process.digitisation_step = cms.Path(process.pdigi)

process.datamixing_step = cms.Path(process.pdatamix)

process.L1simulation_step = cms.Path(process.SimL1Emulator)

process.digi2raw_step = cms.Path(process.DigiToRaw)

process.endjob_step = cms.Path(process.endOfProcess)

process.RAWDEBUGoutput_step = cms.EndPath(process.RAWDEBUGoutput)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.datamixing_step,process.L1simulation_step,process.digi2raw_step)

process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWDEBUGoutput_step])

# special treatment in case of production filter sequence
for path in process.paths: 
  getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq


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

def customiseBeamRaw(process):
    customiseBeamSpot(process)
    customiseSiStripConditions(process)
    customiseDummyVtx(process)
    return process


process = customiseCloneMatchRaw(process)


# Automatic addition of the customisation function from SimGeneral.DataMixingModule.DataMixer_DataConditions_3_8_X_data2010
def customise(process):

#
# IOV set based on GlobalTag GR_R_35X_V8B
#
# placeholder !!!!!! replace with the actual run number of
# the real run to be overlaid

    process.source.firstRun = cms.untracked.uint32(132599)

    process.ecalConditions1 = cms.ESSource("PoolDBESSource",                                          
         process.CondDBSetup,                                                                         
         timetype = cms.string('runnumber'),                                                          
         toGet = cms.VPSet(                                                                           
             cms.PSet(                                                                                
        record = cms.string('EcalADCToGeVConstantRcd'),                                               
        tag = cms.string('EcalADCToGeVConstant_v6_offline')
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalChannelStatusRcd'),                                                  
        tag = cms.string('EcalChannelStatus_v04_offline')                                   
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalGainRatiosRcd'),                                                     
        tag = cms.string('EcalGainRatio_TestPulse2009_offline')                                      
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalIntercalibConstantsRcd'),                                            
        tag = cms.string('EcalIntercalibConstants_v6_offline')                                 
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalIntercalibErrorsRcd'),                                               
        tag = cms.string('EcalIntercalibErrors_mc')                                                   
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalMappingElectronicsRcd'),                                             
        tag = cms.string('EcalMappingElectronics_EEMap')                                              
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalPedestalsRcd'),                                                      
        tag = cms.string('EcalPedestals_2009runs_hlt')                                                
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalTBWeightsRcd'),                                                      
        tag = cms.string('EcalTBWeights_EBEE_v01_offline')                                     
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalTimeCalibConstantsRcd'),                                             
        tag = cms.string('EcalTimeCalibConstants_v02_offline')
        ),                                                                                            
             cms.PSet(                                                                                
        record = cms.string('EcalWeightXtalGroupsRcd'),                                               
        tag = cms.string('EcalWeightXtalGroups_EBEE_offline')                                   
        ),                                                                                            
             cms.PSet(                                                                   
        record = cms.string('EcalLaserAPDPNRatiosRcd'),                                               
        tag = cms.string('EcalLaserAPDPNRatios_p1p2p3_v2_mc')                                        
        ),                                                                                            
             ),                                                                                       
        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_ECAL'),                            
              authenticationMethod = cms.untracked.uint32(0)                                          
    )                                                                                                 
                                                                                                      

    process.ecalConditions2 = cms.ESSource("PoolDBESSource",
                                           process.CondDBSetup,
                                           timetype = cms.string('runnumber'),
                                           toGet = cms.VPSet(
        cms.PSet(
        record = cms.string('EcalTPGCrystalStatusRcd'),
        tag = cms.string('EcalTPGCrystalStatus_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGFineGrainEBGroupRcd'),
        tag = cms.string('EcalTPGFineGrainEBGroup_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGFineGrainEBIdMapRcd'),
        tag = cms.string('EcalTPGFineGrainEBIdMap_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGFineGrainStripEERcd'),
        tag = cms.string('EcalTPGFineGrainStripEE_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGFineGrainTowerEERcd'),
        tag = cms.string('EcalTPGFineGrainTowerEE_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGLinearizationConstRcd'),
        tag = cms.string('EcalTPGLinearizationConst_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGLutGroupRcd'),
        tag = cms.string('EcalTPGLutGroup_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGLutIdMapRcd'),
        tag = cms.string('EcalTPGLutIdMap_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGPedestalsRcd'),
        tag = cms.string('EcalTPGPedestals_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGPhysicsConstRcd'),
        tag = cms.string('EcalTPGPhysicsConst_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGSlidingWindowRcd'),
        tag = cms.string('EcalTPGSlidingWindow_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGTowerStatusRcd'),
        tag = cms.string('EcalTPGTowerStatus_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGWeightGroupRcd'),
        tag = cms.string('EcalTPGWeightGroup_v2_hlt')
        ),
        cms.PSet(
        record = cms.string('EcalTPGWeightIdMapRcd'),
        tag = cms.string('EcalTPGWeightIdMap_v2_hlt')
        ),
        ),
        connect = cms.string('frontier://FrontierProd/CMS_COND_34X_ECAL'),
               authenticationMethod = cms.untracked.uint32(0)
    )

    process.es_prefer_ecal1 = cms.ESPrefer("PoolDBESSource","ecalConditions1")                        
    process.es_prefer_ecal2 = cms.ESPrefer("PoolDBESSource","ecalConditions2")                        

                                                                                                      
    process.hcalConditions = cms.ESSource("PoolDBESSource",                                           
                                          process.CondDBSetup,                          
                                          timetype = cms.string('runnumber'),                         
                                          toGet = cms.VPSet(                                          
        cms.PSet(                                                                                     
        record = cms.string('HcalChannelQualityRcd'),                                                 
        tag = cms.string('HcalChannelQuality_v2.10_offline')                                          
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalElectronicsMapRcd'),                                                 
        tag = cms.string('HcalElectronicsMap_v7.03_hlt')                                              
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalGainsRcd'),                                                          
        tag = cms.string('HcalGains_v2.32_offline')                                                   
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalL1TriggerObjectsRcd'),                                               
        tag = cms.string('HcalL1TriggerObjects_v1.00_hlt')                                            
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalLUTCorrsRcd'),                                                       
        tag = cms.string('HcalLUTCorrs_v1.01_hlt')                                                    
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalPedestalsRcd'),                                                      
        tag = cms.string('HcalPedestals_ADC_v9.12_offline')                                        
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalPedestalWidthsRcd'),                                                 
        tag = cms.string('HcalPedestalWidths_ADC_v7.01_hlt')                                          
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalPFCorrsRcd'),                                                        
        tag = cms.string('HcalPFCorrs_v2.00_express')                                                 
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalQIEDataRcd'),                                                        
        tag = cms.string('HcalQIEData_NormalMode_v7.00_hlt')                                          
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalRespCorrsRcd'),                                                      
        tag = cms.string('HcalRespCorrs_v1.02_express')                                               
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalTimeCorrsRcd'),                                                      
        tag = cms.string('HcalTimeCorrs_v1.00_express')                                               
        ),                                                                                            
        cms.PSet(                                                                                     
        record = cms.string('HcalZSThresholdsRcd'),                                                   
        tag = cms.string('HcalZSThresholds_v1.01_hlt')                                                
        ),                                                                                            
        ),                                                                                            
             connect = cms.string('frontier://FrontierProd/CMS_COND_31X_HCAL'),                       
                      authenticationMethod = cms.untracked.uint32(0)                                  
    )                                                                                                 
                                                                                                      
    process.es_prefer_hcal = cms.ESPrefer("PoolDBESSource","hcalConditions")                          
                                                                                                      
    try: 
        process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]  
    except:
        return(process)
 
    return(process)


process = customise(process)


# End of customisation functions

process.mixData.input.fileNames = process.source.fileNames



