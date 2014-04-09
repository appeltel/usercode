import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process('GEN')

# ================= var parsing ======================

options = VarParsing.VarParsing ('standard')

options.output = 'genJetSpectrum.root'
options.maxEvents = 100

options.register('processType',
                 "NSD",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Pythia process type with pT_hat range")

options.register('sqrtS',
                 5020.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Center-of-mass energy")

options.register('ptHatLow',
                 50,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Minimum pt-hat")

options.register('ptHatHigh',
                 80,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum pt-hat")

options.parseArguments()

# ================= load fragments =====================

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

# ================ options and files ====================

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.output)
)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ============= Pythia setting  ================================
#
# Uncomment one of these settings corresponding to the
# tune that you want to generate
#

from Configuration.Generator.PythiaUEZ2Settings_cfi import *
#from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
#from PythiaUEAMBT2Settings_cfi import *

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(5020.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
         pythiaUESettingsBlock,
         processParameters = cms.vstring(
                 'MSEL=1         ! High Pt QCD'
                 ),
         parameterSets = cms.vstring('pythiaUESettings',
                                     'processParameters',
                                     #'pthatLow',
                                     #'pthatHigh')
                                    )
         )
)

process.gen_step = cms.Path(process.generator
                            * process.genParticles )

# update the process parameters and c.o.m energy
from customiseGEN2_cfi import *
updatePy6ProcParameters(process.generator,options.processType,options.ptHatLow,options.ptHatHigh,options.sqrtS)

print process.generator.PythiaParameters.processParameters

# ============= Gen jet ================================
process.ak3GenJets = process.ak5GenJets.clone( rParam = 0.3 )

process.genjet_step = cms.Path(process.genJetParticles 
                               * process.ak3GenJets )

# =============== Analysis =============================
process.ak3GenJetSpectrum = cms.EDAnalyzer('GenJetCrossCheckAnalyzer',
    genJetSrc = cms.InputTag("ak3GenJets"),
    genParticleSrc = cms.InputTag("genParticles"),
    etaMin = cms.double(-1.0),
    etaMax = cms.double(1.0),
    jetRadius = cms.double(0.3),
    pthatMin = cms.double(options.ptHatLow),
    pthatMax = cms.double(options.ptHatHigh),
    ptBins = cms.vdouble( 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000 ),
    pythiaProcess = cms.string(options.processType )    
)

process.ak3GenJetSpectrum_n22_n12 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(-2.2),
    etaMax = cms.double(-1.2)
)

process.ak3GenJetSpectrum_n12_n07 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(-1.2),
    etaMax = cms.double(-0.7)
)

process.ak3GenJetSpectrum_n07_n03 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(-0.7),
    etaMax = cms.double(-0.3)
)

process.ak3GenJetSpectrum_p12_p22 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(1.2),
    etaMax = cms.double(2.2)
)


process.ak3GenJetSpectrum_p07_p12 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(0.7),
    etaMax = cms.double(1.2)
)

process.ak3GenJetSpectrum_p03_p07 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(0.3),
    etaMax = cms.double(0.7)
)


process.ak3GenJetSpectrum_n03_p03 = process.ak3GenJetSpectrum.clone(
    etaMin = cms.double(-0.3),
    etaMax = cms.double(0.3)
)

process.ana_step = cms.Path(
    process.ak3GenJetSpectrum * 
    process.ak3GenJetSpectrum_n22_n12 *
    process.ak3GenJetSpectrum_n12_n07 *
    process.ak3GenJetSpectrum_n07_n03 *
    process.ak3GenJetSpectrum_n03_p03 *
    process.ak3GenJetSpectrum_p03_p07 *
    process.ak3GenJetSpectrum_p07_p12 *
    process.ak3GenJetSpectrum_p12_p22 
)

process.schedule = cms.Schedule(process.gen_step,
                                process.genjet_step,
                                process.ana_step)
