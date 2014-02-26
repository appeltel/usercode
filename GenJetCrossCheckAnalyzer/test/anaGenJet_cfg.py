import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process('GEN')

# ================= var parsing ======================

options = VarParsing.VarParsing ('standard')

options.output = 'genJetSpectrum.root'
options.maxEvents = 100

options.register('processType',
                 "NSD_50_to_80",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Pythia process type with pT_hat range")

options.register('sqrtS',
                 5020.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Center-of-mass energy")

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
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.output)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ============= Pythia setting  ================================
from Configuration.Generator.PythiaUEZ2Settings_cfi import *
#from Configuration.Generator.PythiaUEZ2starSettings_cfi import *

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(5020.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
         pythiaUESettingsBlock,
         processParameters = cms.vstring(
                 'MSEL=1         ! High Pt QCD',
                 'CKIN(3)=50     ! Pt hat lower cut',
                 'CKIN(4)=80     ! Pt hat upper cut'
                 ),
         parameterSets = cms.vstring('pythiaUESettings',
                                     'processParameters')
         )
)

process.gen_step = cms.Path(process.generator
                            * process.genParticles )

# update the process parameters and c.o.m energy
from customiseGEN_cfi import *
updatePy6ProcParameters(process.generator,options.processType,options.sqrtS)

# ============= Gen jet ================================
process.ak3GenJets = process.ak5GenJets.clone( rParam = 0.3 )

process.genjet_step = cms.Path(process.genJetParticles 
                               * process.ak3GenJets )

# =============== Analysis =============================
process.ak3GenJetSpectrum = cms.EDAnalyzer('GenJetCrossCheckAnalyzer',
    genJetSrc = cms.InputTag("ak3GenJets"),
    etaMin = cms.double(-1.0),
    etaMax = cms.double(1.0),
    ptBins = cms.vdouble( 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000 ),
    pythiaProcess = cms.string(options.processType )    
)

process.ana_step = cms.Path(process.ak3GenJetSpectrum)

process.schedule = cms.Schedule(process.gen_step,
                                process.genjet_step,
                                process.ana_step)
