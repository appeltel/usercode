import FWCore.ParameterSet.Config as cms

process = cms.Process("PIXELFEDCUT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

process.GlobalTag.globaltag = 'START39_V7HI::All'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/6A231157-B5FA-DF11-BB60-003048322C3E.root'
    )
)

process.fiducialPixelRecHits = cms.EDProducer('PixelFiducialRemover',
    src = cms.InputTag("siPixelRecHits"),
    sides = cms.vint32( 4 ),
    disks = cms.vint32( 2 ),
    blades = cms.vint32( 4,5,6 )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.siPixelRecHits * process.fiducialPixelRecHits)

process.e = cms.EndPath(process.out)
