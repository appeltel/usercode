import FWCore.ParameterSet.Config as cms

process = cms.Process('GENPARTANA')


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('particleSpecies.root')
)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/46B060AD-C9FA-DF11-B18A-485B3919F121.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/1A80A9F9-C7FA-DF11-9587-0023541E7797.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/1A3F7B7E-D5FA-DF11-BA8A-001A92971C6B.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/120C00AB-C8FA-DF11-85CA-0023546BA20F.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/0A456C69-9FFA-DF11-B353-003048322C3A.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0002/00D94EBD-9DFA-DF11-9C83-001E68864975.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FED91C44-9AFA-DF11-B174-0002C90B3976.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FED5DA39-8BFA-DF11-9345-001A9227D3C1.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FEB4E23E-90FA-DF11-B2D1-00304865C254.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FEA0CDF0-98FA-DF11-8FCC-00304865C456.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FE8694B5-9FFA-DF11-98AD-00238BCE45FA.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FE75EC30-9CFA-DF11-8F0F-00238BCE44F0.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FE40A0D0-91FA-DF11-816C-001E6849D21C.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FCE32635-8DFA-DF11-AAA8-00304865C45A.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FCB8F72D-8FFA-DF11-A57E-00304867FEAB.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FC954A43-95FA-DF11-A896-0002C90B7F3C.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FC6567FD-93FA-DF11-9CF7-00304865C29A.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FC0A55B1-8EFA-DF11-9811-842B2B185476.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FC05F6AA-94FA-DF11-9765-00304865C49C.root',
        '/store/himc/Fall10/Hydjet_Bass_MinBias_2760GeV/GEN-SIM-RECODEBUG/START39_V7HI-v1/0000/FC00D543-95FA-DF11-9FF6-0002C90B7428.root'
  )
)


process.GlobalTag.globaltag = 'START39_V7HI::All'

process.genPart = cms.EDAnalyzer('GenParticleAnalyzer',
   genSrc = cms.InputTag("hiGenParticles"),
   etaMax = cms.double(0.8)
)

process.p = cms.Path( process.genPart )
