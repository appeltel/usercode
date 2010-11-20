import FWCore.ParameterSet.Config as cms

process = cms.Process('SIMPLETRACKANA')


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('t0.root')
)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    #nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( 
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_9_1_0h7.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_8_1_fwd.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_7_1_zan.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_6_1_4in.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_5_1_KKX.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_4_1_yae.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_3_1_yDz.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_2_1_NL4.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_20_1_gxS.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_1_1_kjk.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_19_1_vXo.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_18_1_MHQ.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_17_1_DAo.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_16_1_cRx.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_15_1_CX5.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_14_1_LU5.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_13_1_phV.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_12_1_gq5.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_11_1_8Cx.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalReco392patch5/2ce78842a994017e4b204f1c1bc6bb28/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_LocalReco392patch5_10_1_tVF.root' 
  )
)


process.GlobalTag.globaltag = 'GR10_P_V12::All'

process.simpleTrack = cms.EDAnalyzer('SimpleTrackAnalyzer',
   trackSrc = cms.InputTag("hiSelectedTracks"),
   vertexSrc = cms.InputTag("hiSelectedVertex")
)

process.p = cms.Path( process.simpleTrack )
