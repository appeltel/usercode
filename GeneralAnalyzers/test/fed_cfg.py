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
    fileName = cms.string('zs.root')
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
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_9_1_qgq.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_8_1_f2e.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_7_1_SKN.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_6_1_nki.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_5_1_ci8.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_4_1_1sa.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_3_1_PEm.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_2_1_VpU.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_20_1_sta.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_1_1_9pj.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_19_1_J9j.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_18_1_WsI.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_17_1_Vvy.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_16_1_Bk9.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_15_1_WXH.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_14_1_AfY.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_13_1_0xa.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_12_1_UhS.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_11_1_b25.root',
       '/store/caf/user/icali/HIExpressPhysics/ExpressLocalRecoFEDZS/9332484fceede81b8cf0f2dafa6b3bd8/RecoExpress151126-151353_HLT_HIMinBiasHfOrBSC_FEDZS_10_1_6MD.root'
  )
)


process.GlobalTag.globaltag = 'GR10_P_V12::All'

process.simpleTrack = cms.EDAnalyzer('SimpleTrackAnalyzer',
   trackSrc = cms.InputTag("hiSelectedTracks"),
   vertexSrc = cms.InputTag("hiSelectedVertex")
)

process.p = cms.Path( process.simpleTrack )
