import FWCore.ParameterSet.Config as cms

process = cms.Process("ClusterAna")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.GlobalTag.globaltag = 'GR10_P_V12::All'

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFhits"),
    #nonDefaultGlauberModel = cms.string("Hydjet_Bass"),
    #nonDefaultGlauberModel = cms.string("AMPT_Organ"),
    centralitySrc = cms.InputTag("hiCentrality")
    )

from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)


import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltMinBiasHF = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltMinBiasHF.HLTPaths = ["HLT_HIMinBiasHF"]

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5CCEC6E1-31EC-DF11-8B01-0016177CA778.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5CCA2FD6-46EC-DF11-93DB-00304879BAB2.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5C72F3AA-42EC-DF11-BECB-001D09F2932B.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5A503498-2EEC-DF11-A1CC-001D09F241F0.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5A3C1D39-31EC-DF11-A102-001D09F29169.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5A30F667-35EC-DF11-90A1-001D09F25401.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/58025397-2EEC-DF11-A13B-003048F118AA.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/56EA6285-40EC-DF11-BF68-0030487CD77E.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/56BB5F98-32EC-DF11-8F99-001D09F26C5C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/566B5331-2DEC-DF11-8318-001D09F282F5.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/566B2631-2DEC-DF11-9427-001D09F2305C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/566AD6EC-3FEC-DF11-9A3B-0019DB2F3F9A.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/54F7263E-31EC-DF11-A03A-0030487CD716.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/54844A71-43EC-DF11-A7E0-003048D37580.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/545E654F-48EC-DF11-B1F4-0030487C7392.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/52FD6699-32EC-DF11-9E86-001D09F2A465.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/52D53098-32EC-DF11-A21B-001D09F25041.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/5283D39C-28EC-DF11-97D0-0030487A1990.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/52464A5D-29EC-DF11-B97A-001D09F28D54.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/52282C17-3DEC-DF11-8CA5-0030486733D8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/508BFE3B-3FEC-DF11-AF77-001D09F24DDF.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/507AA126-36EC-DF11-B909-001D09F24353.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/4ED6B1F1-41EC-DF11-8D55-0030487CD6D2.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/4EC96B79-30EC-DF11-9B65-001D09F29169.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v2/000/150/476/4EB2BF86-40EC-DF11-92AE-0030487CD76A.root'
    )
)

process.clusterAna = cms.EDAnalyzer('SiStripClusterAnalyzer',
    clusterSrc = cms.InputTag('siStripClusters'),
    etaMin = cms.double(-99.0),
    etaMax = cms.double(99.0)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('clusterWidths.root')
)


process.p = cms.Path(process.hltMinBiasHF*process.clusterAna)

