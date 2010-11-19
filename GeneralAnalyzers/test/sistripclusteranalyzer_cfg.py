import FWCore.ParameterSet.Config as cms

process = cms.Process("ClusterAna")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

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
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasHFOrBSC = process.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_Core"]

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/CC18ADCB-36EE-DF11-85C7-001D09F29849.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/CC17F10E-3DEE-DF11-9E28-000423D94E70.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/CA02DEB1-28EE-DF11-BD3F-0030487CD840.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/C610E6FF-27EE-DF11-9FD7-003048F118E0.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/C22D7A87-39EE-DF11-8680-0019DB2F3F9A.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/C06D8CE0-3FEE-DF11-A3E8-0030487C7828.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/B87D5502-2FEE-DF11-911F-001617E30CC8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/B4062DE0-3FEE-DF11-A0ED-0030487CD812.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/A8AE2A43-20EE-DF11-BC1F-001D09F24259.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/A641636E-30EE-DF11-90F3-001D09F2B2CF.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/A2E5D188-39EE-DF11-91D2-001617C3B79A.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/A0B70B6C-29EE-DF11-B1EB-0030487A18D8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/9E8E46B5-2FEE-DF11-A76B-0030487C8CB6.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/984A589A-2DEE-DF11-91D1-003048F118C2.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/966E805D-43EE-DF11-B956-000423D33970.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/96163800-28EE-DF11-84E2-001617C3B5D8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/94151CB0-28EE-DF11-BE74-001617DBD316.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/8C981B6C-29EE-DF11-BB68-0030487CD6F2.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/8C618DB3-42EE-DF11-B1C4-003048F11114.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/8C007E27-3FEE-DF11-8FB0-000423D996C8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/8A657101-28EE-DF11-A079-001D09F23C73.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/88F9007C-45EE-DF11-95EF-000423D98B6C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/881A6C13-36EE-DF11-9E6B-001D09F24DDA.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/7C68EBA1-3BEE-DF11-A5A0-001617E30F48.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/76CF8D1E-2AEE-DF11-A2AE-001D09F252F3.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/7625E32F-2CEE-DF11-8078-001617E30D12.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/727DD35B-43EE-DF11-8CCC-001617E30D40.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/70E678CB-44EE-DF11-97DC-003048F1110E.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/701E7FB2-42EE-DF11-8303-003048F01E88.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/6CB9C177-35EE-DF11-A3C0-000423D98B6C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/686C0E43-41EE-DF11-8E9A-001617E30D00.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/64615C55-27EE-DF11-97F2-001D09F2924F.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/5AD31018-34EE-DF11-88B8-001617DBD472.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/5263771A-2AEE-DF11-B96C-001D09F24448.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/4E3939B4-42EE-DF11-ACCA-000423D9890C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/4C22B146-33EE-DF11-85EE-0030487C8CB6.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/4898CB57-22EE-DF11-B1E6-001D09F24DDA.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/483D008D-40EE-DF11-8D80-000423D996C8.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/463F8F73-3EEE-DF11-9AB4-001617E30F48.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/42E514B2-28EE-DF11-AC70-0030487CD6DA.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/42501E43-41EE-DF11-8778-001617E30F4C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/3ABBD7E5-2CEE-DF11-9FDB-001617E30D4A.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/384C3D7C-24EE-DF11-A76B-001D09F24303.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/34B534D3-2AEE-DF11-A03A-0016177CA778.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/3496E012-36EE-DF11-A13B-001D09F231C9.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/32D78347-33EE-DF11-A160-0030487A3DE0.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/32BE7B00-28EE-DF11-80E6-001617DBD224.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/32384D00-28EE-DF11-B359-001617C3B6E2.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2AD13FDA-38EE-DF11-9846-003048D2BB90.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2ABEF8A3-32EE-DF11-A198-001D09F292D1.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2A5E8842-41EE-DF11-A262-000423D98B6C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2A14D410-44EE-DF11-949B-003048F118DE.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/2660FB66-29EE-DF11-BFB7-0030487CAEAC.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/229D3579-35EE-DF11-8736-003048D37538.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/221699A6-26EE-DF11-BA47-001D09F252F3.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/1EBC9184-53EE-DF11-975B-0030487C7E18.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/1C673F67-29EE-DF11-9BB9-0030487CD7EE.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/12DCB7CB-36EE-DF11-8916-001D09F24493.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/0EAB1AEF-3AEE-DF11-820F-001617E30F48.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/0A0F139A-2DEE-DF11-9F82-003048F118DE.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/0818CFB2-42EE-DF11-9864-003048F1C58C.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/065F7412-44EE-DF11-AD82-003048F1183E.root',
        '/store/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/150/887/049B1FD3-2AEE-DF11-8EC2-001617E30CC8.root'

    )
)

process.clusterAna = cms.EDAnalyzer('SiStripClusterAnalyzer',
    stripClusterSrc = cms.InputTag('siStripClusters'),
    pixelClusterSrc = cms.InputTag('siPixelClusters'),
    trackSrc = cms.InputTag('hiSelectedTracks'),
    etaMin = cms.double(-99.0),
    etaMax = cms.double(99.0),
    onlyCount = cms.bool(False)

)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('clusterWidths.root')
)

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")


process.path = cms.Path(process.hltMinBiasHFOrBSC*process.collisionEventSelection*process.clusterAna)
