#!/bin/sh

cd $CMSSW_BASE/src

# event selection and centrality 
cvs co -r pPbProd_v1 DataFormats/HeavyIonEvent
cvs co -r pPbProd_v1 RecoHI/HiCentralityAlgos
cvs co -r pPbProd_v1 HeavyIonsAnalysis/Configuration

# forest
cvs co -d HiForest -r HiForest_V2_03_01 UserCode/CmsHi/HiForest/V2

scram b
