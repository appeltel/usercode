#!/bin/sh

cd $CMSSW_BASE/src

# event selection and centrality 
cvs co -r pPbProd_v1 DataFormats/HeavyIonEvent
cvs co -r pPbProd_v1 RecoHI/HiCentralityAlgos
cvs co -r pPbProd_v1 HeavyIonsAnalysis/Configuration

scram b
