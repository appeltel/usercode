#!/bin/sh

cd $CMSSW_BASE/src

# event selection (replace with tag once available)
cvs co HeavyIonsAnalysis/Configuration


scram b
