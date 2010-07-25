#!/bin/sh

cd $CMSSW_BASE/src

cvs co UserCode/Appeltel/JEX
mv UserCode/Appeltel .
rm -r UserCode/Appeltel

scram b
