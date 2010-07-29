#!/bin/sh

cd $CMSSW_BASE/src

cvs co UserCode/Appeltel/hackedSiStripDigitizer
cvs co SimTracker/SiStripDigitizer

mv UserCode/Appeltel/hackedSiStripDigitizer/interface/SiStripDigitizer.h SimTracker/SiStripDigitizer/interface/SiStripDigitizer.h
mv UserCode/Appeltel/hackedSiStripDigitizer/plugins/SiStripDigitizer.cc SimTracker/SiStripDigitizer/plugins/SiStripDigitizer.cc
mv UserCode/Appeltel/hackedSiStripDigitizer/interface/SiStripDigitizerAlgorithm.h SimTracker/SiStripDigitizer/interface/SiStripDigitizerAlgorithm.h
mv UserCode/Appeltel/hackedSiStripDigitizer/src/SiStripDigitizerAlgorithm.h SimTracker/SiStripDigitizer/src/SiStripDigitizerAlgorithm.h

scram b

