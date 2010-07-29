#!/bin/sh

cd $CMSSW_BASE/src

cvs co UserCode/Appeltel/hackedSiStripDigitizer
cvs co SimTracker/SiStripDigitizer

mv UserCode/Appeltel/hackedSiStripDigitizer/SiStripDigitizer.h SimTracker/SiStripDigitizer/interface/SiStripDigitizer.h
mv UserCode/Appeltel/hackedSiStripDigitizer/SiStripDigitizer.cc SimTracker/SiStripDigitizer/plugins/SiStripDigitizer.cc
mv UserCode/Appeltel/hackedSiStripDigitizer/SiStripDigitizerAlgorithm.h SimTracker/SiStripDigitizer/interface/SiStripDigitizerAlgorithm.h
mv UserCode/Appeltel/hackedSiStripDigitizer/SiStripDigitizerAlgorithm.cc SimTracker/SiStripDigitizer/src/SiStripDigitizerAlgorithm.cc

scram b

