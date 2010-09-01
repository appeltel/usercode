#!/bin/sh

dir=RecoLocalTracker/SiStripZeroSuppression
dir2=`pwd`

cd $CMSSW_BASE/src
cvs co -r HEAD $dir
cd $dir2
cp SiStripAPVRestorer.h $CMSSW_BASE/src/$dir/interface
cp SiStripRawProcessingFactory.h $CMSSW_BASE/src/$dir/interface
cp SiStripRawProcessingAlgorithms.h $CMSSW_BASE/src/$dir/interface
cp SiStripRawProcessingFactory.cc $CMSSW_BASE/src/$dir/src
cp SiStripAPVRestorer.cc $CMSSW_BASE/src/$dir/src
cp SiStripZeroSuppression.cc $CMSSW_BASE/src/$dir/plugins
cp customiseCMN.py $CMSSW_BASE/src/$dir/python
cp SiStripZeroSuppression_cfi.py $CMSSW_BASE/src/$dir/python
