#!/bin/sh

dir=RecoLocalTracker/SiStripZeroSuppression
dir2=`pwd`

cd $CMSSW_BASE/src
cvs co -r HEAD $dir
cd $dir2
cp *CMNSubtractor.h $CMSSW_BASE/src/$dir/interface
cp *CMNSubtractor.cc SiStripRawProcessingFactory.cc $CMSSW_BASE/src/$dir/src
cp customiseCMN.py $CMSSW_BASE/src/$dir/python

