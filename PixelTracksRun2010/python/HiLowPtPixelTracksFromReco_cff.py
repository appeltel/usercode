import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HILowPtWideCutPixelTracks_cfi import *
from Appeltel.PixelTracksRun2010.HIInvertedLowPtPixelTracks_cfi import *

lowPtPixelTrackReco = cms.Sequence(
    siPixelRecHits * 
    hiLowPtPixelTracks *
    hiLowPtWideCutPixelTracks *
    hiInvertedLowPtPixelTracks 
)


