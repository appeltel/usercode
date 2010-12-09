#ifndef MergedTrackCorrections_UnmergedTrackCorrector
#define MergedTrackCorrections_UnmergedTrackCorrector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
class UnmergedTrackCorrector
{
  public:

   UnmergedTrackCorrector();
   ~UnmergedTrackCorrector();

   void loadMap();
   void loadMap(const char *);

//   float weight( reco::Track &, int );
//   float fakeRate( reco::Track &, int );
//   float effic( reco::Track &, int );

   float weight( float , float , int );
   float fakeRate( float , float , int );
   float effic( float , float , int );

  private:

   double  m_[2][19][3][6];

};


#endif

