#ifndef MergedTrackCorrections_MergedTrackCorrector
#define MergedTrackCorrections_MergedTrackCorrector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
class MergedTrackCorrector
{
  public:

   MergedTrackCorrector();
   ~MergedTrackCorrector();

   void loadEffic();
   void loadEffic(const char *);
   void loadFakes();
   void loadFakes(const char *);

//   float weight( reco::Track &, int );
//   float fakeRate( reco::Track &, int );
//   float effic( reco::Track &, int );

   float weight( float , float , int );
   float fakeRate( float , float , int );
   float effic( float , float , int );

  private:

   float effic_[3][19][3][4];
   float fake_[3][19][3][4];

};


#endif

