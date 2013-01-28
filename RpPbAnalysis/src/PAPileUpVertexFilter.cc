/** \class PAPileUpVertexFilter
 *
 *  
 *  This class is an EDFilter for rejecting  pPb events with multiple collisions
 *  (pileup). This is performed by looking at the characteristics of the 
 *  reconstructed vertices.
 *
 *  $Date: 2013/01/27 16:44:31 $
 *  $Revision: 1.2 $
 *
 *  \author E. Appelt - Vanderbilt University
 *
 */


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <vector>
#include <algorithm>

#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class PAPileUpVertexFilter : public edm::EDFilter {
    public:
       explicit PAPileUpVertexFilter(const edm::ParameterSet&);
       ~PAPileUpVertexFilter();
       virtual void endJob() ;

       virtual bool filter(edm::Event&, const edm::EventSetup&);

    private:

       edm::InputTag vtxSrc_;
       double dxyCut_;
       int trkCut_;
       std::vector<double> dzCutByNtrk_;
       
       

};


PAPileUpVertexFilter::PAPileUpVertexFilter(const edm::ParameterSet& iConfig) :
vtxSrc_(iConfig.getParameter<edm::InputTag>("vtxSrc")),
dxyCut_(iConfig.getParameter<double>("dxyCut")),
trkCut_(iConfig.getParameter<int>("trkCut")),
dzCutByNtrk_(iConfig.getParameter<std::vector<double> >("dzCutByNtrk"))
{
}

PAPileUpVertexFilter::~PAPileUpVertexFilter()
{
}

bool
PAPileUpVertexFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   bool accepted = true;

   Handle<std::vector<reco::Vertex> > vcol;
   iEvent.getByLabel(vtxSrc_, vcol);

   std::vector<reco::Vertex> vsorted = *vcol;
   // sort the vertcies by number of tracks in descending order
   std::sort( vsorted.begin(), vsorted.end(), 
              []( reco::Vertex a, reco::Vertex b) 
   {
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   });

   // check additional PVs
   for (unsigned int i =1; i<vsorted.size(); i++)
   {
     double dz = fabs( vsorted[i].z() - vsorted[0].z() );
     double dx = fabs( vsorted[i].x() - vsorted[0].x() );
     double dy = fabs( vsorted[i].y() - vsorted[0].y() );
     double dxy  = sqrt ( dx*dx + dy*dy );
     double nTrk = vsorted[i].tracksSize();

     // only filter for small dxy
     if ( dxy > dxyCut_ ) continue;

     // filter if above absolute max number of tracks for 2nd PV
     if ( nTrk > trkCut_ )
       accepted = false;

     // for smaller nTrk, filter on dz by number of tracks
     for( unsigned int j=0; j<dzCutByNtrk_.size() ; j++)
     {
       if ( nTrk == (int)j+1 && dz > dzCutByNtrk_[j] ) 
         accepted = false;
     }   

   }

   return accepted;

}

void
PAPileUpVertexFilter::endJob()
{
}

DEFINE_FWK_MODULE(PAPileUpVertexFilter);
