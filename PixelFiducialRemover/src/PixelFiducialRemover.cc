// -*- C++ -*-
//
// Package:    PixelFiducialRemover
// Class:      PixelFiducialRemover
// 
/**\class PixelFiducialRemover PixelFiducialRemover.cc Appeltel/PixelFiducialRemover/src/PixelFiducialRemover.cc

 Description: Performs fiducial cuts on a collection of pixel rechits

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric A. Appelt
//         Created:  Fri Feb  4 14:41:27 EST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include <DataFormats/SiPixelDetId/interface/PXFDetId.h>
#include <DataFormats/SiPixelDetId/interface/PixelEndcapName.h>
#include <DataFormats/DetId/interface/DetId.h>

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

//
// class declaration
//

class PixelFiducialRemover : public edm::EDProducer {
   public:
      explicit PixelFiducialRemover(const edm::ParameterSet&);
      ~PixelFiducialRemover();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      edm::InputTag src_;
      std::vector<int> sides_;
      std::vector<int> disks_;
      std::vector<int> blades_;
      
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PixelFiducialRemover::PixelFiducialRemover(const edm::ParameterSet& iConfig)
  :
  src_( iConfig.getParameter<edm::InputTag>( "src" ) ),
  sides_( iConfig.getParameter<std::vector<int> >( "sides") ),
  disks_( iConfig.getParameter<std::vector<int> >( "disks") ),
  blades_( iConfig.getParameter<std::vector<int> >( "blades" ) ) 
{
   produces<SiPixelRecHitCollection>();
}


PixelFiducialRemover::~PixelFiducialRemover()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PixelFiducialRemover::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<SiPixelRecHitCollectionNew> input;
   iEvent.getByLabel( src_, input );

   std::auto_ptr<SiPixelRecHitCollectionNew> output(new SiPixelRecHitCollectionNew);

   SiPixelRecHitCollectionNew::const_iterator module = input->begin();
   for( ; module != input->end(); ++module )
   {

     DetId id = module->detId();
     bool keepModule = true;
     if ( (unsigned int)(id.subdetId()) == PixelSubdetector::PixelEndcap )
     {
//       PXFDetId pixel_mod( module->rawId() );

       std::vector<int>::const_iterator side;
       for( side = sides_.begin(); side != sides_.end(); ++side )
       {
         std::vector<int>::const_iterator disk;
         for( disk = disks_.begin(); disk != disks_.end(); ++disk )
         {
           std::vector<int>::const_iterator blade;
           for( blade = blades_.begin(); blade != blades_.end(); ++blade )
           {
             if( ( *blade == PixelEndcapName(id.rawId()).bladeName() ) &&
                 ( *disk == PixelEndcapName(id.rawId()).diskName() ) &&
                 ( *side == PixelEndcapName(id.rawId()).halfCylinder() ) )
             {
               keepModule = false;
             }
           }
         }
       }
     }

     if( keepModule ) 
     {
       SiPixelRecHitCollectionNew::FastFiller recHitsOnDetUnit(*output,id.rawId());
       edmNew::DetSet<SiPixelRecHit>::const_iterator hitIt = module->begin(), hitEnd = module->end();	
       for ( ; hitIt != hitEnd; hitIt++) 
       {
         recHitsOnDetUnit.push_back(*hitIt);
       }
     }

   }

   iEvent.put(output);
}

// ------------ method called once each job just before starting event loop  ------------
void 
PixelFiducialRemover::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelFiducialRemover::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelFiducialRemover);
