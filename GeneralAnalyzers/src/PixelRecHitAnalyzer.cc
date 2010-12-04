// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include "TMath.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include <DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h>
#include <DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h>
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"



//
// class declaration
//

class PixelRecHitAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelRecHitAnalyzer(const edm::ParameterSet&);
      ~PixelRecHitAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag pixelRecHitSrc_;
      edm::InputTag vertexSrc_;

      CentralityProvider * centrality_;

      int eventsProcessed_;
 
      TH1I* centbins_;
      TH2F* clusterWidthEta_;
      TH2F* centralityOccupancy_;

};

PixelRecHitAnalyzer::PixelRecHitAnalyzer(const edm::ParameterSet& iConfig):
pixelRecHitSrc_(iConfig.getParameter<edm::InputTag>("pixelRecHitSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{
  edm::Service<TFileService> fs;

  centralityOccupancy_ = fs->make<TH2F>("centOccupancy","Good Pixel Clusters by Centrality Bin",40,0,40,500,0,70000);
  clusterWidthEta_ = fs->make<TH2F>("clusterWidthEta","Cluster Size Along Beam Line vs Eta",100,-3,3, 100, 0, 20);
  centbins_ = fs->make<TH1I>("centbins","centbins", 40, 0, 39);
  // safety
  centrality_ = 0;
}


PixelRecHitAnalyzer::~PixelRecHitAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelRecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Get Centrality bin for the event

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);

   edm::ESHandle<TrackerGeometry> geom;
   iSetup.get<TrackerDigiGeometryRecord>().get(geom);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   double rver_x =0.0, rver_y =0.0, rver_z =0.0;
  
   if(vertex->size() > 0) 
   {
     rver_x = vertex->begin()->x();
     rver_y = vertex->begin()->y();
     rver_z = vertex->begin()->z();
   }

   int pixelCount = 0;

   edm::Handle<SiPixelRecHitCollection> pixelRecHits;
   iEvent.getByLabel(pixelRecHitSrc_, pixelRecHits);

   for(SiPixelRecHitCollection::DataContainer::const_iterator iPRH = pixelRecHits->data().begin(); iPRH != pixelRecHits->data().end(); ++iPRH ) {
     const PixelGeomDetUnit* pgdu = dynamic_cast<const PixelGeomDetUnit*>(geom->idToDetUnit((*iPRH).geographicalId()));
     PXBDetId pdetId = PXBDetId((*iPRH).geographicalId());
     //int layer=pdetId.layer();
     GlobalPoint globalPosition = pgdu->toGlobal((*iPRH).localPosition());
     float x = fabs(globalPosition.x() - rver_x);
     float y = fabs(globalPosition.y() - rver_y);
     double R = TMath::Sqrt(x*x + y*y);
     double eta_corr = -log(tan(atan2(R,(globalPosition.z()-rver_z))/2.));
     int clustSizeY = (iPRH)->cluster()->sizeY();
     if(pgdu->subDetector() == GeomDetEnumerators::PixelBarrel) 
     {
       clusterWidthEta_->Fill( eta_corr, clustSizeY );
       if(fabs(eta_corr)<3.0 && fabs(eta_corr)>=2.5 && clustSizeY>8.5)   
         pixelCount++;
       if(fabs(eta_corr)<2.5 && fabs(eta_corr)>=2.0 && clustSizeY>5.5)   
         pixelCount++;
       if(fabs(eta_corr)<2.0 && fabs(eta_corr)>=1.5 && clustSizeY>3.5)   
         pixelCount++;
       if(fabs(eta_corr)<1.5 && fabs(eta_corr)>=1.0 && clustSizeY>2.5)   
         pixelCount++;
       if(fabs(eta_corr)<1.0 && fabs(eta_corr)>=0.5 && clustSizeY>1.5)   
         pixelCount++;
       if(fabs(eta_corr)<0.5)   
         pixelCount++;
     }
   }

   centralityOccupancy_->Fill(bin, pixelCount);

   eventsProcessed_++;

}


// ------------ method called once each job just before starting event loop  ------------
void 
PixelRecHitAnalyzer::beginJob()
{
   eventsProcessed_ = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelRecHitAnalyzer::endJob() {
 
  std::cout << eventsProcessed_ << " events were processed.\n";

}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelRecHitAnalyzer);

