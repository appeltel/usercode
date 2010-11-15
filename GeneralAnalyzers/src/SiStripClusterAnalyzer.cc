// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

//
// class declaration
//

class SiStripClusterAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripClusterAnalyzer(const edm::ParameterSet&);
      ~SiStripClusterAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag stripClusterSrc_;
      edm::InputTag pixelClusterSrc_;

      double etaMin_;
      double etaMax_;

      bool onlyCount_;

      int eventsProcessed_;

      CentralityProvider * centrality_;
 
      TH1I* clusterwidth_[10];

      TH1I* centbins_;

      TH2F* widthDetId_;

      TH2F* clusterCountPixelStrip_;

};

SiStripClusterAnalyzer::SiStripClusterAnalyzer(const edm::ParameterSet& iConfig):
stripClusterSrc_(iConfig.getParameter<edm::InputTag>("stripClusterSrc")),
pixelClusterSrc_(iConfig.getParameter<edm::InputTag>("pixelClusterSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
onlyCount_(iConfig.getParameter<bool>("onlyCount"))
{

  edm::Service<TFileService> fs;

  for( int i = 0; i<10; i++)
  {
    clusterwidth_[i] = fs->make<TH1I>(Form("clusterwidth%d",i),Form("clusterwidth%d",i), 250, 1, 250 );
  }
  centbins_ = fs->make<TH1I>("centbins","centbins", 40, 0, 39);

  widthDetId_ = fs->make<TH2F>("widthDetId","Cluster Widths vs Detector ID", 250, 1., 250., 500, 0., 66000.0); 

  clusterCountPixelStrip_ = fs->make<TH2F>("clusterCountPixelStrip","Pixel Clusters versus Strip Clusters",
                                          200, 0., 100000., 200, 0., 50000.);

  // safety
  centrality_ = 0;
}


SiStripClusterAnalyzer::~SiStripClusterAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiStripClusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   
   // Get Centrality bin for the event

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);

   int pixelCount = 0;
   int stripCount = 0;

   edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
   iEvent.getByLabel(stripClusterSrc_, clusters);

   stripCount += clusters->dataSize();

   if( ! onlyCount_ )
   {
     edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters->begin();
     for ( ; it != clusters->end(); ++it )
     {
       for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
       {
          int firststrip = clus->firstStrip();
          int laststrip = firststrip + clus->amplitudes().size();
          clusterwidth_[bin/4]->Fill( laststrip-firststrip);     

          widthDetId_->Fill( laststrip-firststrip, it->detId() % 65536 );


       }
     }
   }

   edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pclusters;
   iEvent.getByLabel(pixelClusterSrc_, pclusters);

   pixelCount += pclusters->dataSize();


   clusterCountPixelStrip_->Fill( pixelCount, stripCount );   

   eventsProcessed_++;

}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripClusterAnalyzer::beginJob()
{
   eventsProcessed_ = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripClusterAnalyzer::endJob() {
 
  std::cout << eventsProcessed_ << " events were processed.\n";

}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripClusterAnalyzer);

