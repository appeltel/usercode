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

#include <DataFormats/TrackReco/interface/Track.h> 
#include <DataFormats/TrackReco/interface/TrackFwd.h> 

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
      edm::InputTag trackSrc_;

      double etaMin_;
      double etaMax_;

      bool onlyCount_;

      int eventsProcessed_;

      CentralityProvider * centrality_;
 
      TH1I* clusterwidth_[10];
      TH1I* clusterADCtot_[10];

      TH1I* clusterADCtotTIB_;
      TH1I* clusterADCtotTOB_;
      TH1I* clusterADCtotTID_;
      TH1I* clusterADCtotTEC_;

      TH1I* centbins_;

      TH2F* widthDetId_;

      TH2I* trackClusterModules_;

      TH2F* clusterCountPixelStrip_;

};

SiStripClusterAnalyzer::SiStripClusterAnalyzer(const edm::ParameterSet& iConfig):
stripClusterSrc_(iConfig.getParameter<edm::InputTag>("stripClusterSrc")),
pixelClusterSrc_(iConfig.getParameter<edm::InputTag>("pixelClusterSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
onlyCount_(iConfig.getParameter<bool>("onlyCount"))
{

  char centStrings[10][256] = { "0-10%", "10-20%", "20-30%", "30-40%",
       "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };


  edm::Service<TFileService> fs;

  for( int i = 0; i<10; i++)
  {
    clusterwidth_[i] = fs->make<TH1I>(Form("clusterwidth%d",i),
        Form("Cluster Widths for Centrality Range %s",centStrings[i]), 250, 1, 250 );
    clusterADCtot_[i] = fs->make<TH1I>(Form("clusterADCtot%d",i),
        Form("Uncorrected Cluster ADC Count for Centrality Range %s",centStrings[i]), 
        500, 0, 3000 );
  }

  clusterADCtotTIB_ = fs->make<TH1I>("clusterADCtotTIB",
        "Uncorrected Cluster ADC Count for the TIB", 500, 0, 3000 );
  clusterADCtotTOB_ = fs->make<TH1I>("clusterADCtotTOB",
        "Uncorrected Cluster ADC Count for the TOB", 500, 0, 3000 );
  clusterADCtotTID_ = fs->make<TH1I>("clusterADCtotTID",
        "Uncorrected Cluster ADC Count for the TID", 500, 0, 3000 );
  clusterADCtotTEC_ = fs->make<TH1I>("clusterADCtotTEC",
        "Uncorrected Cluster ADC Count for the TEC", 500, 0, 3000 );

  centbins_ = fs->make<TH1I>("centbins","centbins", 40, 0, 39);

  widthDetId_ = fs->make<TH2F>("widthDetId","Cluster Widths vs Detector ID", 250, 1., 250., 500, 0., 66000.0); 

  clusterCountPixelStrip_ = fs->make<TH2F>("clusterCountPixelStrip","Pixel Clusters versus Strip Clusters",
                                          200, 0., 70000., 200, 0., 400000.);

  trackClusterModules_ =fs->make<TH2I>("trackClusterModules",
      "Tracks on Module versus Clusters on Module", 101,0,100,26,0,25);

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

     // Loop over tracks and make a map of the number of track hits on each module

     edm::Handle<reco::TrackCollection> tracks;
     iEvent.getByLabel(trackSrc_, tracks);
     
     std::map<int,int> tracksOnModule;
 
     for( reco::TrackCollection::const_iterator tk  = tracks->begin();
          tk != tracks->end(); ++tk )
     {
       trackingRecHit_iterator hit = tk->recHitsBegin();
       for( ; hit != tk->recHitsEnd(); ++hit )
       {
         std::map<int,int>::iterator iter = tracksOnModule.find( (*hit)->geographicalId() );
         if( iter != tracksOnModule.end() )
         {
           iter->second = iter->second + 1;
         } else {
           tracksOnModule.insert( std::pair<int,int>( (*hit)->geographicalId() , 1 ) );
         }

       }
     }  


     // Loop over clusters

     edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters->begin();
     for ( ; it != clusters->end(); ++it )
     {
       int clustersOnModule = 0;

       for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
       {
          clustersOnModule++;

          int firststrip = clus->firstStrip();
          int laststrip = firststrip + clus->amplitudes().size();
          clusterwidth_[bin/4]->Fill( laststrip-firststrip);     

          widthDetId_->Fill( laststrip-firststrip, it->detId() % 65536 );

          int totalADC = 0;
          std::vector<uint8_t>::const_iterator adc;
          for( adc = clus->amplitudes().begin(); adc != clus->amplitudes().end(); ++adc )
            totalADC += *adc;
          clusterADCtot_[bin/4]->Fill( totalADC );

          SiStripDetId sid( it->detId() );

    
          switch ( sid.subdetId() ) 
          {
            case SiStripDetId::TIB:
              clusterADCtotTIB_->Fill( totalADC );
              break;
            case SiStripDetId::TOB:
              clusterADCtotTOB_->Fill( totalADC );
              break;
            case SiStripDetId::TID:
              clusterADCtotTID_->Fill( totalADC );
              break;
            case SiStripDetId::TEC:
              clusterADCtotTEC_->Fill( totalADC );
              break;
          } 
       }

       std::map<int,int>::iterator iter = tracksOnModule.find( it->detId() );
       if( iter != tracksOnModule.end() )
       {
         trackClusterModules_->Fill( clustersOnModule, iter->second); 
       } else {
         trackClusterModules_->Fill( clustersOnModule, 0);
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

