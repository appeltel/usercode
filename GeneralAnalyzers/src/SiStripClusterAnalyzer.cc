// system include files
#include <memory>

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

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"


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

      edm::InputTag clusterSrc_;

      double etaMin_;
      double etaMax_;

      CentralityProvider * centrality_;
 
      TH1I* clusterwidth_[10];

      TH1I* centbins_;

};

SiStripClusterAnalyzer::SiStripClusterAnalyzer(const edm::ParameterSet& iConfig):
clusterSrc_(iConfig.getParameter<edm::InputTag>("clusterSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax"))
{

  edm::Service<TFileService> fs;

  for( int i = 0; i<10; i++)
  {
    clusterwidth_[i] = fs->make<TH1I>(Form("clusterwidth%d",i),Form("clusterwidth%d",i), 250, 1, 250 );
  }
  centbins_ = fs->make<TH1I>("centbins","centbins", 40, 0, 39);

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

   edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
   edm::InputTag clusLabel("siStripClusters");
   iEvent.getByLabel(clusLabel, clusters);

   edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters->begin();
   for ( ; it != clusters->end(); ++it )
   {
     for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
     {
        int firststrip = clus->firstStrip();
        int laststrip = firststrip + clus->amplitudes().size();
        clusterwidth_[bin/4]->Fill( laststrip-firststrip);     

     }
   }

   


}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripClusterAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripClusterAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripClusterAnalyzer);

