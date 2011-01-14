//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


class GenParticleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenParticleAnalyzer(const edm::ParameterSet&);
      ~GenParticleAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      TH1I* centbins_;
      TH1F* chargedHadrons_[10];
      TH1F* pions_[10];
      TH1F* kaons_[10];
      TH1F* protons_[10];
      TH1F* antiprotons_[10];

      edm::InputTag genSrc_;
      double etaMax_;

      CentralityProvider * centrality_;



};

GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig):
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
etaMax_(iConfig.getParameter<double>("etaMax"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;

  char centString[10][256] = { "0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                               "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };

  for ( int i = 0; i<10; i++ )
  {
    chargedHadrons_[i] = fs->make<TH1F>(Form("chargedHadrons%d",i),
      Form("Charged Hadron Spectrum, |#eta| < %f %s",etaMax_, centString[i]),30,0.,3.); 
    pions_[i] = fs->make<TH1F>(Form("pions%d",i),
      Form("Charged Pion Spectrum, |#eta| < %f %s",etaMax_, centString[i]),30,0.,3.); 
    kaons_[i] = fs->make<TH1F>(Form("kaons%d",i),
      Form("Charged Kaon Spectrum, |#eta| < %f %s",etaMax_, centString[i]),30,0.,3.); 
    protons_[i] = fs->make<TH1F>(Form("protons%d",i),
      Form("Proton Spectrum, |#eta| < %f %s",etaMax_, centString[i]),30,0.,3.); 
    antiprotons_[i] = fs->make<TH1F>(Form("antiprotons%d",i),
      Form("Antiproton Spectrum, |#eta| < %f %s",etaMax_, centString[i]),30,0.,3.); 


  }

  centbins_ = fs->make<TH1I>("centbins","Centrality Distribution", 40, 0, 39);


  // safety
  centrality_ = 0;
}

GenParticleAnalyzer::~GenParticleAnalyzer()
{
}

void
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   Handle<reco::GenParticleCollection> gcol;
   iEvent.getByLabel(genSrc_, gcol);

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);



   for( reco::GenParticleCollection::const_iterator particle = gcol->begin();
        particle != gcol->end(); ++particle )
   {

     if ( particle->status() == 1 && 
          particle->charge() != 0 &&
          fabs(particle->eta()) <= etaMax_ )
     {
       chargedHadrons_[bin/4]->Fill( particle->pt() );

       if ( abs(particle->pdgId()) == 211 )
         pions_[bin/4]->Fill( particle->pt() );
       if ( abs(particle->pdgId()) == 321 )
         kaons_[bin/4]->Fill( particle->pt() );
       if ( particle->pdgId() == 2212 )
         protons_[bin/4]->Fill( particle->pt() );
       if ( particle->pdgId() == -2212 )
         antiprotons_[bin/4]->Fill( particle->pt() );
     }
 
    

   }

}

void
GenParticleAnalyzer::beginJob()
{
}

void
GenParticleAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(GenParticleAnalyzer);
