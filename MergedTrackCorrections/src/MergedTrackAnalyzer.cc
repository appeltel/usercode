//
// Original Author:  Eric A. Appelt
//         Created:  Fri Nov 12 04:59:50 EST 2010


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

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "Appeltel/MergedTrackCorrections/interface/MergedTrackCorrector.h"

//
// class declaration
//

class MergedTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MergedTrackAnalyzer(const edm::ParameterSet&);
      ~MergedTrackAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag trackSrc_;
      edm::InputTag vertexSrc_;
      edm::InputTag genSrc_;

      bool doGenerator_;

      double etaCut_;
      double ptMin_;

      CentralityProvider * centrality_;

      MergedTrackCorrector * corrector_;

      TH1F* trackpt_; 
      TH1F* trackphi_;
      TH1F* tracketa_;
      TH1F* trackd0_; 
      TH1F* trackd0err_; 
      TH1F* trackdz_; 
      TH1F* trackdzerr_; 
      TH1F* trackchi2_;
      TH1F* tracknhit_;
 
      TH1F* trackptcent_[10];
      TH1F* trackphicent_[10];
      TH1F* tracketacent_[10];
      TH1F* trackd0cent_[10]; 
      TH1F* trackd0errcent_[10]; 
      TH1F* trackdzcent_[10]; 
      TH1F* trackdzerrcent_[10]; 
      TH1F* trackchi2cent_[10];
      TH1F* tracknhitcent_[10];


      TH1F* genptcent_[10];
      TH1F* genphicent_[10];
      TH1F* genetacent_[10];
      TH1I* centbins_;

};

MergedTrackAnalyzer::MergedTrackAnalyzer(const edm::ParameterSet& iConfig):
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
doGenerator_(iConfig.getParameter<bool>("doGenerator")),
etaCut_(iConfig.getParameter<double>("etaCut")),
ptMin_(iConfig.getParameter<double>("ptMin"))
{

  edm::Service<TFileService> fs;

  corrector_  = new MergedTrackCorrector();

  corrector_->loadEffic();
  corrector_->loadFakes();

  char centStrings[10][256] = { "0-10%", "10-20%", "20-30%", "30-40%",
       "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };

  trackpt_  = fs->make<TH1F>("trackpt", "Merged Track p_{T} Spectrum",    200,  0.,10.);
  trackphi_ = fs->make<TH1F>("trackphi", "Merged Track #phi Distribution", 200, -3.15, 3.15);
  tracketa_ = fs->make<TH1F>("tracketa", "Merged Track #eta Distribution", 200, -2.5, 2.5);
  trackd0_ = fs->make<TH1F>("trackd0", "Merged Track Transverse DCA Distribution", 200, -0.25, 0.25);
  trackd0err_ = fs->make<TH1F>("trackd0err", "Merged Track Transverse DCA Significance", 200, 0, 5);
  trackdz_ = fs->make<TH1F>("trackdz", "Merged Track Longitudinal DCA Distribution", 200, -0.25, 0.25);
  trackdzerr_ = fs->make<TH1F>("trackdzerr", "Merged Track Longitudinal DCA Significance", 200, 0, 5);
  trackchi2_ = fs->make<TH1F>("trackchi2", "Merged Track #chi^2 / n.d.o.f. Distribution",200,0.,100.);
  tracknhit_ = fs->make<TH1F>("tracknhit", "Merged Track Valid Hit Distribution",30.,0.,30.);


  for( int i = 0; i<10; i++)
  {
    trackptcent_[i] = fs->make<TH1F>(Form("trackptcent%d",i),
        Form("Merged Track p_{T} Spectrum Centrality Range %s",centStrings[i]), 200, 0., 10. );
    trackphicent_[i] = fs->make<TH1F>(Form("trackphicent%d",i),
        Form("Merged Track #phi Distribution Centrality Range %s",centStrings[i]), 200, -3.15, 3.15 );
    tracketacent_[i] = fs->make<TH1F>(Form("tracketacent%d",i),
        Form("Merged Track #eta Distribution Centrality Range %s",centStrings[i]), 200, -2.5, 2.5 );
    trackd0cent_[i] = fs->make<TH1F>(Form("trackd0cent%d",i),
        Form("Merged Track Transverse DCA  Distribution Centrality Range %s",centStrings[i]), 200, -0.25, 0.25 );
    trackd0errcent_[i] = fs->make<TH1F>(Form("trackd0errcent%d",i),
        Form("Merged Track Transverse DCA  Significance Centrality Range %s",centStrings[i]), 200, 0, 5 );
    trackdzcent_[i] = fs->make<TH1F>(Form("trackdzcent%d",i),
        Form("Merged Track Longitudinal DCA  Distribution Centrality Range %s",centStrings[i]), 200, -0.25, 0.25 );
    trackdzerrcent_[i] = fs->make<TH1F>(Form("trackdzerrcent%d",i),
        Form("Merged Track Longitudinal DCA  Significance Centrality Range %s",centStrings[i]), 200, 0, 5 );
    trackchi2cent_[i] = fs->make<TH1F>(Form("trackchi2cent%d",i), 
        Form("Merged Track #chi^2 / n.d.o.f. Distribution Centrality Range %s",centStrings[i]),200,0.,100.);
    tracknhitcent_[i] = fs->make<TH1F>(Form("tracknhitcent%d",i), 
        Form("Merged Track Valid Hit Distribution Centrality Range %s",centStrings[i]),30.,0.,30.);
    genptcent_[i] = fs->make<TH1F>(Form("genptcent%d",i),
        Form("Charged GenParticle p_{T} Spectrum Centrality Range %s",centStrings[i]), 200, 0., 10. );
    genphicent_[i] = fs->make<TH1F>(Form("genphicent%d",i),
        Form("Charged GenParticle #phi Distribution Centrality Range %s",centStrings[i]), 200, -3.15, 3.15 );
    genetacent_[i] = fs->make<TH1F>(Form("genetacent%d",i),
        Form("Charged GenParticle #eta Distribution Centrality Range %s",centStrings[i]), 200, -2.5, 2.5 );

  }
  centbins_ = fs->make<TH1I>("centbins","Centrality Distribution", 40, 0, 39);

  // safety
  centrality_ = 0;
}


MergedTrackAnalyzer::~MergedTrackAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MergedTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<reco::Track> > tracks;
   iEvent.getByLabel(trackSrc_, tracks);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   
   // Get Centrality bin for the event

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);

   // find the vertex point and error

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
   if(vertex->size()>0) {
     vtxPoint=vertex->begin()->position();
     vzErr=vertex->begin()->zError();
     vxErr=vertex->begin()->xError();
     vyErr=vertex->begin()->yError();
   }


   // loop over tracks, fill histos

   std::vector<reco::Track>::const_iterator tk = tracks->begin();
   for ( ; tk != tracks->end(); ++tk )
   {
     if ( fabs( tk->eta() ) <= etaCut_ && tk->pt() >= ptMin_ )
     {

       double weight = corrector_->weight(tk->pt(), tk->eta(), bin);

       trackpt_->Fill( tk->pt(), weight );
       tracketa_->Fill( tk->eta(), weight );
       trackphi_->Fill( tk->phi(), weight );
       trackchi2_->Fill( tk->normalizedChi2(), weight );
       tracknhit_->Fill( tk->numberOfValidHits(), weight );

       trackptcent_[bin/4]->Fill( tk->pt(), weight );
       tracketacent_[bin/4]->Fill( tk->eta(), weight );
       trackphicent_[bin/4]->Fill( tk->phi(), weight );
       trackchi2cent_[bin/4]->Fill( tk->normalizedChi2(), weight );
       tracknhitcent_[bin/4]->Fill( tk->numberOfValidHits(), weight );

       double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
       d0 = -1.*tk->dxy(vtxPoint);
       dz = tk->dz(vtxPoint);
       d0sigma = sqrt(tk->d0Error()*tk->d0Error()+vxErr*vyErr);
       dzsigma = sqrt(tk->dzError()*tk->dzError()+vzErr*vzErr);
   
       trackd0_->Fill( d0, weight );
       trackdz_->Fill( dz, weight );
       trackdzerr_->Fill( fabs(dz/dzsigma), weight );
       trackd0err_->Fill( fabs(d0/d0sigma), weight ); 

       trackd0cent_[bin/4]->Fill( d0, weight );
       trackdzcent_[bin/4]->Fill( dz, weight );
       trackdzerrcent_[bin/4]->Fill( fabs(dz/dzsigma), weight );
       trackd0errcent_[bin/4]->Fill( fabs(d0/d0sigma), weight ); 

       trackptcent_[bin/4]->Fill ( tk->pt(), weight );

     }
  }

  if ( doGenerator_ )
  {

    Handle<reco::GenParticleCollection > gens;
    iEvent.getByLabel(genSrc_, gens);

    reco::GenParticleCollection::const_iterator q;
    for( q = gens->begin(); q != gens->end(); ++q )
    {
       if ( q->charge() != 0 && q->status() == 1 )
       {
         genptcent_[bin/4]->Fill( q->pt() );
         genetacent_[bin/4]->Fill( q->eta() );
         genphicent_[bin/4]->Fill( q->phi() );
       }

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MergedTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MergedTrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MergedTrackAnalyzer);

