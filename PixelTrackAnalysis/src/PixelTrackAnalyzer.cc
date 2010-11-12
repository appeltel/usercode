// -*- C++ -*-
//
// Package:    PixelTrackAnalysis
// Class:      PixelTrackAnalysis
// 
/**\class PixelTrackAnalysis PixelTrackAnalysis.cc Appeltel/PixelTrackAnalysis/src/PixelTrackAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric A. Appelt
//         Created:  Fri Nov 12 04:59:50 EST 2010
// $Id: PixelTrackAnalyzer.cc,v 1.1 2010/11/12 10:50:53 appeltel Exp $
//
//


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


#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


//
// class declaration
//

class PixelTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PixelTrackAnalyzer(const edm::ParameterSet&);
      ~PixelTrackAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag pixelSrc_;
      edm::InputTag vertexSrc_;

      double etaCut_;

      CentralityProvider * centrality_;



      TH1F* pixelpt_; 
      TH1F* pixelphi_;
      TH1F* pixeleta_;
      TH1F* pixeld0_; 
      TH1F* pixeld0err_; 
      TH1F* pixeldz_; 
      TH1F* pixeldzerr_; 
 
      TH1F* pixelptcent_[10];

      TH1I* centbins_;

};

PixelTrackAnalyzer::PixelTrackAnalyzer(const edm::ParameterSet& iConfig):
pixelSrc_(iConfig.getParameter<edm::InputTag>("pixelSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
etaCut_(iConfig.getParameter<double>("etaCut"))
{

  edm::Service<TFileService> fs;

  pixelpt_  = fs->make<TH1F>("pixelpt", "pixelpt",    200,  0.,10.);
  pixelphi_ = fs->make<TH1F>("pixelphi", "pixelphi", 200, -3.15, 3.15);
  pixeleta_ = fs->make<TH1F>("pixeleta", "pixeleta", 200, -2.5, 2.5);
  pixeld0_ = fs->make<TH1F>("pixeld0", "pixeld0", 200, -0.25, 0.25);
  pixeld0err_ = fs->make<TH1F>("pixeld0err", "pixeld0err", 200, 0, 5);
  pixeldz_ = fs->make<TH1F>("pixeldz", "pixeldz", 200, -0.25, 0.25);
  pixeldzerr_ = fs->make<TH1F>("pixeldzerr", "pixeldzerr", 200, 0, 5);
  for( int i = 0; i<10; i++)
  {
    pixelptcent_[i] = fs->make<TH1F>(Form("pixelptcent%d",i),Form("pixelptcent%d",i), 200, 0., 10. );
  }
  centbins_ = fs->make<TH1I>("centbins","centbins", 40, 0, 39);

  // safety
  centrality_ = 0;
}


PixelTrackAnalyzer::~PixelTrackAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PixelTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<reco::Track> > tracks;
   iEvent.getByLabel(pixelSrc_, tracks);

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
     if ( fabs( tk->eta() ) <= etaCut_ )
     {
       pixelpt_->Fill( tk->pt() );
       pixeleta_->Fill( tk->eta() );
       pixelphi_->Fill( tk->phi() );

       double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
       d0 = -1.*tk->dxy(vtxPoint);
       dz = tk->dz(vtxPoint);
       d0sigma = sqrt(tk->d0Error()*tk->d0Error()+vxErr*vyErr);
       dzsigma = sqrt(tk->dzError()*tk->dzError()+vzErr*vzErr);
   
       pixeld0_->Fill( d0 );
       pixeldz_->Fill( dz );
       pixeldzerr_->Fill( fabs(dz/dzsigma) );
       pixeld0err_->Fill( fabs(d0/d0sigma) ); 

       pixelptcent_[bin/4]->Fill ( tk->pt() );

     }
  }


}


// ------------ method called once each job just before starting event loop  ------------
void 
PixelTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelTrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelTrackAnalyzer);
