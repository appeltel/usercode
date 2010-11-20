// -*- C++ -*-
//
// Package:    SimpleTrackAnalyzer
// Class:      SimpleTrackAnalyzer
// 
/**\class SimpleTrackAnalyzer SimpleTrackAnalyzer.cc Appeltel/GeneralAnalyzers/src/SimpleTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
// $Id$
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
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"


class SimpleTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SimpleTrackAnalyzer(const edm::ParameterSet&);
      ~SimpleTrackAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      TH1I * hit_;
      TH1F * dxyerr_;
      TH1F * pt_;
      TH1F * phi_;
      TH1F * eta_;
      TH2F * hit_dxyerr_;

      TH1I * chit_[10];
      TH1F * cdxyerr_[10];
      TH1F * cpt_[10];
      TH1F * cphi_[10];
      TH1F * ceta_[10];
      TH2F * chit_dxyerr_[10];
      TH1I* centbins_;


      edm::InputTag trackSrc_;
      edm::InputTag vertexSrc_;

      CentralityProvider * centrality_;



};

SimpleTrackAnalyzer::SimpleTrackAnalyzer(const edm::ParameterSet& iConfig):
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;

  hit_  = fs->make<TH1I>("hit", "Tracks by Number of Valid Hits",    31,  0,30);
  hit_->GetXaxis()->SetTitle("Number of Valid Hits");

  pt_ = fs->make<TH1F>("pt", "p_{T} Spectrum",100,0,20);
  pt_->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  eta_ = fs->make<TH1F>("eta", "Pseudorapidity Distribution",50,-2.5,2.5);
  eta_->GetXaxis()->SetTitle("#eta");

  phi_ = fs->make<TH1F>("phi", "Azimuthal Distribution",100,-3.15,3.15);
  phi_->GetXaxis()->SetTitle("#phi");

  dxyerr_ = fs->make<TH1F>("dxyerr", "Transverse DCA Significance",100,-6,6);
  dxyerr_->GetXaxis()->SetTitle("dxy / #sigma_{dxy}");


  hit_dxyerr_ = fs->make<TH2F>("hit_dxyerr", "Transverse DCA Significance versus Number of Valid Hits",
                               30,0,30,100,-6,6);
  hit_dxyerr_->GetYaxis()->SetTitle("dxy / #sigma_{dxy}");
  hit_dxyerr_->GetXaxis()->SetTitle("Number of Valid Hits");

  char centString[10][256] = { "0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                               "50-60%", "60-70%", "70-80%", "80-90%", "90-100%" };

  for ( int i = 0; i<10; i++ )
  {
    chit_[i]  = fs->make<TH1I>("hit",Form("Tracks by Number of Valid Hits Centrality %s",centString[i]),
                               31,  0,30);
    chit_[i]->GetXaxis()->SetTitle("Number of Valid Hits");

    cpt_[i] = fs->make<TH1F>("pt",Form("p_{T} Spectrum Centrality %s",centString[i]),100,0,20);
    cpt_[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    ceta_[i] = fs->make<TH1F>("eta",Form("Pseudorapidity Distribution Centrality %s",centString[i]),50,-2.5,2.5);
    ceta_[i]->GetXaxis()->SetTitle("#eta");

    cphi_[i] = fs->make<TH1F>("phi",Form("Azimuthal Distribution Centrality %s",centString[i]),100,-3.15,3.15);
    cphi_[i]->GetXaxis()->SetTitle("#phi");

    cdxyerr_[i] = fs->make<TH1F>("dxyerr",Form("Transverse DCA Significance Centrality %s",centString[i]),100,-6,6);
    cdxyerr_[i]->GetXaxis()->SetTitle("dxy / #sigma_{dxy}");

    chit_dxyerr_[i] = fs->make<TH2F>("hit_dxyerr", 
                  Form("Transverse DCA Significance versus Number of Valid Hits Centrality %s",centString[i]),
                               30,0,30,100,-6,6);
    chit_dxyerr_[i]->GetYaxis()->SetTitle("dxy / #sigma_{dxy}");
    chit_dxyerr_[i]->GetXaxis()->SetTitle("Number of Valid Hits");

  }

  centbins_ = fs->make<TH1I>("centbins","Centrality Distribution", 40, 0, 39);


  // safety
  centrality_ = 0;
}

SimpleTrackAnalyzer::~SimpleTrackAnalyzer()
{
}

void
SimpleTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   Handle<reco::TrackCollection> tcol;
   iEvent.getByLabel(trackSrc_, tcol);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);


   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
   //double c = centrality_->centralityValue();
   int bin = centrality_->getBin();

   centbins_->Fill(bin);

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
   if(vertex->size()>0) {
     vtxPoint=vertex->begin()->position();
     vzErr=vertex->begin()->zError();
     vxErr=vertex->begin()->xError();
     vyErr=vertex->begin()->yError();
   }



   for( reco::TrackCollection::const_iterator track = tcol->begin();
        track != tcol->end(); ++track )
   {

     double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
     dxy = track->dxy(vtxPoint);
     dz = track->dz(vtxPoint);
     dxysigma = sqrt(track->d0Error()*track->d0Error()+vxErr*vyErr);
     dzsigma = sqrt(track->dzError()*track->dzError()+vzErr*vzErr);

     hit_->Fill( track->numberOfValidHits() );
     pt_->Fill( track->pt() );
     eta_->Fill( track->eta() );
     phi_->Fill( track->phi() );
     dxyerr_->Fill( dxysigma );
     hit_dxyerr_->Fill ( track->numberOfValidHits(), dxysigma );

     chit_[bin/4]->Fill( track->numberOfValidHits() );
     cpt_[bin/4]->Fill( track->pt() );
     ceta_[bin/4]->Fill( track->eta() );
     cphi_[bin/4]->Fill( track->phi() );
     cdxyerr_[bin/4]->Fill( dxysigma );
     chit_dxyerr_[bin/4]->Fill ( track->numberOfValidHits(), dxysigma );
    

   }

}

void
SimpleTrackAnalyzer::beginJob()
{
}

void
SimpleTrackAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(SimpleTrackAnalyzer);
