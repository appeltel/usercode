// -*- C++ -*-
//
// Package:    GenJetCrossCheckAnalyzer
// Class:      GenJetCrossCheckAnalyzer
// 
/**\class GenJetCrossCheckAnalyzer GenJetCrossCheckAnalyzer.cc Appeltel/GenJetCrossCheckAnalyzer/src/GenJetCrossCheckAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Tue Feb 25 14:25:23 CST 2014
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "HepPDT/ParticleID.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>
#include <TH2.h>

//
// class declaration
//

class GenJetCrossCheckAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenJetCrossCheckAnalyzer(const edm::ParameterSet&);
      ~GenJetCrossCheckAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      bool isDSEvent( const edm::Event&, const edm::EventSetup& );

      // ----------member data ---------------------------

      edm::InputTag genJetSrc_;
      edm::InputTag genParticleSrc_;
      bool doCMatrix_;
      bool jetsByAbsRapidity_;
      double etaMin_;
      double etaMax_;
      double jetRadius_;
      double pthatMin_;
      double pthatMax_;
      std::vector<double> ptBins_;
      std::string pythiaProcess_;

      std::map<std::string,TH1F*> hist_;
      std::map<std::string,TH2F*> hist2D_;

};


//
// constructors and destructor
//

GenJetCrossCheckAnalyzer::GenJetCrossCheckAnalyzer(const edm::ParameterSet& iConfig):
genJetSrc_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
genParticleSrc_(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
doCMatrix_(iConfig.getParameter<bool>("doCMatrix")),
jetsByAbsRapidity_(iConfig.getParameter<bool>("jetsByAbsRapidity")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
jetRadius_(iConfig.getParameter<double>("jetRadius")),
pthatMin_(iConfig.getParameter<double>("pthatMin")),
pthatMax_(iConfig.getParameter<double>("pthatMax")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
pythiaProcess_(iConfig.getParameter<std::string>("pythiaProcess"))
{
   edm::Service<TFileService> fs;

   hist_["spectrum"] = fs->make<TH1F>("spectrum",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist_["pspectrum"] = fs->make<TH1F>("pspectrum",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist_["spectrum_DS"] = fs->make<TH1F>("spectrum_DS",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist_["pspectrum_DS"] = fs->make<TH1F>("pspectrum_DS",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist_["spectrum_fine"] = fs->make<TH1F>("spectrum_fine",";p_{T};counts",
                           1000,0.,1000.);
   hist_["pspectrum_fine"] = fs->make<TH1F>("pspectrum_fine",";p_{T};counts",
                           1000,0.2,200.);

   hist_["qscale"] = fs->make<TH1F>("qscale",";p_{T}-hat;counts",
                           1000,0.,1000.);

   hist_["events"] = fs->make<TH1F>("events",";;events",1,0.,2.);
   hist_["events_DS"] = fs->make<TH1F>("events_DS",";;double sided events",1,0.,2.);

   hist2D_["cmatrix"] = fs->make<TH2F>("cmatrix",";p_{T} Jet;p_{T} h^{#pm}",
                           1000,0.,1000.,1000,0.,200.);

   // logarthmically binned spectra, delta pt = 1%

   std::vector<double> logbins; logbins.reserve(2000); logbins.push_back(0);
   for( double bin = 0.1; bin < 1000; bin += bin*0.1)
     logbins.push_back(bin);
   hist_["spectrum_log"] = fs->make<TH1F>("spectrum_log",";counts;p_{T}",
                           logbins.size()-1, &logbins[0]);
   hist_["pspectrum_log"] = fs->make<TH1F>("pspectrum_log",";counts;p_{T}",
                           logbins.size()-1, &logbins[0]);

  
   std::cout << "Pythia process: " << pythiaProcess_ << std::endl;

}


GenJetCrossCheckAnalyzer::~GenJetCrossCheckAnalyzer()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenJetCrossCheckAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Pythia does not respect max pt-hat for MB process  
   // where the minimum is 0. In this case we need to 
   // remove events over the pt-hat maximum by hand.
   // We skip and do not count such events.
   // Here it is only coded for the MB_0_to_20 process
   edm::Handle<GenEventInfoProduct> genEvtInfo;
   iEvent.getByLabel("generator", genEvtInfo);
   
   if( pthatMin_ < 1. )
   {
     if( genEvtInfo->qScale() > pthatMax_ ) return;
   }

   // Check if the event is DS
   bool isDS = isDSEvent(iEvent,iSetup);

   Handle<reco::GenJetCollection> gcol;
   iEvent.getByLabel(genJetSrc_,gcol);

   Handle<reco::GenParticleCollection> pcol;
   iEvent.getByLabel(genParticleSrc_,pcol);

   hist_["events"]->Fill(1);
   if(isDS) hist_["events_DS"]->Fill(1);
   hist_["qscale"]->Fill(genEvtInfo->qScale());

   // genjet spectra
   for( const auto & jet : *gcol )
   {
     bool isInRange = false;

     if( jetsByAbsRapidity_ == false &&  jet.eta() <= etaMax_ && jet.eta() >= etaMin_ )
       isInRange = true;

     if( jetsByAbsRapidity_ == true &&  fabs(jet.y()) <= etaMax_ && fabs(jet.y()) >= etaMin_ )
       isInRange = true;


     if( isInRange )
     {
       hist_["spectrum"]->Fill(jet.pt());
       if(isDS) hist_["spectrum_DS"]->Fill(jet.pt());
       hist_["spectrum_fine"]->Fill(jet.pt());
       hist_["spectrum_log"]->Fill(jet.pt());
     }
   }

   // charged particle spectra and convolution matrix
   for( const auto & h : *pcol )
   {
     // skip decayed and neutral particles
     if( h.status() != 1 || h.charge() == 0  ) continue;


     if( h.eta() <= etaMax_ && h.eta() >= etaMin_ )
     {
       hist_["pspectrum"]->Fill(h.pt());
       if(isDS) hist_["pspectrum_DS"]->Fill(h.pt());
       hist_["pspectrum_fine"]->Fill(h.pt());
       hist_["pspectrum_log"]->Fill(h.pt());
     
       // associate with the highest-pt jet for which 
       // the track is found in the jet cone
       if( doCMatrix_) 
       {
         double maxPtJet = 0.;
         for( const auto & jet : *gcol )
         {
           if( jet.eta() <= etaMax_ && jet.eta() >= etaMin_ )
           {
             double dr = deltaR(jet,h);
             if( dr < jetRadius_ && jet.pt() > maxPtJet)
              maxPtJet = jet.pt();
           }
         }
         hist2D_["cmatrix"]->Fill( maxPtJet, h.pt());
       }
     }
   }

}


bool 
GenJetCrossCheckAnalyzer::isDSEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     using namespace edm;

     bool posDS = false; bool negDS = false;

     edm::ESHandle<ParticleDataTable> particleDataTable_;
     iSetup.getData(particleDataTable_);

     Handle<reco::GenParticleCollection> gcol;
     iEvent.getByLabel(genParticleSrc_, gcol);
     for( const auto & gen : * gcol )
     {
       // see if genpartice counts for DS
       HepPDT::ParticleID particleID(gen.pdgId());
       if (particleID.isValid())
       {
         ParticleData const * particleData = particleDataTable_->particle(particleID);
         if (particleData)
         { 
           double tau =  particleDataTable_->particle(particleID)->lifetime();
           if ( tau  > 1e-18 || tau == 0.0 )
           {
             if( gen.energy() > 3.0 && gen.eta() > 3.0 && gen.eta() < 5.0 ) posDS = true;
             if( gen.energy() > 3.0 && gen.eta() < -3.0 && gen.eta() > -5.0 ) negDS = true;
           }
         }
       }
     }
     if( posDS && negDS ) return true;
     else return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenJetCrossCheckAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenJetCrossCheckAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenJetCrossCheckAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenJetCrossCheckAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenJetCrossCheckAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenJetCrossCheckAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenJetCrossCheckAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetCrossCheckAnalyzer);
