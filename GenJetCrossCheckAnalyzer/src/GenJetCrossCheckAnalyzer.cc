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
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>

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

      // ----------member data ---------------------------

      edm::InputTag genJetSrc_;
      double etaMin_;
      double etaMax_;
      std::vector<double> ptBins_;
      std::string pythiaProcess_;

      std::map<std::string,TH1F*> hist_;

};


//
// constructors and destructor
//

GenJetCrossCheckAnalyzer::GenJetCrossCheckAnalyzer(const edm::ParameterSet& iConfig):
genJetSrc_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
pythiaProcess_(iConfig.getParameter<std::string>("pythiaProcess"))
{
   edm::Service<TFileService> fs;

   hist_["spectrum"] = fs->make<TH1F>("spectrum",";counts;p_{T}",
                           ptBins_.size()-1, &ptBins_[0]);

   hist_["events"] = fs->make<TH1F>("events",";;events",1,0.,2.);

   std::cout << "Pyhtia process: " << pythiaProcess_ << std::endl;

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
   
   std::string mbProcess = "MB_0_to_20";
   if( pythiaProcess_ == mbProcess )
   {
     edm::Handle<GenEventInfoProduct> genEvtInfo;
     iEvent.getByLabel("generator", genEvtInfo);
     if( genEvtInfo->qScale() > 20 ) std::cout << "qScale = " << genEvtInfo->qScale() << std::endl;
     if( genEvtInfo->qScale() > 20 ) return;
   }

   // Analyze GenJets

   Handle<reco::GenJetCollection> gcol;
   iEvent.getByLabel(genJetSrc_,gcol);

   hist_["events"]->Fill(1);

   for( const auto & jet : *gcol )
   {
     if( jet.eta() <= etaMax_ && jet.eta() >= etaMin_ )
       hist_["spectrum"]->Fill(jet.pt());
   }
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
