// -*- C++ -*-
//
// Package:    ImpactParameterFilter
// Class:      ImpactParameterFilter
// 
/**\class ImpactParameterFilter ImpactParameterFilter.cc Appeltel/ImpactParameterFilter/src/ImpactParameterFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Mon Jun 30 13:40:49 CDT 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/HeavyIon.h"

//
// class declaration
//

class ImpactParameterFilter : public edm::EDFilter {
   public:
      explicit ImpactParameterFilter(const edm::ParameterSet&);
      ~ImpactParameterFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      std::vector<std::string> hepmcSrc_;
      double bMin_;
      double bMax_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ImpactParameterFilter::ImpactParameterFilter(const edm::ParameterSet& iConfig):
hepmcSrc_(iConfig.getParameter<std::vector<std::string> >("generators")),
bMin_(iConfig.getParameter<double>("bMin")),
bMax_(iConfig.getParameter<double>("bMax"))
{
   //now do what ever initialization is needed

}


ImpactParameterFilter::~ImpactParameterFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ImpactParameterFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   bool inRange = false;

   for(size_t ihep = 0; ihep < hepmcSrc_.size(); ++ihep)
   {
     edm::Handle<edm::HepMCProduct> hepmc;
     iEvent.getByLabel(hepmcSrc_[ihep],hepmc);
     const HepMC::HeavyIon* hi = hepmc->GetEvent()->heavy_ion();
     double b = hi->impact_parameter();
     if( b < bMax_ && b > bMin_ ) inRange = true;
   } 

   return inRange;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ImpactParameterFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ImpactParameterFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ImpactParameterFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ImpactParameterFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ImpactParameterFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ImpactParameterFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ImpactParameterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ImpactParameterFilter);
