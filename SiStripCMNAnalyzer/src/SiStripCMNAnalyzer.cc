// -*- C++ -*-
//
// Package:    SiStripCMNAnalyzer
// Class:      SiStripCMNAnalyzer
// 
/**\class SiStripCMNAnalyzer SiStripCMNAnalyzer.cc Appeltel/SiStripCMNAnalyzer/src/SiStripCMNAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Mon Jul 26 10:37:24 CDT 2010
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingFactory.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingAlgorithms.h"

//#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
//#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
//#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
//#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"



#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TGraph.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

class SiStripCMNAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripCMNAnalyzer(const edm::ParameterSet&);
      ~SiStripCMNAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      inline float median( std::vector<int16_t>& sample);
      inline float percentile ( std::vector<int16_t>& sample, double pct);

      void processRaw(const edm::InputTag&, const edm::DetSetVector<SiStripRawDigi>&, std::vector<edm::DetSet<SiStripDigi> >& );
      std::auto_ptr<SiStripRawProcessingAlgorithms> algorithms;
      edm::InputTag inputTag;

      TH1F* allMedians_;
      TH1F* all25ths_;
      TH1F* TIBMedians_;
      TH1F* TIB25ths_;
      TH1F* TOBMedians_;
      TH1F* TOB25ths_;
      TH1F* TIDMedians_;
      TH1F* TID25ths_;
      TH1F* TECMedians_;
      TH1F* TEC25ths_;
      TH2F* medvs25_;

      // These are displays of ten interesting APV25s, with the median and 25th in the 129th and 130th points.

      TGraph* stripCount_[10];

      bool firstEvent;

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
SiStripCMNAnalyzer::SiStripCMNAnalyzer(const edm::ParameterSet& iConfig)
  :  algorithms(SiStripRawProcessingFactory::create(iConfig.getParameter<edm::ParameterSet>("Algorithms"))),
     inputTag(iConfig.getParameter<edm::InputTag>("RawDigiProducersList"))
{

  edm::Service<TFileService> fs;
  allMedians_=fs->make<TH1F>("allMedians","Median ADC count of APV25", 100, 0., 500.);
  all25ths_=fs->make<TH1F>("all25ths","25th percentile ADC count of APV25", 100, 0., 500.);
  TIBMedians_=fs->make<TH1F>("TIBMedians","Median ADC count of APV25 for the TIB", 100, 0., 500.);
  TIB25ths_=fs->make<TH1F>("TIB25ths","25th percentile ADC count of APV25 for the TIB", 100, 0., 500.);
  TOBMedians_=fs->make<TH1F>("TOBMedians","Median ADC count of APV25 for the TOB", 100, 0., 500.);
  TOB25ths_=fs->make<TH1F>("TOB25ths","25th percentile ADC count of APV25 for the TOB", 100, 0., 500.);
  TIDMedians_=fs->make<TH1F>("TIDMedians","Median ADC count of APV25 for the TID", 100, 0., 500.);
  TID25ths_=fs->make<TH1F>("TID25ths","25th percentile ADC count of APV25 for the TID", 100, 0., 500.);
  TECMedians_=fs->make<TH1F>("TECMedians","Median ADC count of APV25 for the TEC", 100, 0., 500.);
  TEC25ths_=fs->make<TH1F>("TEC25ths","25th percentile ADC count of APV25 for the TEC", 100, 0., 500.);
  medvs25_=fs->make<TH2F>("medvs25","Median ADC count versus 25th percentile",50,0.,500.,50,0.,500.);

  for( int i = 0; i<10; i++)
  {
    stripCount_[i] = fs->make<TGraph>(130);
    stripCount_[i]->SetName(Form("stripCount%d",i));

  }

  firstEvent = true;
}


SiStripCMNAnalyzer::~SiStripCMNAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiStripCMNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;




   algorithms->initialize(iSetup);
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > input;
   iEvent.getByLabel(inputTag,input);

   float maxDiff = 0;
   int currentGraph = 0;

   for ( edm::DetSetVector<SiStripRawDigi>::const_iterator 
         rawDigis = input->begin(); rawDigis != input->end(); rawDigis++) 
   {
    
    SiStripDetId sid( rawDigis->detId() );

    std::vector<int16_t> processedRawDigis(rawDigis->size());
    algorithms->subtractorPed->subtract( *rawDigis, processedRawDigis);
 
    // MEDIAN ALGORITHM

    std::vector<int16_t> tmp;  tmp.reserve(128);  
    std::vector<int16_t>::iterator  
        strip( processedRawDigis.begin() ), 
        end(   processedRawDigis.end()   ),
        endAPV;
  
    while( strip < end ) {
      endAPV = strip+128; tmp.clear();
      tmp.insert(tmp.end(),strip,endAPV);
      const float offset = median(tmp);
      const float per25 = percentile(tmp,25.);

      allMedians_->Fill( offset );
      all25ths_->Fill( per25 );

      switch ( sid.subdetId() ) 
      {
        case SiStripDetId::TIB:
          TIBMedians_->Fill(offset);
          TIB25ths_->Fill(per25);
          break;
        case SiStripDetId::TOB:
          TOBMedians_->Fill(offset);
          TOB25ths_->Fill(per25);
          break;
        case SiStripDetId::TID:
          TIDMedians_->Fill(offset);
          TID25ths_->Fill(per25);
          break;
        case SiStripDetId::TEC:
          TECMedians_->Fill(offset);
          TEC25ths_->Fill(per25);
          break;
      }      
      medvs25_->Fill( offset, per25);

      // make plots of the 10 wonkiest strips in the first event

      if ( offset-per25 > maxDiff && firstEvent)
      {
        maxDiff = offset-per25;
        for( int i =0; i<128; i++)
        {
          stripCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(*(strip+i)));
        }
        stripCount_[currentGraph]->SetPoint(128,150.,(double)offset);
        stripCount_[currentGraph]->SetPoint(139,150.,(double)per25);
        currentGraph++;
        currentGraph %= 10;
      }


      strip += 128;
    }
    

    
  }

  firstEvent = false;

}


// ------------ method called once each job just before starting event loop  ------------
void 
SiStripCMNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripCMNAnalyzer::endJob() {
}


inline float SiStripCMNAnalyzer::median( std::vector<int16_t>& sample) {
  std::vector<int16_t>::iterator mid = sample.begin() + sample.size()/2;
  std::nth_element(sample.begin(), mid, sample.end());
  if( sample.size() & 1 ) //odd size
    return *mid;
  return ( *std::max_element(sample.begin(), mid) + *mid ) / 2.;
}

inline float SiStripCMNAnalyzer::percentile( std::vector<int16_t>& sample, double pct) {
  std::vector<int16_t>::iterator mid = sample.begin() + int(sample.size()*pct/100.0);
  std::nth_element(sample.begin(), mid, sample.end());
  return *mid;
} 


//define this as a plug-in
DEFINE_FWK_MODULE(SiStripCMNAnalyzer);
