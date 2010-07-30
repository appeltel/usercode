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

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"
#include "CLHEP/Random/RandomEngine.h"


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
      inline float integrateADC ( std::vector<int16_t>& sample, double offset);

      void processRaw(const edm::InputTag&, const edm::DetSetVector<SiStripRawDigi>&, std::vector<edm::DetSet<SiStripDigi> >& );
      std::auto_ptr<SiStripRawProcessingAlgorithms> algorithms;

      edm::InputTag inputTag;
      edm::InputTag inputTagNoise;
      edm::InputTag inputTagSignal;
      bool doNoiseAndSignal;


      CLHEP::HepRandomEngine* rndEngine;
      CLHEP::RandGaussQ* gaussDistribution;
      GaussianTailNoiseGenerator* genNoise;


      float cmnTIBRMS_;
      float cmnTOBRMS_;
      float cmnTIDRMS_;
      float cmnTECRMS_;


      TH1F* allMedians_;
      TH1F* all25ths_;
      TH1F* TIBMedians_;
      TH1F* TIB25ths_;
      TH1F* TIBGens_;
      TH1F* TOBMedians_;
      TH1F* TOB25ths_;
      TH1F* TOBGens_;
      TH1F* TIDMedians_;
      TH1F* TID25ths_;
      TH1F* TIDGens_;
      TH1F* TECMedians_;
      TH1F* TEC25ths_;
      TH1F* TECGens_;
      TH2F* medvs25_;

      // These are displays of ten interesting APV25s, with the median and 25th in the 129th and 130th points.

      TGraph* stripCount_[10];
      TGraph* stripSignalCount_[10];

      // Histograms of integrated ADC  on each APV

      TH1F* totalADCMedian_;
      TH1F* totalADCPer25_;

      // for comparison to actual integrated ADC from signal only with perfect baseline

      TH1F* totalADCSignal_;
      TH1F* totalADCMedianDiff_;
      TH1F* totalADCPer25Diff_;



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
     inputTag(iConfig.getParameter<edm::InputTag>("RawDigiProducersList")),
     inputTagNoise(iConfig.getParameter<edm::InputTag>("RawDigiProducersListNoise")),
     inputTagSignal(iConfig.getParameter<edm::InputTag>("RawDigiProducersListSignal")),
     doNoiseAndSignal(iConfig.getParameter<bool>("doNoiseAndSignal"))
{

  edm::Service<TFileService> fs;
  allMedians_=fs->make<TH1F>("allMedians","Median ADC count of APV25", 500, 0., 200.);
  all25ths_=fs->make<TH1F>("all25ths","25th percentile ADC count of APV25", 500, 0., 200.);
  medvs25_=fs->make<TH2F>("medvs25","Median ADC count versus 25th percentile",500,0.,500.,500,0.,200.);

  TIBMedians_=fs->make<TH1F>("TIBMedians","Median ADC count of APV25 for the TIB", 500, 0., 200.);
  TIB25ths_=fs->make<TH1F>("TIB25ths","25th percentile ADC count of APV25 for the TIB", 500, 0., 200.);
  TIBGens_=fs->make<TH1F>("TIBGens","Expected CMN Distribution for the TIB", 500, 0., 200.);

  TOBMedians_=fs->make<TH1F>("TOBMedians","Median ADC count of APV25 for the TOB", 500, 0., 200.);
  TOB25ths_=fs->make<TH1F>("TOB25ths","25th percentile ADC count of APV25 for the TOB", 500, 0., 200.);
  TOBGens_=fs->make<TH1F>("TOBGens","Expected CMN Distribution for the TOB", 500, 0., 200.);

  TIDMedians_=fs->make<TH1F>("TIDMedians","Median ADC count of APV25 for the TID", 500, 0., 200.);
  TID25ths_=fs->make<TH1F>("TID25ths","25th percentile ADC count of APV25 for the TID", 500, 0., 200.);
  TIDGens_=fs->make<TH1F>("TIDGens","Expected CMN Distribution for the TID", 500, 0., 200.);

  TECMedians_=fs->make<TH1F>("TECMedians","Median ADC count of APV25 for the TEC", 500, 0., 200.);
  TEC25ths_=fs->make<TH1F>("TEC25ths","25th percentile ADC count of APV25 for the TEC", 500, 0., 200.);
  TECGens_=fs->make<TH1F>("TECGens","Expected CMN Distribution for the TEC", 500, 0., 200.);

  totalADCMedian_=fs->make<TH1F>("totalADCMedian","Integrated ADC count for each APV25 after Median CMN subtraction",500,0.,10000.);
  totalADCPer25_=fs->make<TH1F>("totalADCPer25","Integrated ADC count for each APV25 after 25th Percentile CMN subtraction",500,0.,10000.);
  if(doNoiseAndSignal)
  {
    totalADCSignal_=fs->make<TH1F>("totalADCSignal","Integrated ADC count for each APV25 from true SimHits only",500,0.,10000.);
    totalADCMedianDiff_=fs->make<TH1F>("totalADCMedianDiff",
      "Difference between median subtracted integrated ADC count and ADC count from signal",500,-3000.,2000.);
    totalADCPer25Diff_=fs->make<TH1F>("totalADCPer25Diff",
      "Difference between 25th percentile subtracted integrated ADC count and ADC count from signal",500,-3000.,2000.);
  }    
  

  for( int i = 0; i<10; i++)
  {
    stripCount_[i] = fs->make<TGraph>(130);
    stripCount_[i]->SetName(Form("stripCount%d",i));
    if ( doNoiseAndSignal )
    {
      stripSignalCount_[i] = fs->make<TGraph>(128);
      stripSignalCount_[i]->SetName(Form("stripSignalCount%d",i));
    }
  }

  firstEvent = true;

  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << "SiStripCMNAnalyzer requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }
  
  rndEngine       = &(rng->getEngine());
  genNoise = new GaussianTailNoiseGenerator(*rndEngine);
  gaussDistribution = new CLHEP::RandGaussQ(*rndEngine);

  cmnTIBRMS_ = 5.92;
  cmnTOBRMS_ = 1.08;
  cmnTIDRMS_ = 3.08;
  cmnTECRMS_ = 2.44;

}


SiStripCMNAnalyzer::~SiStripCMNAnalyzer()
{

  delete genNoise;
  delete gaussDistribution;
 
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
   
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > inputNoise;
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > inputSignal;
   if ( doNoiseAndSignal )
   {
     iEvent.getByLabel(inputTagNoise,inputNoise);
     iEvent.getByLabel(inputTagSignal,inputSignal);
   }

   float maxDiff = 0;
   int currentGraph = 0;

   for ( edm::DetSetVector<SiStripRawDigi>::const_iterator 
         rawDigis = input->begin(); rawDigis != input->end(); rawDigis++) 
   {
     
    edm::DetSetVector<SiStripRawDigi>::const_iterator  rawSignalDigis;

    if ( doNoiseAndSignal )
    {
      // need to replace this with an Association Map!
      rawSignalDigis = inputSignal->begin();
      for (rawSignalDigis = inputSignal->begin(); rawSignalDigis != inputSignal->end(); rawSignalDigis++)
        if( rawDigis->detId() == rawSignalDigis->detId() ) break;
    }  


    SiStripDetId sid( rawDigis->detId() );

    std::vector<int16_t> processedRawDigis(rawDigis->size());
    algorithms->subtractorPed->subtract( *rawDigis, processedRawDigis);
 
    std::vector<int16_t> signalDigis(rawDigis->size());
    if ( doNoiseAndSignal )
    {
      if ( rawSignalDigis != inputSignal->end() && rawSignalDigis->size() == rawDigis->size() )
      {
        std::vector<int16_t>::iterator sD = signalDigis.begin();
        edm::DetSet<SiStripRawDigi>::const_iterator rSD = rawSignalDigis->begin();
        while ( rSD != rawSignalDigis->end() )
        {
          *sD = rSD->adc();
          ++sD;
          ++rSD;
        }
      } 
      else 
      {
        for( std::vector<int16_t>::iterator sD = signalDigis.begin(); sD != signalDigis.end(); sD++ )
          *sD = 0;
      }
    }


    // MEDIAN ALGORITHM

    std::vector<int16_t> tmp;  tmp.reserve(128);  
    std::vector<int16_t> tmpSignal;  tmpSignal.reserve(128);  
    std::vector<int16_t>::iterator  
        strip( processedRawDigis.begin() ), 
        signal(signalDigis.begin() ),
        end(   processedRawDigis.end()   ),
        endAPV, 
        endAPVSignal;
  
    while( strip < end ) {
      endAPV = strip+128;
      endAPVSignal = signal+128; 
      tmp.clear(); tmpSignal.clear();
      tmp.insert(tmp.end(),strip,endAPV);
      tmpSignal.insert(tmpSignal.end(),signal,endAPVSignal);
      const float offset = median(tmp);
      const float per25 = percentile(tmp,25.);

      // Fill offset plots for all algorthims for all detectors and for TIB, TOB, TID, TEC

      allMedians_->Fill( offset );
      all25ths_->Fill( per25 );

      switch ( sid.subdetId() ) 
      {
        case SiStripDetId::TIB:
          TIBMedians_->Fill(offset);
          TIB25ths_->Fill(per25);
          TIBGens_->Fill( gaussDistribution->fire(128.,cmnTIBRMS_));
          break;
        case SiStripDetId::TOB:
          TOBMedians_->Fill(offset);
          TOB25ths_->Fill(per25);
          TOBGens_->Fill( gaussDistribution->fire(128.,cmnTOBRMS_));
          break;
        case SiStripDetId::TID:
          TIDMedians_->Fill(offset);
          TID25ths_->Fill(per25);
          TIDGens_->Fill( gaussDistribution->fire(128.,cmnTIDRMS_));
          break;
        case SiStripDetId::TEC:
          TECMedians_->Fill(offset);
          TEC25ths_->Fill(per25);
          TECGens_->Fill( gaussDistribution->fire(128.,cmnTECRMS_));
          break;
      }      
      medvs25_->Fill( offset, per25);

      // Integrated ADC after baseline subtraction plots
      float intMed = integrateADC(tmp, offset );
      float intPer25 = integrateADC(tmp, per25 );

      totalADCMedian_->Fill( intMed );
      totalADCPer25_->Fill( intPer25 );
      if ( doNoiseAndSignal )
      {
         float intSignal = integrateADC(tmpSignal, 0. );
         totalADCSignal_->Fill( intSignal );
         totalADCMedianDiff_->Fill( intMed - intSignal );
         totalADCPer25Diff_->Fill( intPer25 - intSignal );
      }

      // make plots of the 10 wonkiest strips in the first event

      if ( offset-per25 > maxDiff && firstEvent)
      {
        maxDiff = offset-per25;
        for( int i =0; i<128; i++)
        {
          stripCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(*(strip+i)));
          if(doNoiseAndSignal) stripSignalCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(*(signal+i)));
        }
        stripCount_[currentGraph]->SetPoint(128,150.,(double)offset);
        stripCount_[currentGraph]->SetPoint(129,150.,(double)per25);
        currentGraph++;
        currentGraph %= 10;
      }


      strip += 128;
      signal += 128;
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

inline float SiStripCMNAnalyzer::integrateADC ( std::vector<int16_t>& sample, double offset ) {
  float sum = 0.;
  for( std::vector<int16_t>::iterator it = sample.begin(); it != sample.end(); ++it)
  { 
    sum += *it-offset > 0. ? *it-offset : 0. ;
  }
  return sum;
}



//define this as a plug-in
DEFINE_FWK_MODULE(SiStripCMNAnalyzer);
