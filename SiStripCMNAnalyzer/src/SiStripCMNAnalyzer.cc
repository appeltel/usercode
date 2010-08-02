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
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingFactory.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingAlgorithms.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithm.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithmFactory.h"


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

   typedef   edmNew::DetSetVector<SiStripCluster> ClusterCollection;
   typedef   edmNew::DetSet<SiStripCluster>::const_iterator ClusIter;

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

      void subtractMedian( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results);
      void subtractPer25( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results);

      int digiCount( edm::DetSet<SiStripDigi>& module, int APV );

      void processRaw(const edm::InputTag&, const edm::DetSetVector<SiStripRawDigi>&, std::vector<edm::DetSet<SiStripDigi> >& );

      // As a dumb hack, clusters for the gallery plot are stored in a TGraph with  
      // (x,y) = (first strip, last strip) or (0,0) if there is no cluster.

      void makeClusterGraph( edmNew::DetSetVector<SiStripCluster>& clusters, TGraph * graph, int APV );

      std::auto_ptr<SiStripRawProcessingAlgorithms> algorithms;
      std::auto_ptr<StripClusterizerAlgorithm> clusterAlgorithm;

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

      // Occupancy Comparison

      TH2F* medianOccupancy_;
      TH2F* medianTrueOccupancy_;
      TH2F* per25Occupancy_;
      TH2F* per25TrueOccupancy_;
 
      //
      // For each gallery there  are displays of ten interesting APV25s, 
      // with the median and 25th in the 129th and 130th points.

      TGraph* galACount_[10];
      TGraph* galASignalCount_[10];
      TGraph* galBCount_[10];
      TGraph* galBSignalCount_[10];
      TGraph* galDCount_[10];
      TGraph* galDSignalCount_[10];
      TGraph* galECount_[10];
      TGraph* galESignalCount_[10];
  
      // These graphs are really just storage for clusters for the APV25s graphed above.
      // This is a stupid way to do things and they should be replaced by nTuples
      // or something. I am doing this as a path of least resistance given my current knowledge of ROOT.

      TGraph* galAMedianClusters_[10];
      TGraph* galAPer25Clusters_[10]; 
      TGraph* galBMedianClusters_[10];
      TGraph* galBPer25Clusters_[10]; 
      TGraph* galDMedianClusters_[10];
      TGraph* galDPer25Clusters_[10]; 
      TGraph* galEMedianClusters_[10];
      TGraph* galEPer25Clusters_[10]; 
 
      // Histograms of integrated ADC  on each APV

      TH1F* totalADCMedian_;
      TH1F* totalADCPer25_;

      // for comparison to actual integrated ADC from signal only with perfect baseline

      TH1F* totalADCSignal_;
      TH1F* totalADCMedianDiff_;
      TH1F* totalADCPer25Diff_;

      // Cluster Comparison
 
      TH1F* medianTotalClus_;
      TH1F* per25TotalClus_;
      TH1F* signalTotalClus_;


      int galAcount;
      int galBcount;
      int galDcount;
      int galEcount;

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
     clusterAlgorithm( StripClusterizerAlgorithmFactory::create(iConfig.getParameter<edm::ParameterSet>("Clusterizer")) ), 
     inputTag(iConfig.getParameter<edm::InputTag>("RawDigiProducersList")),
     inputTagNoise(iConfig.getParameter<edm::InputTag>("RawDigiProducersListNoise")),
     inputTagSignal(iConfig.getParameter<edm::InputTag>("RawDigiProducersListSignal")),
     doNoiseAndSignal(iConfig.getParameter<bool>("doNoiseAndSignal"))
{

  edm::Service<TFileService> fs;
  allMedians_=fs->make<TH1F>("allMedians","Median ADC count of APV25", 500, 0., 200.);
  all25ths_=fs->make<TH1F>("all25ths","25th percentile ADC count of APV25", 500, 0., 200.);
  medvs25_=fs->make<TH2F>("medvs25","Median ADC count versus 25th percentile",500,0.,500.,500,0.,200.);

  medianOccupancy_=fs->make<TH2F>("medianOccupancy","Median ADC count of APV25 versus Reconstructed # of Unsuppressed Strips",
                                  128.,0.,128.,500,0.,500.);
  per25Occupancy_=fs->make<TH2F>("per25Occupancy","25th Percentile ADC count of APV25 versus Reconstructed # of Unsuppressed Strips",
                                  128.,0.,128.,500,0.,500.);
  if(doNoiseAndSignal)
  {
    medianTrueOccupancy_=fs->make<TH2F>("medianTrueOccupancy","Median ADC count of APV25 versus MC Correct # of Unsuppressed Strips",
                                    128.,0.,128.,500,0.,500.);
    per25TrueOccupancy_=fs->make<TH2F>("per25TrueOccupancy","25th Percentile ADC count of APV25 versus MC Correct # of Unsuppressed Strips",
                                    128.,0.,128.,500,0.,500.);
  }

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
  
  medianTotalClus_=fs->make<TH1F>("medianTotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);
  per25TotalClus_=fs->make<TH1F>("per25TotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);
  if(doNoiseAndSignal)
    signalTotalClus_=fs->make<TH1F>("signalTotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);


  for( int i = 0; i<10; i++)
  {
    galACount_[i] = fs->make<TGraph>(130);
    galACount_[i]->SetName(Form("galACount%d",i));
    galAMedianClusters_[i] = fs->make<TGraph>(128);
    galAMedianClusters_[i]->SetName(Form("galAMedianClusters%d",i));
    galAPer25Clusters_[i] = fs->make<TGraph>(128);
    galAPer25Clusters_[i]->SetName(Form("galAPer25Clusters%d",i));
    galBCount_[i] = fs->make<TGraph>(130);
    galBCount_[i]->SetName(Form("galBCount%d",i));
    galBMedianClusters_[i] = fs->make<TGraph>(128);
    galBMedianClusters_[i]->SetName(Form("galBMedianClusters%d",i));
    galBPer25Clusters_[i] = fs->make<TGraph>(128);
    galBPer25Clusters_[i]->SetName(Form("galBPer25Clusters%d",i));
    galDCount_[i] = fs->make<TGraph>(130);
    galDCount_[i]->SetName(Form("galDCount%d",i));
    galDMedianClusters_[i] = fs->make<TGraph>(128);
    galDMedianClusters_[i]->SetName(Form("galDMedianClusters%d",i));
    galDPer25Clusters_[i] = fs->make<TGraph>(128);
    galDPer25Clusters_[i]->SetName(Form("galDPer25Clusters%d",i));
    galECount_[i] = fs->make<TGraph>(130);
    galECount_[i]->SetName(Form("galECount%d",i));
    galEMedianClusters_[i] = fs->make<TGraph>(128);
    galEMedianClusters_[i]->SetName(Form("galEMedianClusters%d",i));
    galEPer25Clusters_[i] = fs->make<TGraph>(128);
    galEPer25Clusters_[i]->SetName(Form("galEPer25Clusters%d",i));
    if ( doNoiseAndSignal )
    {
      galASignalCount_[i] = fs->make<TGraph>(128);
      galASignalCount_[i]->SetName(Form("galASignalCount%d",i));
      galBSignalCount_[i] = fs->make<TGraph>(128);
      galBSignalCount_[i]->SetName(Form("galBSignalCount%d",i));
      galDSignalCount_[i] = fs->make<TGraph>(128);
      galDSignalCount_[i]->SetName(Form("galDSignalCount%d",i));
      galESignalCount_[i] = fs->make<TGraph>(128);
      galESignalCount_[i]->SetName(Form("galESignalCount%d",i));
    }
  }

  galAcount = 0;
  galBcount = 0;
  galDcount = 0;
  galEcount = 0;

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


void
SiStripCMNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // GENERAL INITIALIZATION

   algorithms->initialize(iSetup);
   clusterAlgorithm->initialize(iSetup);
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > input;
   iEvent.getByLabel(inputTag,input);
   
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > inputNoise;
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > inputSignal;
   if ( doNoiseAndSignal )
   {
     iEvent.getByLabel(inputTagNoise,inputNoise);
     iEvent.getByLabel(inputTagSignal,inputSignal);
   }

   int medianTotClus =0;
   int per25TotClus =0;
   int signalTotClus =0;

   //----------------------------------
   //
   // MAIN LOOP OVER MODULES
   //
   //----------------------------------

   for ( edm::DetSetVector<SiStripRawDigi>::const_iterator 
         rawDigis = input->begin(); rawDigis != input->end(); rawDigis++) 
   {
     
     // ---
     // Set up some vectors to store digis at each stage of processing
     // and other variables for stuff

     std::vector<int16_t> pedestalSubDigis(rawDigis->size());

     std::vector<int16_t> medianSubtractedDigis(rawDigis->size());
     std::vector<int16_t> per25SubtractedDigis(rawDigis->size());
     std::vector<int16_t> signalSubtractedDigis( rawDigis->size()); // Digi Truth
     std::vector<float> medianOffset;
     std::vector<float> per25Offset;

     edm::DetSet<SiStripDigi> medianZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> per25ZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> signalZSDigis(rawDigis->id); // Digi Truth
     
     SiStripDetId sid( rawDigis->detId() );

     // The clusterizer wants to see entire DetSetVectors

     edm::DetSetVector<SiStripDigi> medianZSVec;
     edm::DetSetVector<SiStripDigi> per25ZSVec;
     edm::DetSetVector<SiStripDigi> signalZSVec;

     edmNew::DetSetVector<SiStripCluster> medianClusterVec;
     edmNew::DetSetVector<SiStripCluster> per25ClusterVec;
     edmNew::DetSetVector<SiStripCluster> signalClusterVec;

     // ---
     // Match real Digis with the signal only Digis
     // if they are available
     //

     edm::DetSetVector<SiStripRawDigi>::const_iterator  rawSignalDigis;
     if ( doNoiseAndSignal ) // Slow - should replace this with an Association Map!
     {
       rawSignalDigis = inputSignal->begin();
       for (rawSignalDigis = inputSignal->begin(); rawSignalDigis != inputSignal->end(); rawSignalDigis++)
         if( rawDigis->detId() == rawSignalDigis->detId() ) break;
     }   

    // ---
    // Perform pedestal subtraction
    //

    algorithms->subtractorPed->subtract( *rawDigis, pedestalSubDigis);

    // ---
    // Perform CMN subtraction and store algorithm results by APV
    //

    subtractMedian( pedestalSubDigis, medianSubtractedDigis, medianOffset );
    subtractPer25( pedestalSubDigis, per25SubtractedDigis, per25Offset );

    // for the signal digis, there is no baseline to remove. 
    // still need to convert to <int16_t> and fill with zeros if something is funky.
    // or if there are no signals 

    if ( doNoiseAndSignal ) {
      if ( rawSignalDigis != inputSignal->end() && rawSignalDigis->size() == rawDigis->size() ) {
        std::vector<int16_t>::iterator sD = signalSubtractedDigis.begin();
        edm::DetSet<SiStripRawDigi>::const_iterator rSD = rawSignalDigis->begin();
        while ( rSD != rawSignalDigis->end() ) {
          *sD = rSD->adc();
          ++sD;
          ++rSD;
        }
      } else { 
        for( std::vector<int16_t>::iterator sD = signalSubtractedDigis.begin(); sD != signalSubtractedDigis.end(); sD++ )
          *sD = 0;
      }
    }

    // ---
    // Now perform zero subtraction on the module for each CMN remover algorithm 
    //

    algorithms->suppressor->suppress( medianSubtractedDigis, medianZSDigis );
    algorithms->suppressor->suppress( per25SubtractedDigis, per25ZSDigis );
    algorithms->suppressor->suppress( signalSubtractedDigis, signalZSDigis );

    // ---
    // Run the clusterizer on the ZS digis
    //

    medianZSVec.insert(medianZSDigis);
    per25ZSVec.insert(per25ZSDigis);
    signalZSVec.insert(signalZSDigis);

    clusterAlgorithm->clusterize(medianZSVec, medianClusterVec);
    clusterAlgorithm->clusterize(per25ZSVec, per25ClusterVec);
    clusterAlgorithm->clusterize(signalZSVec, signalClusterVec);

    medianTotClus += medianClusterVec.dataSize();
    per25TotClus += per25ClusterVec.dataSize();
    signalTotClus += signalClusterVec.dataSize();

    //---------------------------------------------------
    //
    // ----- LOOP OVER EACH APV AND FILL HISTOGRAMS -----
    //
    //---------------------------------------------------

    for ( int APV = 0; APV < (int)rawDigis->size() / 128 ; APV++ ) 
    { 
      // ---
      // Get the pedestal subtracted ADCs for plotting
      //

      std::vector<int16_t> adc;  adc.reserve(128);  
      std::vector<int16_t> adcSignal; adcSignal.reserve(128);
      { 
        std::vector<int16_t>::const_iterator beginAPV, endAPV;
        beginAPV = pedestalSubDigis.begin(); beginAPV += 128*APV;
        endAPV = beginAPV + 128;
        adc.insert(adc.end(),beginAPV, endAPV);
      }
      if ( doNoiseAndSignal) 
      {
        std::vector<int16_t>::const_iterator beginAPV, endAPV;
        beginAPV = signalSubtractedDigis.begin(); beginAPV += 128*APV;
        endAPV = beginAPV + 128;
        adcSignal.insert(adcSignal.end(),beginAPV, endAPV);
      }

      // ---
      // Fill calclated baseline plots for all algorthims for all detectors and for TIB, TOB, TID, TEC
      //

      allMedians_->Fill( medianOffset[APV] );
      all25ths_->Fill( per25Offset[APV] );

      switch ( sid.subdetId() ) 
      {
        case SiStripDetId::TIB:
          TIBMedians_->Fill(medianOffset[APV]);
          TIB25ths_->Fill(per25Offset[APV]);
          TIBGens_->Fill( gaussDistribution->fire(128.,cmnTIBRMS_));
          break;
        case SiStripDetId::TOB:
          TOBMedians_->Fill(medianOffset[APV]);
          TOB25ths_->Fill(per25Offset[APV]);
          TOBGens_->Fill( gaussDistribution->fire(128.,cmnTOBRMS_));
          break;
        case SiStripDetId::TID:
          TIDMedians_->Fill(medianOffset[APV]);
          TID25ths_->Fill(per25Offset[APV]);
          TIDGens_->Fill( gaussDistribution->fire(128.,cmnTIDRMS_));
          break;
        case SiStripDetId::TEC:
          TECMedians_->Fill(medianOffset[APV]);
          TEC25ths_->Fill(per25Offset[APV]);
          TECGens_->Fill( gaussDistribution->fire(128.,cmnTECRMS_));
          break;
      }      
      medvs25_->Fill( medianOffset[APV], per25Offset[APV]);

      // ---
      // Make occupancy histos
      //
      medianOccupancy_->Fill( digiCount(medianZSDigis, APV), medianOffset[APV] );
      per25Occupancy_->Fill( digiCount(per25ZSDigis, APV), per25Offset[APV] );
      if(doNoiseAndSignal)
      {
        medianTrueOccupancy_->Fill( digiCount(signalZSDigis, APV), medianOffset[APV] );
        per25TrueOccupancy_->Fill( digiCount(signalZSDigis, APV), per25Offset[APV] );
      }

      // ---
      // Integrate ADC counts and fill histos
      //

      float intMed = integrateADC(adc, medianOffset[APV] );
      float intPer25 = integrateADC(adc, per25Offset[APV] );

      totalADCMedian_->Fill( intMed );
      totalADCPer25_->Fill( intPer25 );
      if ( doNoiseAndSignal )
      {
         float intSignal = integrateADC(adcSignal, 0. );
         totalADCSignal_->Fill( intSignal );
         totalADCMedianDiff_->Fill( intMed - intSignal );
         totalADCPer25Diff_->Fill( intPer25 - intSignal );
      }
  
      // ---
      // Make the gallery plot
      //
   
      // gallery A
      if ( (medianOffset[APV] >= 116 && medianOffset[APV] < 140) && 
           (per25Offset[APV] >= 116 && per25Offset[APV] < 140) )
      {
         int currentGraph = galAcount % 10;
         bool doGallery = true;
         if ( galAcount / 10 > 0 )
           if ( gaussDistribution->fire(0.,1.) < galAcount / 10 )
              doGallery = false;
       
         if ( doGallery ) {
           for( int i =0; i<128; i++) {
             galACount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adc[i]));
             if(doNoiseAndSignal) galASignalCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adcSignal[i]));
           }
           galACount_[currentGraph]->SetPoint(128,150.,(double)medianOffset[APV]);
           galACount_[currentGraph]->SetPoint(129,150.,(double)per25Offset[APV]);
           makeClusterGraph( medianClusterVec, galAMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galAPer25Clusters_[currentGraph], APV );
           galAcount++;
         }
      }

      // gallery B
      if ( (medianOffset[APV] >= 2 && medianOffset[APV] < 100) &&
           (per25Offset[APV] >= 2 && per25Offset[APV] < 100) )
      {
         int currentGraph = galBcount % 10;
         bool doGallery = true;
         if ( galBcount / 10 > 0 )
           if ( gaussDistribution->fire(0.,1.) < galBcount / 10 )
              doGallery = false;

         if ( doGallery ) {
           for( int i =0; i<128; i++) {
             galBCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adc[i]));
             if(doNoiseAndSignal) galBSignalCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adcSignal[i]));
           }
           galBCount_[currentGraph]->SetPoint(128,150.,(double)medianOffset[APV]);
           galBCount_[currentGraph]->SetPoint(129,150.,(double)per25Offset[APV]);
           makeClusterGraph( medianClusterVec, galBMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galBPer25Clusters_[currentGraph], APV );
           galBcount++;
         }
      }

      // gallery D
      if ( (medianOffset[APV] >= 160 ) &&
           (per25Offset[APV] >= 160 ) )
      {
         int currentGraph = galDcount % 10;
         bool doGallery = true;
         if ( galDcount / 10 > 0 )
           if ( gaussDistribution->fire(0.,1.) < galDcount / 10 )
              doGallery = false;

         if ( doGallery ) {
           for( int i =0; i<128; i++) {
             galDCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adc[i]));
             if(doNoiseAndSignal) galDSignalCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adcSignal[i]));
           }
           galDCount_[currentGraph]->SetPoint(128,150.,(double)medianOffset[APV]);
           galDCount_[currentGraph]->SetPoint(129,150.,(double)per25Offset[APV]);
           makeClusterGraph( medianClusterVec, galDMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galDPer25Clusters_[currentGraph], APV );
           galDcount++;
         }
      }

      // gallery E
      if ( (medianOffset[APV] >= 160) &&
           (per25Offset[APV] >= 116 && per25Offset[APV] < 140) )
      {
         int currentGraph = galEcount % 10;
         bool doGallery = true;
         if ( galEcount / 10 > 0 )
           if ( gaussDistribution->fire(0.,1.) < galEcount / 10 )
              doGallery = false;

         if ( doGallery ) {
           for( int i =0; i<128; i++) {
             galECount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adc[i]));
             if(doNoiseAndSignal) galESignalCount_[currentGraph]->SetPoint(i,(double)i+1,(double)(adcSignal[i]));
           }
           galECount_[currentGraph]->SetPoint(128,150.,(double)medianOffset[APV]);
           galECount_[currentGraph]->SetPoint(129,150.,(double)per25Offset[APV]);
           makeClusterGraph( medianClusterVec, galEMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galEPer25Clusters_[currentGraph], APV );
           galEcount++;
         }
      }


    }  // END HISTOGRAMMING LOOP OVER APVS IN THE MODULE
    


  } // END LOOP OVER MODULES

    
  // ---
  // Get the total number of clusters
  //
   
  medianTotalClus_->Fill( medianTotClus );
  per25TotalClus_->Fill( per25TotClus );
  if(doNoiseAndSignal) signalTotalClus_->Fill( signalTotClus );


}


void 
SiStripCMNAnalyzer::beginJob()
{
}

void 
SiStripCMNAnalyzer::endJob() {
}


void SiStripCMNAnalyzer::subtractMedian( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results)
{
  std::vector<int16_t> tmp;  tmp.reserve(128);  
    std::vector<int16_t>::iterator  
    strip( in.begin() ),
    stripOut ( out.begin() ), 
    end(   in.end()   ),
    endAPV;
  
  while( strip < end ) {
    endAPV = strip+128; tmp.clear();
    tmp.insert(tmp.end(),strip,endAPV);
    const float offset = median(tmp);
    results.push_back( offset );

    while (strip < endAPV) {
      *stripOut = *strip-offset;
      strip++; stripOut++;
    }
  }
}    

void SiStripCMNAnalyzer::subtractPer25( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results)
{
  std::vector<int16_t> tmp;  tmp.reserve(128);
    std::vector<int16_t>::iterator 
    strip( in.begin() ),  
    stripOut( out.begin() ), 
    end(   in.end()   ),
    endAPV;
  
  while( strip < end ) {
    endAPV = strip+128; tmp.clear();
    tmp.insert(tmp.end(),strip,endAPV);
    const float offset = percentile(tmp, 25.);
    results.push_back( offset );

    while (strip < endAPV) {
      *stripOut = *strip-offset;
      strip++; stripOut++;
    }
  }
}




inline float SiStripCMNAnalyzer::median( std::vector<int16_t>& sample ) {
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


int SiStripCMNAnalyzer::digiCount( edm::DetSet<SiStripDigi>& module, int APV )
{
   int total = 0;
  
   for( edm::DetSet<SiStripDigi>::const_iterator di = module.begin(); di != module.end(); ++di )
   {
     if ( (int)di->strip() >= APV*128 && (int)di->strip() < (APV+1)*128 )
       total++;
   }
   return total;
}

void SiStripCMNAnalyzer::makeClusterGraph( edmNew::DetSetVector<SiStripCluster>& clusters, TGraph * graph, int APV )
{
   bool firstDet = true;
   int position = 0;
 
   std::cout << "\n\n\n\n***Graphing cluster Vector with size " << clusters.size() << " and datasize " << clusters.dataSize() << "\n";


   // zero out the old points
   for(int i=0 ; i < 128; i++)
     graph->SetPoint( i, 0., 0. );

   edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters.begin();
   for ( ; it != clusters.end(); ++it )
   {
     // there should only be one DetSet
     if(firstDet == false ) break;
     firstDet = false;
     for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
     {
       if( position > 127) break;  // no more than 128 clusters on a strip!
       int firststrip = clus->firstStrip();
       int laststrip = firststrip + clus->amplitudes().size();
       std::cout << "   Cluster starting at " << firststrip << " ending at " << laststrip << "\n";
       // only fill clusters that are fully or partially on the current APV
       // and make strip numbers relative to the APV
       if( ( firststrip >= APV*128 && firststrip < (APV+1)*128 ) ||
           ( laststrip >= APV*128 && laststrip < (APV+1)*128 )   )
       {
         graph->SetPoint( position, firststrip - APV*128, laststrip - APV*128 );
         position++;
       }
     }

   }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripCMNAnalyzer);
