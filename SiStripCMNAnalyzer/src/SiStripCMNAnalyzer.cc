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


#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"
#include "CLHEP/Random/RandomEngine.h"

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"

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
      inline float pairMedian( std::vector<std::pair<float,float> >& sample);
      inline float percentile ( std::vector<int16_t>& sample, double pct);
      inline float integrateADC ( std::vector<int16_t>& sample, double offset);

      void subtractMedian( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results);
      void subtractPer25( std::vector<int16_t>& in, std::vector<int16_t>& out, std::vector<float>& results);
      void subtractIterMed( std::vector<int16_t>& in, std::vector<int16_t>& out, 
                            std::vector<float>& results, uint32_t detId);

      int digiCount( edm::DetSet<SiStripDigi>& module, int APV );

      void processRaw(const edm::InputTag&, const edm::DetSetVector<SiStripRawDigi>&, std::vector<edm::DetSet<SiStripDigi> >& );

      // As a dumb hack, clusters for the gallery plot are stored in a TGraph with  
      // (x,y) = (first strip, last strip) or (0,0) if there is no cluster.

      void makeClusterGraph( edmNew::DetSetVector<SiStripCluster>& clusters, TGraph * graph, int APV );

      void fillClusterWidths( edmNew::DetSetVector<SiStripCluster>& clusters, TH1F * hist, int APV );
      int countClusters( edmNew::DetSetVector<SiStripCluster>&clusters, int APV);

      std::auto_ptr<SiStripRawProcessingAlgorithms> algorithms;
      std::auto_ptr<StripClusterizerAlgorithm> clusterAlgorithm;

      edm::InputTag inputTag;
      edm::InputTag inputTagNoise;
      edm::InputTag inputTagSignal;
      bool doNoiseAndSignal;

      CLHEP::HepRandomEngine* rndEngine;
      CLHEP::RandGaussQ* gaussDistribution;
      GaussianTailNoiseGenerator* genNoise;

      edm::ESHandle<SiStripNoises> noiseHandle;
      edm::ESHandle<SiStripQuality> qualityHandle;
      uint32_t noise_cache_id, quality_cache_id;


      float cmnTIBRMS_;
      float cmnTOBRMS_;
      float cmnTIDRMS_;
      float cmnTECRMS_;


      TH1F* allMedians_;
      TH1F* all25ths_;
      TH1F* allIterMeds_;
      TH1F* TIBMedians_;
      TH1F* TIB25ths_;
      TH1F* TIBIterMeds_;
      TH1F* TIBGens_;
      TH1F* TOBMedians_;
      TH1F* TOB25ths_;
      TH1F* TOBIterMeds_;
      TH1F* TOBGens_;
      TH1F* TIDMedians_;
      TH1F* TID25ths_;
      TH1F* TIDIterMeds_;
      TH1F* TIDGens_;
      TH1F* TECMedians_;
      TH1F* TEC25ths_;
      TH1F* TECIterMeds_;
      TH1F* TECGens_;
      TH2F* medvs25_;

      // Occupancy Comparison

      TH2F* medianOccupancy_;
      TH2F* medianTrueOccupancy_;
      TH2F* per25Occupancy_;
      TH2F* per25TrueOccupancy_;
      TH2F* iterMedOccupancy_;
      TH2F* iterMedTrueOccupancy_;

      // Cluster width

      TH1F* medianClusWidth_;
      TH1F* iterMedClusWidth_;
      TH1F* per25ClusWidth_;

      // Clusters per APV

      TH1F* medianAPVClusNum_;
      TH1F* iterMedAPVClusNum_;
      TH1F* per25APVClusNum_;
 
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
      TGraph* galAIterMedClusters_[10]; 
      TGraph* galBMedianClusters_[10];
      TGraph* galBPer25Clusters_[10]; 
      TGraph* galBIterMedClusters_[10]; 
      TGraph* galDMedianClusters_[10];
      TGraph* galDPer25Clusters_[10]; 
      TGraph* galDIterMedClusters_[10]; 
      TGraph* galEMedianClusters_[10];
      TGraph* galEPer25Clusters_[10]; 
      TGraph* galEIterMedClusters_[10]; 
 
      // Histograms of integrated ADC  on each APV

      TH1F* totalADCMedian_;
      TH1F* totalADCPer25_;
      TH1F* totalADCIterMed_;

      // for comparison to actual integrated ADC from signal only with perfect baseline

      TH1F* totalADCSignal_;
      TH1F* totalADCMedianDiff_;
      TH1F* totalADCPer25Diff_;
      TH1F* totalADCIterMedDiff_;

      // Cluster Comparison
 
      TH1F* medianTotalClus_;
      TH1F* per25TotalClus_;
      TH1F* iterMedTotalClus_;
      TH1F* signalTotalClus_;


      int galAcount;
      int galBcount;
      int galDcount;
      int galEcount;

};

SiStripCMNAnalyzer::SiStripCMNAnalyzer(const edm::ParameterSet& iConfig)
  :  algorithms(SiStripRawProcessingFactory::create(iConfig.getParameter<edm::ParameterSet>("Algorithms"))),
     clusterAlgorithm( StripClusterizerAlgorithmFactory::create(iConfig.getParameter<edm::ParameterSet>("Clusterizer")) ), 
     inputTag(iConfig.getParameter<edm::InputTag>("RawDigiProducersList")),
     inputTagNoise(iConfig.getParameter<edm::InputTag>("RawDigiProducersListNoise")),
     inputTagSignal(iConfig.getParameter<edm::InputTag>("RawDigiProducersListSignal")),
     doNoiseAndSignal(iConfig.getParameter<bool>("doNoiseAndSignal"))
{

  edm::Service<TFileService> fs;
  allMedians_=fs->make<TH1F>("allMedians","Median ADC count of APV25", 500, 0., 255.);
  allIterMeds_=fs->make<TH1F>("allIterMeds","Iterated Median ADC count of APV25", 500, 0., 255.);
  all25ths_=fs->make<TH1F>("all25ths","25th percentile ADC count of APV25", 500, 0., 255.);
  medvs25_=fs->make<TH2F>("medvs25","Median ADC count versus 25th percentile",500,0.,500.,500,0.,255.);

  medianOccupancy_=fs->make<TH2F>("medianOccupancy","Median ADC count of APV25 versus Reconstructed # of Unsuppressed Strips",
                                  128.,0.,128.,500,0.,500.);
  iterMedOccupancy_=fs->make<TH2F>("iterMedOccupancy","Iterated Median ADC count of APV25 versus Reconstructed # of Unsuppressed Strips",
                                  128.,0.,128.,500,0.,500.);
  per25Occupancy_=fs->make<TH2F>("per25Occupancy","25th Percentile ADC count of APV25 versus Reconstructed # of Unsuppressed Strips",
                                  128.,0.,128.,500,0.,500.);
  if(doNoiseAndSignal)
  {
    medianTrueOccupancy_=fs->make<TH2F>("medianTrueOccupancy","Median ADC count of APV25 versus MC Correct # of Unsuppressed Strips",
                                    128.,0.,128.,500,0.,500.);
    iterMedTrueOccupancy_=fs->make<TH2F>("iterMedTrueOccupancy","Iterated Median ADC count of APV25 versus MC Correct # of Unsuppressed Strips",
                                    128.,0.,128.,500,0.,500.);
    per25TrueOccupancy_=fs->make<TH2F>("per25TrueOccupancy","25th Percentile ADC count of APV25 versus MC Correct # of Unsuppressed Strips",
                                    128.,0.,128.,500,0.,500.);
  }

  TIBMedians_=fs->make<TH1F>("TIBMedians","Median ADC count of APV25 for the TIB", 500, 0., 255.);
  TIBIterMeds_=fs->make<TH1F>("TIBIterMeds","Iterated Median ADC count of APV25 for the TIB", 500, 0., 255.);
  TIB25ths_=fs->make<TH1F>("TIB25ths","25th percentile ADC count of APV25 for the TIB", 500, 0., 255.);
  TIBGens_=fs->make<TH1F>("TIBGens","Expected CMN Distribution for the TIB", 500, 0., 255.);

  TOBMedians_=fs->make<TH1F>("TOBMedians","Median ADC count of APV25 for the TOB", 500, 0., 255.);
  TOBIterMeds_=fs->make<TH1F>("TOBIterMeds","Iterated Median ADC count of APV25 for the TOB", 500, 0., 255.);
  TOB25ths_=fs->make<TH1F>("TOB25ths","25th percentile ADC count of APV25 for the TOB", 500, 0., 255.);
  TOBGens_=fs->make<TH1F>("TOBGens","Expected CMN Distribution for the TOB", 500, 0., 255.);

  TIDMedians_=fs->make<TH1F>("TIDMedians","Median ADC count of APV25 for the TID", 500, 0., 255.);
  TIDIterMeds_=fs->make<TH1F>("TIDIterMeds","Iterated Median ADC count of APV25 for the TID", 500, 0., 255.);
  TID25ths_=fs->make<TH1F>("TID25ths","25th percentile ADC count of APV25 for the TID", 500, 0., 255.);
  TIDGens_=fs->make<TH1F>("TIDGens","Expected CMN Distribution for the TID", 500, 0., 255.);


  TECMedians_=fs->make<TH1F>("TECMedians","Median ADC count of APV25 for the TEC", 500, 0., 255.);
  TECIterMeds_=fs->make<TH1F>("TECIterMeds","Iterated Median ADC count of APV25 for the TEC", 500, 0., 255.);
  TEC25ths_=fs->make<TH1F>("TEC25ths","25th percentile ADC count of APV25 for the TEC", 500, 0., 255.);
  TECGens_=fs->make<TH1F>("TECGens","Expected CMN Distribution for the TEC", 500, 0., 255.);

  totalADCMedian_=fs->make<TH1F>("totalADCMedian","Integrated ADC count for each APV25 after Median CMN subtraction",500,0.,10000.);
  totalADCIterMed_=fs->make<TH1F>("totalADCIterMed","Integrated ADC count for each APV25 after Iterated Median CMN subtraction",500,0.,10000.);
  totalADCPer25_=fs->make<TH1F>("totalADCPer25","Integrated ADC count for each APV25 after 25th Percentile CMN subtraction",500,0.,10000.);
  if(doNoiseAndSignal)
  {
    totalADCSignal_=fs->make<TH1F>("totalADCSignal","Integrated ADC count for each APV25 from true SimHits only",500,0.,10000.);
    totalADCMedianDiff_=fs->make<TH1F>("totalADCMedianDiff",
      "Difference between median subtracted integrated ADC count and ADC count from signal",500,-3000.,2000.);
    totalADCIterMedDiff_=fs->make<TH1F>("totalADCIterMedDiff",
      "Difference between iterated median subtracted integrated ADC count and ADC count from signal",500,-3000.,2000.);
    totalADCPer25Diff_=fs->make<TH1F>("totalADCPer25Diff",
      "Difference between 25th percentile subtracted integrated ADC count and ADC count from signal",500,-3000.,2000.);
  }    
  
  medianTotalClus_=fs->make<TH1F>("medianTotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);
  iterMedTotalClus_=fs->make<TH1F>("iterMedTotalClus","Number of clusters found using iterated median CMN subtraction",500,0,500000.);
  per25TotalClus_=fs->make<TH1F>("per25TotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);
  if(doNoiseAndSignal)
    signalTotalClus_=fs->make<TH1F>("signalTotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);

  medianClusWidth_ = fs->make<TH1F>("medianClusWidth","Cluster Width - Median CMN Subtraction",256,0.,256.);
  iterMedClusWidth_ = fs->make<TH1F>("iterMedClusWidth","Cluster Width - Iterated Median CMN Subtraction",256,0.,256.);
  per25ClusWidth_ = fs->make<TH1F>("per25ClusWidth","Cluster Width - 25th Percentile CMN Subtraction",256,0.,256.);

  medianAPVClusNum_ = fs->make<TH1F>("medianAPVClusNum","Number of Clusters on an APV25 - Median CMN Subtraction",128,0.,128.);
  iterMedAPVClusNum_ = fs->make<TH1F>("iterMedAPVClusNum","Number of Clusters on an APV25 - Iterated Median CMN Subtraction",128,0.,128.);
  per25APVClusNum_ = fs->make<TH1F>("per25APVClusNum","Number of Clusters on an APV25 - 25th Percentile CMN Subtraction",128,0.,128.);
  

  for( int i = 0; i<10; i++)
  {
    galACount_[i] = fs->make<TGraph>(131);
    galACount_[i]->SetName(Form("galACount%d",i));
    galAMedianClusters_[i] = fs->make<TGraph>(128);
    galAMedianClusters_[i]->SetName(Form("galAMedianClusters%d",i));
    galAIterMedClusters_[i] = fs->make<TGraph>(128);
    galAIterMedClusters_[i]->SetName(Form("galAIterMedClusters%d",i));
    galAPer25Clusters_[i] = fs->make<TGraph>(128);
    galAPer25Clusters_[i]->SetName(Form("galAPer25Clusters%d",i));
    galBCount_[i] = fs->make<TGraph>(131);
    galBCount_[i]->SetName(Form("galBCount%d",i));
    galBMedianClusters_[i] = fs->make<TGraph>(128);
    galBMedianClusters_[i]->SetName(Form("galBMedianClusters%d",i));
    galBIterMedClusters_[i] = fs->make<TGraph>(128);
    galBIterMedClusters_[i]->SetName(Form("galBIterMedClusters%d",i));
    galBPer25Clusters_[i] = fs->make<TGraph>(128);
    galBPer25Clusters_[i]->SetName(Form("galBPer25Clusters%d",i));
    galDCount_[i] = fs->make<TGraph>(131);
    galDCount_[i]->SetName(Form("galDCount%d",i));
    galDMedianClusters_[i] = fs->make<TGraph>(128);
    galDMedianClusters_[i]->SetName(Form("galDMedianClusters%d",i));
    galDIterMedClusters_[i] = fs->make<TGraph>(128);
    galDIterMedClusters_[i]->SetName(Form("galDIterMedClusters%d",i));
    galDPer25Clusters_[i] = fs->make<TGraph>(128);
    galDPer25Clusters_[i]->SetName(Form("galDPer25Clusters%d",i));
    galECount_[i] = fs->make<TGraph>(131);
    galECount_[i]->SetName(Form("galECount%d",i));
    galEMedianClusters_[i] = fs->make<TGraph>(128);
    galEMedianClusters_[i]->SetName(Form("galEMedianClusters%d",i));
    galEIterMedClusters_[i] = fs->make<TGraph>(128);
    galEIterMedClusters_[i]->SetName(Form("galEIterMedClusters%d",i));
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

   uint32_t n_cache_id = iSetup.get<SiStripNoisesRcd>().cacheIdentifier();
   uint32_t q_cache_id = iSetup.get<SiStripQualityRcd>().cacheIdentifier();

   if(n_cache_id != noise_cache_id) {
     iSetup.get<SiStripNoisesRcd>().get( noiseHandle );
     noise_cache_id = n_cache_id;
   }
   if(q_cache_id != quality_cache_id) {
     iSetup.get<SiStripQualityRcd>().get( qualityHandle );
     quality_cache_id = q_cache_id;
   }


   int medianTotClus =0;
   int iterMedTotClus =0;
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
     std::vector<int16_t> iterMedSubtractedDigis( rawDigis->size());
     std::vector<int16_t> signalSubtractedDigis( rawDigis->size()); // Digi Truth
     std::vector<float> medianOffset;
     std::vector<float> iterMedOffset;
     std::vector<float> per25Offset;

     edm::DetSet<SiStripDigi> medianZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> per25ZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> iterMedZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> signalZSDigis(rawDigis->id); // Digi Truth
     
     SiStripDetId sid( rawDigis->detId() );

     // The clusterizer wants to see entire DetSetVectors

     edm::DetSetVector<SiStripDigi> medianZSVec;
     edm::DetSetVector<SiStripDigi> per25ZSVec;
     edm::DetSetVector<SiStripDigi> iterMedZSVec;
     edm::DetSetVector<SiStripDigi> signalZSVec;

     edmNew::DetSetVector<SiStripCluster> medianClusterVec;
     edmNew::DetSetVector<SiStripCluster> per25ClusterVec;
     edmNew::DetSetVector<SiStripCluster> iterMedClusterVec;
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
    subtractIterMed( pedestalSubDigis, iterMedSubtractedDigis, iterMedOffset, rawDigis->detId() );
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
    algorithms->suppressor->suppress( iterMedSubtractedDigis, iterMedZSDigis );
    algorithms->suppressor->suppress( per25SubtractedDigis, per25ZSDigis );
    algorithms->suppressor->suppress( signalSubtractedDigis, signalZSDigis );

    // ---
    // Run the clusterizer on the ZS digis
    //

    medianZSVec.insert(medianZSDigis);
    per25ZSVec.insert(per25ZSDigis);
    iterMedZSVec.insert(iterMedZSDigis);
    signalZSVec.insert(signalZSDigis);

    clusterAlgorithm->clusterize(medianZSVec, medianClusterVec);
    clusterAlgorithm->clusterize(iterMedZSVec, iterMedClusterVec);
    clusterAlgorithm->clusterize(per25ZSVec, per25ClusterVec);
    clusterAlgorithm->clusterize(signalZSVec, signalClusterVec);

    medianTotClus += medianClusterVec.dataSize();
    iterMedTotClus += iterMedClusterVec.dataSize();
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
      allIterMeds_->Fill( iterMedOffset[APV] );
      all25ths_->Fill( per25Offset[APV] );

      switch ( sid.subdetId() ) 
      {
        case SiStripDetId::TIB:
          TIBMedians_->Fill(medianOffset[APV]);
          TIBIterMeds_->Fill(iterMedOffset[APV]);
          TIB25ths_->Fill(per25Offset[APV]);
          TIBGens_->Fill( gaussDistribution->fire(128.,cmnTIBRMS_));
          break;
        case SiStripDetId::TOB:
          TOBMedians_->Fill(medianOffset[APV]);
          TOBIterMeds_->Fill(iterMedOffset[APV]);
          TOB25ths_->Fill(per25Offset[APV]);
          TOBGens_->Fill( gaussDistribution->fire(128.,cmnTOBRMS_));
          break;
        case SiStripDetId::TID:
          TIDMedians_->Fill(medianOffset[APV]);
          TIDIterMeds_->Fill(iterMedOffset[APV]);
          TID25ths_->Fill(per25Offset[APV]);
          TIDGens_->Fill( gaussDistribution->fire(128.,cmnTIDRMS_));
          break;
        case SiStripDetId::TEC:
          TECMedians_->Fill(medianOffset[APV]);
          TECIterMeds_->Fill(iterMedOffset[APV]);
          TEC25ths_->Fill(per25Offset[APV]);
          TECGens_->Fill( gaussDistribution->fire(128.,cmnTECRMS_));
          break;
      }      
      medvs25_->Fill( medianOffset[APV], per25Offset[APV]);

      // ---
      // Make occupancy histos
      //
      medianOccupancy_->Fill( digiCount(medianZSDigis, APV), medianOffset[APV] );
      iterMedOccupancy_->Fill( digiCount(iterMedZSDigis, APV), iterMedOffset[APV] );
      per25Occupancy_->Fill( digiCount(per25ZSDigis, APV), per25Offset[APV] );
      if(doNoiseAndSignal)
      {
        medianTrueOccupancy_->Fill( digiCount(signalZSDigis, APV), medianOffset[APV] );
        iterMedTrueOccupancy_->Fill( digiCount(signalZSDigis, APV), iterMedOffset[APV] );
        per25TrueOccupancy_->Fill( digiCount(signalZSDigis, APV), per25Offset[APV] );
      }

      // ---
      // Make cluster width and count histos
      //

      medianAPVClusNum_->Fill( countClusters( medianClusterVec, APV));
      iterMedAPVClusNum_->Fill( countClusters( iterMedClusterVec, APV));
      per25APVClusNum_->Fill( countClusters( per25ClusterVec, APV));
 
      fillClusterWidths( medianClusterVec, medianClusWidth_, APV); 
      fillClusterWidths( iterMedClusterVec, iterMedClusWidth_, APV); 
      fillClusterWidths( per25ClusterVec, per25ClusWidth_, APV); 

      // ---
      // Integrate ADC counts and fill histos
      //

      float intMed = integrateADC(adc, medianOffset[APV] );
      float intIterMed = integrateADC(adc, iterMedOffset[APV] );
      float intPer25 = integrateADC(adc, per25Offset[APV] );

      totalADCMedian_->Fill( intMed );
      totalADCIterMed_->Fill( intIterMed);
      totalADCPer25_->Fill( intPer25 );
      if ( doNoiseAndSignal )
      {
         float intSignal = integrateADC(adcSignal, 0. );
         totalADCSignal_->Fill( intSignal );
         totalADCMedianDiff_->Fill( intMed - intSignal );
         totalADCIterMedDiff_->Fill( intIterMed - intSignal );
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
           galACount_[currentGraph]->SetPoint(130,150.,(double)iterMedOffset[APV]);
           makeClusterGraph( medianClusterVec, galAMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galAPer25Clusters_[currentGraph], APV );
           makeClusterGraph( iterMedClusterVec, galAIterMedClusters_[currentGraph], APV );
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
           galBCount_[currentGraph]->SetPoint(130,150.,(double)iterMedOffset[APV]);
           makeClusterGraph( medianClusterVec, galBMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galBPer25Clusters_[currentGraph], APV );
           makeClusterGraph( iterMedClusterVec, galBIterMedClusters_[currentGraph], APV );
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
           galDCount_[currentGraph]->SetPoint(130,150.,(double)iterMedOffset[APV]);
           makeClusterGraph( medianClusterVec, galDMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galDPer25Clusters_[currentGraph], APV );
           makeClusterGraph( iterMedClusterVec, galDIterMedClusters_[currentGraph], APV );
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
           galECount_[currentGraph]->SetPoint(130,150.,(double)iterMedOffset[APV]);
           makeClusterGraph( medianClusterVec, galEMedianClusters_[currentGraph], APV );
           makeClusterGraph( per25ClusterVec, galEPer25Clusters_[currentGraph], APV );
           makeClusterGraph( iterMedClusterVec, galEIterMedClusters_[currentGraph], APV );
           galEcount++;
         }
      }


    }  // END HISTOGRAMMING LOOP OVER APVS IN THE MODULE
    


  } // END LOOP OVER MODULES

    
  // ---
  // Get the total number of clusters
  //
   
  medianTotalClus_->Fill( medianTotClus );
  iterMedTotalClus_->Fill( iterMedTotClus );
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


void SiStripCMNAnalyzer::subtractIterMed( std::vector<int16_t>& in, std::vector<int16_t>& out, 
                                          std::vector<float>& results, uint32_t detId)
{
   // for now we will do 2-sigma and 3 iterations. This should be given in the configuration.
   
  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId);
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId);

  std::vector<int16_t>::iterator fs,ls,os;
  os = out.begin();
  float offset = 0;  
  std::vector< std::pair<float,float> > subset;
  subset.reserve(128); 


  for( int16_t APV=0; APV< (int)in.size()/128; ++APV)
  {
    subset.clear();
    // fill subset vector with all good strips and their noises
    for (int16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip)
    {
      if ( !qualityHandle->IsStripBad(detQualityRange,istrip) )
      {
        std::pair<float,float> pin((float)in[istrip], (float)noiseHandle->getNoise(istrip,detNoiseRange));
        subset.push_back( pin );
      }
    }
    // caluate offset for all good strips (first iteration)
    if (subset.size() != 0)
      offset = pairMedian(subset);

    // for second, third... iterations, remove strips over threshold
    // and recalculate offset on remaining strips
    for ( int ii = 0; ii<2; ++ii )
    {
      std::vector< std::pair<float,float> >::iterator si = subset.begin();
      while(  si != subset.end() )
      {
        if( si->first-offset > 2*si->second )  
          si = subset.erase(si);
        else
          ++si;
      }

      if ( subset.size() == 0 ) break;
      offset = pairMedian(subset);
    }        

    // remove offset
    results.push_back(offset);
    fs = in.begin()+APV*128;
    ls = in.begin()+(APV+1)*128;
    while (fs < ls)
    {
      *os = *fs-offset;
      ++os;
      ++fs;
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

inline float SiStripCMNAnalyzer::pairMedian( std::vector<std::pair<float,float> >& sample) {
  std::vector<std::pair<float,float> >::iterator mid = sample.begin() + sample.size()/2;
  std::nth_element(sample.begin(), mid, sample.end());
  if( sample.size() & 1 ) //odd size
    return (*mid).first;
  return ( (*std::max_element(sample.begin(), mid)).first + (*mid).first ) / 2.;
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

void SiStripCMNAnalyzer::fillClusterWidths( edmNew::DetSetVector<SiStripCluster>& clusters, TH1F * hist, int APV )
{
   bool firstDet = true;
   edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters.begin();
   for ( ; it != clusters.end(); ++it )
   {
     // there should only be one DetSet
     if(firstDet == false ) break;
     firstDet = false;
     for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
     {
       int firststrip = clus->firstStrip();
       int laststrip = firststrip + clus->amplitudes().size();
       // only fill clusters that begin  on the current APV
       // do not want to double count clusters stretching over multiple APVs
       if(  firststrip >= APV*128 && firststrip < (APV+1)*128 )
       {
         hist->Fill(laststrip-firststrip);
       }
     }
   }
}

int SiStripCMNAnalyzer::countClusters( edmNew::DetSetVector<SiStripCluster>& clusters, int APV)
{
   int count = 0;
   bool firstDet = true;
   edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters.begin();
   for ( ; it != clusters.end(); ++it )
   {
     // there should only be one DetSet
     if(firstDet == false ) break;
     firstDet = false;
     for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
     {
       int firststrip = clus->firstStrip();
       int laststrip = firststrip + clus->amplitudes().size();
       // only fill clusters that begin  on the current APV
       // do not want to double count clusters stretching over multiple APVs
       if( ( firststrip >= APV*128 && firststrip < (APV+1)*128 ) ||
           ( laststrip >= APV*128 && laststrip < (APV+1)*128 )   )
         count++;
     }
   }
   return count;
}


void SiStripCMNAnalyzer::makeClusterGraph( edmNew::DetSetVector<SiStripCluster>& clusters, TGraph * graph, int APV )
{
   bool firstDet = true;
   int position = 0;
 
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
