#ifndef Appeltel_SiStripCMNAnalyzer_SiStripCMNAnalyzer_h
#define Appeltel_SiStripCMNAnalyzer_SiStripCMNAnalyzer_h

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
      TH2F* per25Occupancy_;
      TH2F* iterMedOccupancy_;

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
      TGraph* galBCount_[10];
      TGraph* galDCount_[10];
      TGraph* galECount_[10];
  
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

      // Cluster number comparison
 
      TH1F* medianTotalClus_;
      TH1F* per25TotalClus_;
      TH1F* iterMedTotalClus_;


      int galAcount;
      int galBcount;
      int galDcount;
      int galEcount;

};


#endif
