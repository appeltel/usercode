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
#include <TTree.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// APV readout structure for Gallery plots
//

typedef struct 
{
  int BLFOffset[128];
  int adc[128];
  int clustersMedian[128];
  int clustersIterMed[128];
  int clustersPer25[128];
  int clustersFastLin[128];
  int clustersSplitLin[128];
  int clustersBLF[128];
  double medianOffset;
  double iterMedOffset;
  double per25Offset;
  double fastLinOffset;
  double fastLinSlope;
  double splitLinOffsetA;
  double splitLinOffsetB;
  double splitLinSlopeA;
  double splitLinSlopeB;
  int event;
  int lumi;
  int run;
  uint32_t detID;
  int apv; 
} apvReadout_t;


//
// For BaselineFollower (BLF)
//

typedef std::map< uint16_t, int16_t> DigiMap;
typedef std::map< uint16_t, std::vector < int16_t> > RawDigiMap;
typedef std::map< uint16_t, int16_t>::iterator DigiMapIter;
typedef std::map<uint32_t, std::vector<float> > CMMap;  //detId, Vector of MeanCM per detId


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
      void subtractFastLin( std::vector<int16_t>& in, std::vector<int16_t>& out, 
                            std::vector<float>& offsets, std::vector<float>& slopes );
      void subtractSplitLin( std::vector<int16_t>& in, std::vector<int16_t>& out, 
                             std::vector<float>& offsetsA, std::vector<float>& offsetsB,
                             std::vector<float>& slopesA, std::vector<float>& slopesB );
      void subtractBLF( std::vector<int16_t>& in, std::vector<int16_t>& out,
                        std::vector<std::vector<float> >& results, uint32_t detId);

      int digiCount( edm::DetSet<SiStripDigi>& module, int APV );

      void processRaw(const edm::InputTag&, const edm::DetSetVector<SiStripRawDigi>&, std::vector<edm::DetSet<SiStripDigi> >& );

      void copyAPVReadout( apvReadout_t& src, apvReadout_t& dest );

      void fillClusterValues( edmNew::DetSetVector<SiStripCluster>& clusters, int * array, int APV );

      void fillClusterWidths( edmNew::DetSetVector<SiStripCluster>& clusters, TH1F * hist, int APV );
      int countClusters( edmNew::DetSetVector<SiStripCluster>&clusters, int APV);

      bool galleryCuts( double median, double per25, int gal );

      std::auto_ptr<SiStripRawProcessingAlgorithms> algorithms;
      std::auto_ptr<StripClusterizerAlgorithm> clusterAlgorithm;

      edm::InputTag inputTag;

      CLHEP::HepRandomEngine* rndEngine;
      CLHEP::RandGaussQ* gaussDistribution;
      GaussianTailNoiseGenerator* genNoise;

      edm::ESHandle<SiStripNoises> noiseHandle;
      edm::ESHandle<SiStripQuality> qualityHandle;
      uint32_t noise_cache_id, quality_cache_id;

      // baseline follower configuration variables
      uint32_t nSigmaNoiseDerTh_;          // threshold for rejecting hits strips
      uint32_t consecThreshold_;           // minimum length of flat region
      uint32_t hitStripThreshold_;         // height above median when strip is definitely a hit
      uint32_t nSmooth_;                   // for smoothing and local minimum determination (odd number)
      uint32_t minStripsToFit_;            // minimum strips to try spline algo (otherwise default to median)


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

      apvReadout_t gal_[5][10];
  
      TTree* galleryA_;
      TTree* galleryB_;
      TTree* galleryC_;
      TTree* galleryD_;
      TTree* galleryE_;

      apvReadout_t tmpAPV;

      // In case of something like ZeroBias, one may want to set a minimum
      // number of clusters before filling galleries
      int galleryClusterMin_;

      // Histograms of integrated ADC  on each APV

      TH1F* totalADCMedian_;
      TH1F* totalADCPer25_;
      TH1F* totalADCIterMed_;

      // Cluster number comparison
 
      TH1F* medianTotalClus_;
      TH1F* per25TotalClus_;
      TH1F* iterMedTotalClus_;


      int galCount[5];

};


#endif

