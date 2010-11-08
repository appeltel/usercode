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
// $Id: SiStripCMNAnalyzer.cc,v 1.21 2010/09/09 09:59:11 appeltel Exp $
//
//

#include "Appeltel/SiStripCMNAnalyzer/interface/SiStripCMNAnalyzer.h"

#include <memory>
#include <iostream>

SiStripCMNAnalyzer::SiStripCMNAnalyzer(const edm::ParameterSet& iConfig)
  :  algorithms(SiStripRawProcessingFactory::create(iConfig.getParameter<edm::ParameterSet>("Algorithms"))),
     clusterAlgorithm( StripClusterizerAlgorithmFactory::create(iConfig.getParameter<edm::ParameterSet>("Clusterizer")) ), 
     inputTag(iConfig.getParameter<edm::InputTag>("RawDigiProducersList")),
     galleryClusterMin_(iConfig.exists("galleryClusterMin")?iConfig.getParameter<int>("galleryClusterMin"):0),
     nSigmaNoiseDerTh_(iConfig.getParameter<uint32_t>("nSigmaNoiseDerTh")),
     consecThreshold_(iConfig.getParameter<uint32_t>("consecThreshold")),
     hitStripThreshold_(iConfig.getParameter<uint32_t>("hitStripThreshold")),
     nSmooth_(iConfig.getParameter<uint32_t>("nSmooth")),
     minStripsToFit_(iConfig.getParameter<uint32_t>("minStripsToFit"))

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
  
  medianTotalClus_=fs->make<TH1F>("medianTotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);
  iterMedTotalClus_=fs->make<TH1F>("iterMedTotalClus","Number of clusters found using iterated median CMN subtraction",500,0,500000.);
  per25TotalClus_=fs->make<TH1F>("per25TotalClus","Number of clusters found using median CMN subtraction",500,0,500000.);

  medianClusWidth_ = fs->make<TH1F>("medianClusWidth","Cluster Width - Median CMN Subtraction",256,0.,256.);
  iterMedClusWidth_ = fs->make<TH1F>("iterMedClusWidth","Cluster Width - Iterated Median CMN Subtraction",256,0.,256.);
  per25ClusWidth_ = fs->make<TH1F>("per25ClusWidth","Cluster Width - 25th Percentile CMN Subtraction",256,0.,256.);

  medianAPVClusNum_ = fs->make<TH1F>("medianAPVClusNum","Number of Clusters on an APV25 - Median CMN Subtraction",128,0.,128.);
  iterMedAPVClusNum_ = fs->make<TH1F>("iterMedAPVClusNum","Number of Clusters on an APV25 - Iterated Median CMN Subtraction",128,0.,128.);
  per25APVClusNum_ = fs->make<TH1F>("per25APVClusNum","Number of Clusters on an APV25 - 25th Percentile CMN Subtraction",128,0.,128.);
  

  TString leafstring;
  leafstring = "BLFOffset[128]/I:adc[128]/I:clustersMedian[128]/I:clustersIterMed[128]/I:clustersPer25[128]/I:clustersFastLin[128]/I:clustersSplitLin[128]/I:clustersBLF[128]/I:medianOffset/D:iterMedOffset/D:per25Offset/D:fastLinOffset/D:fastLinSlope/D:splitLinOffsetA/D:splitLinOffsetB/D:splitLinSlopeA/D:splitLinSlopeB/D:event/I:lumi/I:run/I:detID/i:apv/I";

  galleryA_ = fs->make<TTree>("galleryATree","galleryATree");
  galleryA_->Branch("apvReadouts",&tmpAPV, leafstring.Data()); 
  galleryB_ = fs->make<TTree>("galleryBTree","galleryBTree");
  galleryB_->Branch("apvReadouts",&tmpAPV, leafstring.Data()); 
  galleryC_ = fs->make<TTree>("galleryCTree","galleryCTree");
  galleryC_->Branch("apvReadouts",&tmpAPV, leafstring.Data()); 
  galleryD_ = fs->make<TTree>("galleryDTree","galleryDTree");
  galleryD_->Branch("apvReadouts",&tmpAPV, leafstring.Data()); 
  galleryE_ = fs->make<TTree>("galleryETree","galleryETree");
  galleryE_->Branch("apvReadouts",&tmpAPV, leafstring.Data()); 

  for( int i = 0; i < 5; i++) galCount[i] = 0;

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

   int eventNumber = iEvent.id().event();
   int lumiNumber = iEvent.id().luminosityBlock();
   int runNumber = iEvent.id().run();

   int medianTotClus =0;
   int iterMedTotClus =0;
   int per25TotClus =0;

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
     std::vector<int16_t> iterMedSubtractedDigis(rawDigis->size());
     std::vector<int16_t> fastLinSubtractedDigis(rawDigis->size());
     std::vector<int16_t> splitLinSubtractedDigis(rawDigis->size()); 
     std::vector<int16_t> BLFSubtractedDigis(rawDigis->size());
     std::vector<float> medianOffset;
     std::vector<float> iterMedOffset;
     std::vector<float> per25Offset;
     std::vector<float> fastLinOffset;
     std::vector<float> fastLinSlope;
     std::vector<float> splitLinOffsetA;
     std::vector<float> splitLinOffsetB;
     std::vector<float> splitLinSlopeA;
     std::vector<float> splitLinSlopeB;
     std::vector<std::vector<float> > BLFOffset;

     edm::DetSet<SiStripDigi> medianZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> per25ZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> iterMedZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> fastLinZSDigis(rawDigis->id);     
     edm::DetSet<SiStripDigi> splitLinZSDigis(rawDigis->id);
     edm::DetSet<SiStripDigi> BLFZSDigis(rawDigis->id);

     SiStripDetId sid( rawDigis->detId() );

     // The clusterizer wants to see entire DetSetVectors

     edm::DetSetVector<SiStripDigi> medianZSVec;
     edm::DetSetVector<SiStripDigi> per25ZSVec;
     edm::DetSetVector<SiStripDigi> iterMedZSVec;
     edm::DetSetVector<SiStripDigi> fastLinZSVec;
     edm::DetSetVector<SiStripDigi> splitLinZSVec;
     edm::DetSetVector<SiStripDigi> BLFZSVec;

     edmNew::DetSetVector<SiStripCluster> medianClusterVec;
     edmNew::DetSetVector<SiStripCluster> per25ClusterVec;
     edmNew::DetSetVector<SiStripCluster> iterMedClusterVec;
     edmNew::DetSetVector<SiStripCluster> fastLinClusterVec;
     edmNew::DetSetVector<SiStripCluster> splitLinClusterVec;
     edmNew::DetSetVector<SiStripCluster> BLFClusterVec;

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
    subtractFastLin( pedestalSubDigis, fastLinSubtractedDigis, fastLinOffset, fastLinSlope );
    subtractSplitLin( pedestalSubDigis, splitLinSubtractedDigis, splitLinOffsetA, 
                      splitLinOffsetB, splitLinSlopeA, splitLinSlopeB );
    subtractBLF( pedestalSubDigis, BLFSubtractedDigis, BLFOffset, rawDigis->detId() );

    // ---
    // Now perform zero subtraction on the module for each CMN remover algorithm 
    //

    algorithms->suppressor->suppress( medianSubtractedDigis, medianZSDigis );
    algorithms->suppressor->suppress( iterMedSubtractedDigis, iterMedZSDigis );
    algorithms->suppressor->suppress( per25SubtractedDigis, per25ZSDigis );
    algorithms->suppressor->suppress( fastLinSubtractedDigis, fastLinZSDigis );
    algorithms->suppressor->suppress( splitLinSubtractedDigis, splitLinZSDigis );
    algorithms->suppressor->suppress( BLFSubtractedDigis, BLFZSDigis );

    // ---
    // Run the clusterizer on the ZS digis
    //

    medianZSVec.insert(medianZSDigis);
    per25ZSVec.insert(per25ZSDigis);
    iterMedZSVec.insert(iterMedZSDigis);
    fastLinZSVec.insert(fastLinZSDigis);
    splitLinZSVec.insert(splitLinZSDigis);
    BLFZSVec.insert(BLFZSDigis);

    clusterAlgorithm->clusterize(medianZSVec, medianClusterVec);
    clusterAlgorithm->clusterize(iterMedZSVec, iterMedClusterVec);
    clusterAlgorithm->clusterize(per25ZSVec, per25ClusterVec);
    clusterAlgorithm->clusterize(fastLinZSVec, fastLinClusterVec);
    clusterAlgorithm->clusterize(splitLinZSVec, splitLinClusterVec);
    clusterAlgorithm->clusterize(BLFZSVec, BLFClusterVec);

    medianTotClus += medianClusterVec.dataSize();
    iterMedTotClus += iterMedClusterVec.dataSize();
    per25TotClus += per25ClusterVec.dataSize();

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
      { 
        std::vector<int16_t>::const_iterator beginAPV, endAPV;
        beginAPV = pedestalSubDigis.begin(); beginAPV += 128*APV;
        endAPV = beginAPV + 128;
        adc.insert(adc.end(),beginAPV, endAPV);
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
  
      // ---
      // Make the gallery plots
      //
      if (countClusters( medianClusterVec, APV) >= galleryClusterMin_ ) 
      {
        for ( int galNum = 0; galNum < 5; galNum++ )
        {
          if ( galleryCuts( medianOffset[APV], per25Offset[APV], galNum ) )
          {
            int currentGraph = galCount[galNum] % 10;
            bool doGallery = true;
            if ( galCount[galNum] / 10 > 0 )
              if ( gaussDistribution->fire(0.,1.) < galCount[galNum] / 10 )
                 doGallery = false;

            if ( doGallery ) {
              for ( int i = 0; i<128; i++) 
                gal_[galNum][currentGraph].BLFOffset[i] = (int)BLFOffset[APV][i];
              for ( int i = 0; i<128; i++) gal_[galNum][currentGraph].adc[i] = adc[i];
              gal_[galNum][currentGraph].medianOffset = (double)medianOffset[APV];
              gal_[galNum][currentGraph].iterMedOffset = (double)iterMedOffset[APV];
              gal_[galNum][currentGraph].per25Offset = (double)per25Offset[APV];
              gal_[galNum][currentGraph].fastLinOffset = (double)fastLinOffset[APV];
              gal_[galNum][currentGraph].fastLinSlope = (double)fastLinSlope[APV];
              gal_[galNum][currentGraph].splitLinOffsetA = (double)splitLinOffsetA[APV];
              gal_[galNum][currentGraph].splitLinOffsetB = (double)splitLinOffsetB[APV];
              gal_[galNum][currentGraph].splitLinSlopeA = (double)splitLinSlopeA[APV];
              gal_[galNum][currentGraph].splitLinSlopeB = (double)splitLinSlopeB[APV];
              fillClusterValues( medianClusterVec, gal_[galNum][currentGraph].clustersMedian, APV );
              fillClusterValues( per25ClusterVec, gal_[galNum][currentGraph].clustersPer25, APV );
              fillClusterValues( iterMedClusterVec, gal_[galNum][currentGraph].clustersIterMed, APV );
              fillClusterValues( fastLinClusterVec, gal_[galNum][currentGraph].clustersFastLin, APV );
              fillClusterValues( splitLinClusterVec, gal_[galNum][currentGraph].clustersSplitLin, APV );
              fillClusterValues( BLFClusterVec, gal_[galNum][currentGraph].clustersBLF, APV );
              gal_[galNum][currentGraph].event = eventNumber;
              gal_[galNum][currentGraph].lumi = lumiNumber;
              gal_[galNum][currentGraph].run = runNumber;
              gal_[galNum][currentGraph].detID = rawDigis->detId();
              gal_[galNum][currentGraph].apv = APV;
              galCount[galNum]++;

            }
          }
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


}


void 
SiStripCMNAnalyzer::beginJob()
{
}

void 
SiStripCMNAnalyzer::endJob() {

  for( int i = 0; i< 10; i++)
  {
    copyAPVReadout( gal_[0][i], tmpAPV );
    galleryA_->Fill();
    copyAPVReadout( gal_[1][i], tmpAPV );
    galleryB_->Fill();
    copyAPVReadout( gal_[2][i], tmpAPV );
    galleryC_->Fill();
    copyAPVReadout( gal_[3][i], tmpAPV );
    galleryD_->Fill();
    copyAPVReadout( gal_[4][i], tmpAPV );
    galleryE_->Fill();
  }

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

void SiStripCMNAnalyzer::subtractFastLin( std::vector<int16_t>& in, std::vector<int16_t>& out,   
                                          std::vector<float>& offsets, std::vector<float>& slopes )
{
  std::vector<int16_t> tmp;  tmp.reserve(128);
  std::vector<int16_t>::iterator 
    strip( in.begin() ), 
    stripOut( out.begin() ),
    end(   in.end()   ),
    endAPV, high, low;

  while( strip < end ) {
    endAPV = strip+128; tmp.clear();
    tmp.insert(tmp.end(),strip,endAPV);
    const float offset = median(tmp);
    offsets.push_back( offset );
    
    low = strip;   high = strip+64;   tmp.clear(); 
    while( high < endAPV) tmp.push_back( *high++ - *low++ );
    const float slope = median(tmp)/64.;
    slopes.push_back( slope );    

    while (strip < endAPV) {
      *stripOut = *strip - (offset + slope*(65 - (endAPV-strip) ) );
      strip++; stripOut++;
    }

  }
}

void SiStripCMNAnalyzer::subtractSplitLin( std::vector<int16_t>& in, std::vector<int16_t>& out,   
                                           std::vector<float>& offsetsA, std::vector<float>& offsetsB,
                                           std::vector<float>& slopesA, std::vector<float>& slopesB )
{
  std::vector<int16_t> tmp;  tmp.reserve(64);
  std::vector<int16_t>::iterator 
    strip( in.begin() ), 
    stripOut( out.begin() ),
    end(   in.end()   ),
    endSegment, high, low;

  while( strip < end ) {
    endSegment = strip+64; tmp.clear();
    tmp.insert(tmp.end(),strip,endSegment);
    const float offsetA = median(tmp);
    offsetsA.push_back( offsetA );
    
    low = strip;   high = strip+32;   tmp.clear(); 
    while( high < endSegment) tmp.push_back( *high++ - *low++ );
    const float slopeA = median(tmp)/32.;
    slopesA.push_back( slopeA );    

    while (strip < endSegment) {
      *stripOut = *strip - (offsetA + slopeA*(33 - (endSegment-strip) ) );
      strip++; stripOut++;
    }

    endSegment = strip+64; tmp.clear();
    tmp.insert(tmp.end(),strip,endSegment);
    const float offsetB = median(tmp);
    offsetsB.push_back( offsetB );

    low = strip;   high = strip+32;   tmp.clear();
    while( high < endSegment) tmp.push_back( *high++ - *low++ );
    const float slopeB = median(tmp)/32.;
    slopesB.push_back( slopeB );

    while (strip < endSegment) {
      *stripOut = *strip - (offsetB + slopeB*(33 - (endSegment-strip) ) );
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

void SiStripCMNAnalyzer::subtractBLF( std::vector<int16_t>& in, std::vector<int16_t>& out,
                                      std::vector<std::vector<float> >& results, uint32_t detId)
{

  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId);
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId);

  std::vector<int16_t>::iterator fs,ls,os;
  std::vector<float>::iterator bs;
  os = out.begin();
  std::vector<float> adcs;
  adcs.reserve(128);

  results.clear();


  for( int16_t APV=0; APV< (int)in.size()/128; ++APV)
  {

    adcs.clear();
    // fill subset vector with all strips and their noises
    for (int16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip)
      adcs.push_back( (float)in[istrip] );
    

    DigiMap smoothedpoints;

    DigiMap consecpoints;
    DigiMapIter itConsecpoints, itSmoothedpoints;
    consecpoints.erase(consecpoints.begin(), consecpoints.end());
    smoothedpoints.erase(smoothedpoints.begin(), smoothedpoints.end());


    // Find the flat regions to make the smoothed map

    //============================= Height above local minimum ===============================                    
    float adcsLocalMinSubtracted[128];
    for(uint32_t istrip=0; istrip<128; ++istrip) {
      float localmin = 999.9;		
      for(uint16_t jstrip=std::max(0,(int)(istrip-nSmooth_/2)); jstrip<std::min(128,(int)(istrip+nSmooth_/2)); ++jstrip) {
        float nextvalue = adcs[jstrip];
        if(nextvalue < localmin) localmin=nextvalue;			
      }
      adcsLocalMinSubtracted[istrip] = adcs[istrip] - localmin;
    }
 
    //============================= Find regions with stable slopes ========================
    std::vector<uint16_t> nConsStrip;
    nConsStrip.clear();
  
    //Creating maps with all the neighborhood strip and putting in a nCosntStip vector how many we have
    uint16_t consecStrips=0;
    for(uint32_t istrip=0; istrip<128; ++istrip) {    
      int16_t adc = adcs[istrip]; 
    
      if( ( adcsLocalMinSubtracted[istrip] <   nSigmaNoiseDerTh_   * (float)noiseHandle->getNoise(istrip+APV*128,detNoiseRange) )
         /* && ( adc - median <  hitStripThreshold_ ) */ ) 
      {
        consecpoints.insert(consecpoints.end(), std::pair<uint16_t, int16_t >(istrip, adc));
        ++consecStrips;
      
      } else if (consecStrips >0){
        nConsStrip.push_back(consecStrips);
        consecStrips = 0;
      }
    }      		
  
    //to cope with the last flat region of the APV
    if(consecStrips >0) nConsStrip.push_back(consecStrips);
  
  
    //removing from the map the fist and last points in wide flat regions and erasing from the map too small regions
    itConsecpoints = consecpoints.begin();
    float MinSmoothValue=2000., MaxSmoothValue=0.;
    for(std::vector<uint16_t>::iterator itnConsStrip = nConsStrip.begin(); itnConsStrip < nConsStrip.end(); ++itnConsStrip){
    
      consecStrips = *itnConsStrip;
      if(consecStrips >=consecThreshold_){
        ++itConsecpoints;  //skipping first point
        uint16_t nFirstStrip = itConsecpoints->first;
        uint16_t nLastStrip;
        float smoothValue = 0.0;
        float stripCount =1;
        for(uint16_t n =0; n < consecStrips-2; ++n){
	  smoothValue += itConsecpoints->second;
	  if(stripCount == consecThreshold_){
	    smoothValue /= (float)stripCount;
	    nLastStrip = nFirstStrip + stripCount -1;				                    
	    smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
	    smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
	    if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
	    if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
	    nFirstStrip = nLastStrip+1;
	    smoothValue=0;
	    stripCount=0;
	  }
	  ++stripCount;
	  ++itConsecpoints;
        }
        ++itConsecpoints;  //and putting the pointer to the new seies of point 
      
        if(smoothValue>0){
	  --stripCount;
	  smoothValue /= (float)(stripCount);
	  nLastStrip = nFirstStrip + stripCount -1;
	  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
	  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
	  if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
	  if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
        }
      } else{
        for(int n =0; n< consecStrips ; ++n) ++itConsecpoints;
      }
    }

    // Now run the baseline follower over the smooth strips
    std::vector<float> baseline;
    baseline.clear();
    baseline.insert(baseline.begin(),128, 0);  
 
 
    //if not enough points
    if(smoothedpoints.size() >= minStripsToFit_){
     
      DigiMapIter itSmoothedpointsBegin, itSmoothedpointsEnd;
      itSmoothedpointsBegin = smoothedpoints.begin();
      itSmoothedpointsEnd = --(smoothedpoints.end());    
				
      uint16_t firstStripFlat = itSmoothedpointsBegin->first;
      uint16_t lastStripFlat = itSmoothedpointsEnd->first;
      int16_t firstStipFlatADC= itSmoothedpointsBegin->second;
      int16_t lastStipFlatADC= itSmoothedpointsEnd->second;
    
      //adding here the costant line at the extremities 
      baseline.erase(baseline.begin(), baseline.begin()+firstStripFlat);
      baseline.insert(baseline.begin(), firstStripFlat, firstStipFlatADC);
    
      baseline.erase(baseline.begin()+lastStripFlat, baseline.end());
      baseline.insert(baseline.end(), 128 - lastStripFlat, lastStipFlatADC);
    
    
      //IMPORTANT: the itSmoothedpointsEnd should be at least smaller than smoothedpoints.end() -1
      for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != itSmoothedpointsEnd; ++itSmoothedpoints){  
        DigiMapIter itSmoothedpointsNext = itSmoothedpoints;
        ++itSmoothedpointsNext;
        float strip1 = itSmoothedpoints->first;
        float strip2 = itSmoothedpointsNext->first;
        float adc1 = itSmoothedpoints->second;
        float adc2 = itSmoothedpointsNext->second;
      
        baseline[strip1] = adc1;
        baseline[strip2] = adc2;
        float m = (adc2 -adc1)/(strip2 -strip1);
        uint16_t itStrip = strip1 +1;
        float stripadc = adc1 + m; 
        while(itStrip < strip2){
	  baseline[itStrip] = stripadc;
	  ++itStrip;
	  stripadc+=m;
        }
      
      }
    
    }

    // store baseline vector
    results.push_back( baseline );
    //results.push_back( std::vector<float>() );
    //std::vector<float>::iterator derp;
   // for ( derp = baseline.begin(); derp != baseline.end(); ++derp)
   //   results[results.size()-1].push_back(*derp);

    // Now (finally) subtract the baseline
    fs = in.begin()+APV*128;
    ls = in.begin()+(APV+1)*128;
    bs = baseline.begin();
    while (fs < ls)
    {
      *os = *fs-*bs;
      ++os;
      ++fs;
      ++bs;
    }

  } // end loop over APVs

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


void SiStripCMNAnalyzer::fillClusterValues( edmNew::DetSetVector<SiStripCluster>& clusters, int * array, int APV )
{

   // Note: Integer array should be from a apvReadout_t struct and must contain at least 128 values!

   bool firstDet = true;
   int position = 0;
 
   // zero out the old points
   for(int i=0 ; i < 128; i++) array[i] = 0;

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
         array[position] = firststrip - APV*128;
         array[position+1] = laststrip - APV*128;
         position += 2;
       }
     }

   }
}

void SiStripCMNAnalyzer::copyAPVReadout( apvReadout_t & src, apvReadout_t & dest )
{
  for( int i = 0; i < 128; i++)
  {
    dest.BLFOffset[i] = src.BLFOffset[i];
    dest.adc[i] = src.adc[i];
    dest.clustersMedian[i] = src.clustersMedian[i];
    dest.clustersIterMed[i] = src.clustersIterMed[i];
    dest.clustersPer25[i] = src.clustersPer25[i];
    dest.clustersFastLin[i] = src.clustersFastLin[i];
    dest.clustersSplitLin[i] = src.clustersSplitLin[i];
    dest.clustersBLF[i] = src.clustersBLF[i];
  }
  dest.medianOffset = src.medianOffset;
  dest.iterMedOffset = src.iterMedOffset;
  dest.per25Offset = src.per25Offset;
  dest.fastLinOffset = src.fastLinOffset;
  dest.fastLinSlope = src.fastLinSlope;
  dest.splitLinOffsetA = src.splitLinOffsetA;
  dest.splitLinOffsetB = src.splitLinOffsetB;
  dest.splitLinSlopeA = src.splitLinSlopeA;
  dest.splitLinSlopeB = src.splitLinSlopeB;
  dest.event = src.event;
  dest.lumi = src.lumi;
  dest.run = src.run;
  dest.detID = src.detID;
  dest.apv = src.apv;
}  


bool SiStripCMNAnalyzer::galleryCuts( double median, double per25, int gal )
{

  switch (gal)
  {
    case 0:
      if ( ( 116 <= median && median < 140 ) && 
           ( 116 <= per25 && per25 < 140 ) ) 
         return true;
      return false;
    case 1:
      if ( ( 2 <= median && median < 100 ) && 
           ( 2 <= per25 && per25 < 100 ) ) 
         return true;
      return false;
    case 2:
      if ( ( median < 1 ) && 
           ( per25 < 1 ) ) 
         return true;
      return false;
    case 3:
      if ( ( 160 <= median ) && 
           ( 160 <= per25 ) )
         return true;
      return false;
    case 4:
      if ( ( 160 <= median ) && 
           ( 116 <= per25 && per25 < 140 ) ) 
         return true;
      return false;
  }
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripCMNAnalyzer);
