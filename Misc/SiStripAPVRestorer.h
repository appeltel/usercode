#ifndef RECOLOCALTRACKER_SISTRIPZEROSUPPRESSION_SISTRIPAPVRESTORER_H
#define RECOLOCALTRACKER_SISTRIPZEROSUPPRESSION_SISTRIPAPVRESTORER_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"


#include <vector>
//#include <algorithm>
#include <stdint.h>

class SiStripQuality;

class SiStripAPVRestorer {

 friend class SiStripRawProcessingFactory;


 public:
  
  SiStripAPVRestorer( double rt)
  : restoreThreshold_(rt)
  {};
  ~SiStripAPVRestorer() {};
  void init(const edm::EventSetup& es);
  void inspect(const uint32_t&, std::vector<int16_t>&);
  void inspect(const uint32_t&, std::vector<float>&);
  void restore(std::vector<int16_t>&);
  void restore(std::vector<float>&);
  
  void fixAPVsCM(edm::DetSet<SiStripProcessedRawDigi>& );
  
 private:

  template<typename T >void inspect_(const uint32_t&,std::vector<T>&);
  template<typename T >void restore_(std::vector<T>&);

  
  std::vector<bool> apvKills;

  edm::ESHandle<SiStripQuality> qualityHandle;
  uint32_t  quality_cache_id;
  
  double restoreThreshold_;
};

#endif
