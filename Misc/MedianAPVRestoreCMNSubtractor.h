#ifndef RECOLOCALTRACKER_SISTRIPZEROSUPPRESSION_SISTRIPMEDIANAPVRESTORECOMMONMODENOISESUBTRACTION_H
#define RECOLOCALTRACKER_SISTRIPZEROSUPPRESSION_SISTRIPMEDIANAPVRESTORECOMMONMODENOISESUBTRACTION_H
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripCommonModeNoiseSubtractor.h"

class MedianAPVRestoreCMNSubtractor : public SiStripCommonModeNoiseSubtractor {
  
  friend class SiStripRawProcessingFactory;
  
 public:
  
  void subtract(const uint32_t&,std::vector<int16_t>&);
  void subtract(const uint32_t&,std::vector<float>&);
  
 private:
  
  template<typename T> void subtract_(const uint32_t&,std::vector<T>&);
  MedianAPVRestoreCMNSubtractor(){};
  
};
#endif
