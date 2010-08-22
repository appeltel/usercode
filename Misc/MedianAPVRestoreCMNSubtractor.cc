#include "RecoLocalTracker/SiStripZeroSuppression/interface/MedianAPVRestoreCMNSubtractor.h"

void MedianAPVRestoreCMNSubtractor::subtract(const uint32_t& detId,std::vector<int16_t>& digis) {subtract_(detId,digis);}
void MedianAPVRestoreCMNSubtractor::subtract(const uint32_t& detId,std::vector<float>& digis) {subtract_(detId,digis);}

template<typename T> 
inline
void MedianAPVRestoreCMNSubtractor::
subtract_(const uint32_t& detId,std::vector<T>& digis){
  
  std::vector<T> tmp;  tmp.reserve(128);  
  typename std::vector<T>::iterator  
    strip( digis.begin() ), tmpiter,
    end(   digis.end()   ),
    endAPV;
  
  _vmedians.clear();
  
  while( strip < end ) {
    endAPV = strip+128; tmp.clear();
    tmp.insert(tmp.end(),strip,endAPV);
    const float offset = median(tmp);

    int zeroCount = 0;
    for ( tmpiter = tmp.begin(); tmpiter < tmp.end(); tmpiter++ )
      if ( (int) *tmpiter < 1 ) zeroCount++;

    
    if ( zeroCount > 64 )
    {
      _vmedians.push_back(std::make_pair<short,float>((strip-digis.begin())/128,-10));
      int counter = 0;
      while (strip < endAPV) {
        *strip = static_cast<T>(150);
        if (counter == 0) *strip = static_cast<T>(0);
        if (counter == 127) *strip = static_cast<T>(0);
        counter++;
        strip++;
      }
    } else {
      _vmedians.push_back(std::make_pair<short,float>((strip-digis.begin())/128,offset));
      while (strip < endAPV) {
        *strip = static_cast<T>(*strip-offset);
        strip++;
      }
    }

  }
}
