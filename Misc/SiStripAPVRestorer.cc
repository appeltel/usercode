#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripAPVRestorer.h"

#include <cmath>

#include <iostream>

void SiStripAPVRestorer::inspect( const uint32_t& detId,std::vector<int16_t>& digis) {inspect_(detId,digis);}
void SiStripAPVRestorer::inspect( const uint32_t& detId,std::vector<float>& digis) {inspect_(detId,digis);}

void SiStripAPVRestorer::restore( std::vector<int16_t>& digis) {restore_(digis);}
void SiStripAPVRestorer::restore( std::vector<float>& digis) {restore_(digis);}

void SiStripAPVRestorer::init(const edm::EventSetup& es){
  uint32_t q_cache_id = es.get<SiStripQualityRcd>().cacheIdentifier();

  if(q_cache_id != quality_cache_id) {
    es.get<SiStripQualityRcd>().get( qualityHandle );
    quality_cache_id = q_cache_id;
  }
}


template<typename T>
inline
void SiStripAPVRestorer::
inspect_(const uint32_t& detId,std::vector<T>& digis){

  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId);

  typename std::vector<T>::iterator fs;

  apvKills.clear();

  for( uint16_t APV=0; APV< digis.size()/128; ++APV)
  {
    int zeroCount = 0, qualityCount = 0; 
    for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip)
    {
      fs = digis.begin() + istrip;
      if ( !qualityHandle->IsStripBad(detQualityRange,istrip) )
      {
        qualityCount++; 
        if ( (int) *fs < 1 ) zeroCount++;
      }
    }

    if( zeroCount > restoreThreshold_ * qualityCount ) 
      apvKills.push_back( true );
    else 
      apvKills.push_back( false );
    
  }
}

template<typename T>
inline
void SiStripAPVRestorer::
restore_( std::vector<T>& digis ){

  typename std::vector<T>::iterator  
  strip( digis.begin() ), 
  endAPV;

  for( uint16_t APV=0; APV< digis.size()/128; ++APV)
  {
    endAPV = digis.begin() + (APV+1)*128;
    if ( *( apvKills.begin() + APV ) )
    {
      int counter = 0;
      while (strip < endAPV) {
        *strip = static_cast<T>(150);
        if (counter == 0) *strip = static_cast<T>(0);
        if (counter == 127) *strip = static_cast<T>(0);
        counter++;
        strip++;
      }
    }
  }
}

void SiStripAPVRestorer::fixAPVsCM(edm::DetSet<SiStripProcessedRawDigi>& cmdigis )
{

  // cmdigis should be the same size as apvKills
  // otherwise something pathological has happened and we do nothing
  if ( cmdigis.size() != apvKills.size() ) return;

  edm::DetSet<SiStripProcessedRawDigi>::iterator cm_iter = cmdigis.begin();
  std::vector<bool>::const_iterator apvk_iter = apvKills.begin();

  // No way to change the adc value of a SiStripProcessedRawDigi
  // so we just extract the values, clear the DetSet, and
  // replace with the proper values.

  std::vector<float> cmvalues;
  for( ; cm_iter != cmdigis.end(); ++cm_iter  ) cmvalues.push_back( (*cm_iter).adc() );
  cmdigis.clear();

  std::vector<float>::const_iterator cmv_iter = cmvalues.begin();
  while( apvk_iter != apvKills.end() )
  {
    if( *apvk_iter)
      cmdigis.push_back( SiStripProcessedRawDigi( -999.) );
    else
      cmdigis.push_back( SiStripProcessedRawDigi( *cmv_iter ) );
    apvk_iter++;
    cmv_iter++;
  }

}
  
