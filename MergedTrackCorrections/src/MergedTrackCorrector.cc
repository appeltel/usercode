#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <fstream>
using namespace std;


class MergedTrackCorrector
{
  public:
   
   MergedTrackCorrector();
   ~MergedTrackCorrector();

   void loadEffic();
   void loadEffic(const char *);
   void loadFakes();
   void loadFakes(const char *);

//   float weight( reco::Track &, int );
//   float fakeRate( reco::Track &, int );
//   float effic( reco::Track &, int );

   float weight( float , float , int ); 
   float fakeRate( float , float , int ); 
   float effic( float , float , int ); 

  private:
 
   float effic_[3][12][10][4];
   float fake_[3][12][10][4];

};

MergedTrackCorrector::MergedTrackCorrector(){}
MergedTrackCorrector::~MergedTrackCorrector(){}

void MergedTrackCorrector::loadEffic()
{
  loadEffic("Appeltel/MergedTrackCorrections/data/CorrTable_Bass_100k_effic.par");
}

void MergedTrackCorrector::loadEffic( const char * filename )
{

  edm::FileInPath fileInPath(filename);
  ifstream inFile(fileInPath.fullPath().c_str());

  while(inFile.eof() == false)
  {

    int pt, eta, cent;

    inFile >> pt;
    inFile >> eta;
    inFile >> cent;

    for( int par = 0; par<4; par++ )
      inFile >> effic_[pt][eta][cent][par];
  }

  inFile.close();
}

void MergedTrackCorrector::loadFakes()
{
  loadFakes("Appeltel/MergedTrackCorrections/data/CorrTable_Bass_100k_fakes.par");
}

void MergedTrackCorrector::loadFakes( const char * filename )
{

  edm::FileInPath fileInPath(filename);
  ifstream inFile(fileInPath.fullPath().c_str());

  while(inFile.eof() == false)
  {

    int pt, eta, cent;

    inFile >> pt;
    inFile >> eta;
    inFile >> cent;

    for( int par = 0; par<4; par++ )
      inFile >> fake_[pt][eta][cent][par];
  }

  inFile.close();
}



float MergedTrackCorrector::weight( float pt, float eta, int cbin )
{
   return  ( 1. - fakeRate(pt, eta, cbin) )/ effic( pt, eta, cbin );
}

float MergedTrackCorrector::effic( float pt, float eta, int cbin )
{
   int etabin = (int) ( (eta+3.)*2. );
   int ptbin;
   cbin = cbin / 4;
   if( pt < 1.5 )   
     ptbin = 0;
   if( pt >= 1.5 && pt < 1.8 )
     ptbin = 1;
   if( pt >= 1.8)
     ptbin = 2;

   return effic_[ptbin][etabin][cbin][0]*pt*pt*pt + 
          effic_[ptbin][etabin][cbin][1]*pt*pt + 
          effic_[ptbin][etabin][cbin][2]*pt + 
          effic_[ptbin][etabin][cbin][3]; 
}    
  


float MergedTrackCorrector::fakeRate( float pt, float eta, int cbin )
{
   int etabin = (int) ( (eta+3.)*2. );
   int ptbin;
   cbin = cbin / 4;
   if( pt < 1.5 )   
     ptbin = 0;
   if( pt >= 1.5 && pt < 1.8 )
     ptbin = 1;
   if( pt >= 1.8)
     ptbin = 2;

   return fake_[ptbin][etabin][cbin][0]*pt*pt*pt + 
          fake_[ptbin][etabin][cbin][1]*pt*pt + 
          fake_[ptbin][etabin][cbin][2]*pt + 
          fake_[ptbin][etabin][cbin][3]; 
}    
   
  

