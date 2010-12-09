#include "Appeltel/MergedTrackCorrections/interface/UnmergedTrackCorrector.h"

#include <fstream>
#include <strstream>
#include <iostream>

using namespace std;

/********************************************************
 
Fit functions used are

- Efficiency
- Fake rate

[0]*pow(x,[3])/(1+exp([1]*(x+[2]))) + [4]*pow(x,[5])

Eta and Centrality binning are done in the following way,

- 2.4~-1.0 : bin width of 0.2
- 1.0~1.0: bin width of 0.4
- 1.0~2.4: bin width of 0.2

==> 19 bins (0~18)

- 0~10%
- 10~30%
- 30~100%

==> 3 bins (0~3)


********************************************************/


UnmergedTrackCorrector::UnmergedTrackCorrector(){}
UnmergedTrackCorrector::~UnmergedTrackCorrector(){}

void UnmergedTrackCorrector::loadMap()
{
  loadMap("Appeltel/MergedTrackCorrections/data/MAP_validation3D_bass_hiGoodTracks.par");
}

void UnmergedTrackCorrector::loadMap( const char * filename )
{

  edm::FileInPath fileInPath(filename);
  
  ifstream inFile(fileInPath.fullPath().c_str());

  char buffer1[2048];


  while(inFile.getline(buffer1, 2048))
  {

    int ef, dummy, eta, cent;

  istrstream ostr1(buffer1, 2048);

    ostr1 >> ef;
    ostr1 >> dummy;
    ostr1 >> eta;
    ostr1 >> cent;

    for( int par = 0; par<6; par++ )
    {
      ostr1 >> m_[ef][eta][cent][par];
    }
  }


  inFile.close();


}



float UnmergedTrackCorrector::weight( float pt, float eta, int cbin )
{
   return  ( 1. - fakeRate(pt, eta, cbin) )/ effic( pt, eta, cbin );
}

float UnmergedTrackCorrector::fakeRate( float pt, float eta, int centbin )
{

   int etabin = 0;
   if ( eta >= -2.4 && eta < -1.0 ) 
     etabin = floor((eta+2.4)*5);
   if ( eta >= -1.0 && eta < 1.0 ) 
     etabin = floor( (eta+1.0)*2.5 + 7 );
   if ( eta >= 1.0 && eta < 2.4)
     etabin = floor( (eta-1.0)*5 + 12 );

   int cbin = 0;
   if ( centbin < 4)
     cbin = 0;
   if ( centbin >= 4 && centbin < 12 )
     cbin = 1;
   if ( centbin >= 12)
     cbin = 2;

   double par0 = m_[1][etabin][cbin][0]; 
   double par1 = m_[1][etabin][cbin][1]; 
   double par2 = m_[1][etabin][cbin][2]; 
   double par3 = m_[1][etabin][cbin][3]; 
   double par4 = m_[1][etabin][cbin][4]; 
   double par5 = m_[1][etabin][cbin][5]; 
 
   return par0 * pow(pt,par3)/(1.+exp(par1*(pt+par2))) + par4*pow(pt,par5);
      
}    
  
float UnmergedTrackCorrector::effic( float pt, float eta, int centbin )
{

   int etabin = 0;
   if ( eta >= -2.4 && eta < -1.0 ) 
     etabin = floor((eta+2.4)*5);
   if ( eta >= -1.0 && eta < 1.0 ) 
     etabin = floor( (eta+1.0)*2.5 + 7 );
   if ( eta >= 1.0 && eta < 2.4)
     etabin = floor( (eta-1.0)*5 + 12 );

   int cbin = 0;
   if ( centbin < 4)
     cbin = 0;
   if ( centbin >= 4 && centbin < 12 )
     cbin = 1;
   if ( centbin >= 12)
     cbin = 2;

   double par0 = m_[0][etabin][cbin][0]; 
   double par1 = m_[0][etabin][cbin][1]; 
   double par2 = m_[0][etabin][cbin][2]; 
   double par3 = m_[0][etabin][cbin][3]; 
   double par4 = m_[0][etabin][cbin][4]; 
   double par5 = m_[0][etabin][cbin][5]; 
 
   return par0 * pow(pt,par3)/(1.+exp(par1*(pt+par2))) + par4*pow(pt,par5);
      
}    
  
