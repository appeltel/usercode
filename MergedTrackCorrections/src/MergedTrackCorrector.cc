#include "Appeltel/MergedTrackCorrections/interface/MergedTrackCorrector.h"

#include <fstream>
#include <strstream>
#include <iostream>

using namespace std;

/********************************************************
 
Fit functions used are

- Efficiency

- [0]/(1+exp([1]*x+[2])) + [3]*x^2
- pol2
- pol2

- Fake rate

- [0]/(1+exp([1]*x+[2])) + [3]*x^2
- pol2
- pol3


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


MergedTrackCorrector::MergedTrackCorrector(){}
MergedTrackCorrector::~MergedTrackCorrector(){}

void MergedTrackCorrector::loadEffic()
{
  loadEffic("Appeltel/MergedTrackCorrections/data/EFF_validation3D_bass_hiBadgerTracks_v2.par");
}

void MergedTrackCorrector::loadEffic( const char * filename )
{

  edm::FileInPath fileInPath(filename);
  
  ifstream inFile(fileInPath.fullPath().c_str());

  char buffer1[2048];


  while(inFile.getline(buffer1, 2048))
  {

    int pt, eta, cent;

  istrstream ostr1(buffer1, 2048);

    ostr1 >> pt;
    ostr1 >> eta;
    ostr1 >> cent;

    int parmax = (pt == 0) ? 4 : 3;
    for( int par = 0; par<parmax; par++ )
    {
      ostr1 >> effic_[pt][eta][cent][par];
    }
  }


  inFile.close();


}

void MergedTrackCorrector::loadFakes()
{
  loadFakes("Appeltel/MergedTrackCorrections/data/FAK_validation3D_bass_hiBadgerTracks_v2.par");
}

void MergedTrackCorrector::loadFakes( const char * filename )
{

  edm::FileInPath fileInPath(filename);
  ifstream inFile(fileInPath.fullPath().c_str());
  char buffer1[2048];


  while(inFile.getline(buffer1, 2048))
  {

    int pt, eta, cent;
  istrstream ostr1(buffer1, 2048);

    ostr1 >> pt;
    ostr1 >> eta;
    ostr1 >> cent;

    int parmax = (pt==1) ? 3 : 4;
    for( int par = 0; par<parmax; par++ )
      ostr1 >> fake_[pt][eta][cent][par];
  }

  inFile.close();
}



float MergedTrackCorrector::weight( float pt, float eta, int cbin )
{
   return  ( 1. - fakeRate(pt, eta, cbin) )/ effic( pt, eta, cbin );
}

float MergedTrackCorrector::fakeRate( float pt, float eta, int centbin )
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

   int ptbin = 0;
   if( pt < 1.5 )   
     ptbin = 0;
   if( pt >= 1.5 && pt < 1.8 )
     ptbin = 1;
   if( pt >= 1.8)
     ptbin = 2;

   if ( ptbin == 2)
     return fake_[ptbin][etabin][cbin][3]*pt*pt*pt + 
            fake_[ptbin][etabin][cbin][2]*pt*pt + 
            fake_[ptbin][etabin][cbin][1]*pt + 
            fake_[ptbin][etabin][cbin][0]; 

   if ( ptbin == 1)
     return fake_[ptbin][etabin][cbin][2]*pt*pt + 
            fake_[ptbin][etabin][cbin][1]*pt + 
            fake_[ptbin][etabin][cbin][0];

   if ( ptbin == 0)
   {
      float par0 = fake_[ptbin][etabin][cbin][0]; 
      float par1 = fake_[ptbin][etabin][cbin][1]; 
      float par2 = fake_[ptbin][etabin][cbin][2]; 
      float par3 = fake_[ptbin][etabin][cbin][3]; 
      
      return par0 / ( 1. + exp(par1*pt+par2)) + par3*pt*pt;
   }

   return 0.0;
}    
  

float MergedTrackCorrector::effic( float pt, float eta, int centbin )
{

   int etabin = 0;
   if (  eta >= -2.4 && eta < -1.0 ) 
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

   int ptbin = 0;
   if( pt < 1.5 )   
     ptbin = 0;
   if( pt >= 1.5 && pt < 1.8 )
     ptbin = 1;
   if( pt >= 1.8)
     ptbin = 2;


   if ( ptbin >0)
     return effic_[ptbin][etabin][cbin][2]*pt*pt + 
            effic_[ptbin][etabin][cbin][1]*pt + 
            effic_[ptbin][etabin][cbin][0];

   if ( ptbin == 0)
   {
      float par0 = effic_[ptbin][etabin][cbin][0]; 
      float par1 = effic_[ptbin][etabin][cbin][1]; 
      float par2 = effic_[ptbin][etabin][cbin][2]; 
      float par3 = effic_[ptbin][etabin][cbin][3]; 
      
      return par0 / ( 1. + exp(par1*pt+par2)) + par3*pt*pt;
   }

   return 1.0;
}    

   
  

