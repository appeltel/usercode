#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"


int main(int argc, char* argv[]) 
{

  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

//  TFile* inFile = TFile::Open( "file:reco.root" );
//  TFile* inFile = TFile::Open( "file:tier0.root" );
  TFile* inFile = TFile::Open( "file:tier0_fix.root" );

  unsigned int iEvent=0;
  fwlite::Event ev(inFile);
  for(ev.toBegin(); !ev.atEnd(); ++ev, ++iEvent){
    edm::EventBase const & event = ev;

    if( iEvent==2 ) break;

    edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
    edm::InputTag clusLabel("siStripClusters");
    event.getByLabel(clusLabel, clusters);

edmNew::DetSetVector<SiStripCluster>::const_iterator it = clusters->begin();
for ( ; it != clusters->end(); ++it )
{
  for ( edmNew::DetSet<SiStripCluster>::const_iterator clus = it->begin(); clus != it->end(); ++clus)
  {
     int firststrip = clus->firstStrip();
     int laststrip = firststrip + clus->amplitudes().size();
     
     if ( clus->geographicalId() == 436233365 )
     {


       std::cout << "Cluster Found from " << firststrip % 128 << " APV: " << firststrip / 128  
                 << " to " << laststrip % 128 << " APV: " << laststrip / 128 << "\n";
       int derp = 1;
       for( std::vector<uint8_t>::const_iterator ii = clus->amplitudes().begin();
            ii != clus->amplitudes().end(); ++ii)
       { std::cout << derp << " " <<(int) *ii << "\n"; derp++; }
       std::cout << "\n\n";

     }
  }
}

}


  inFile->Close();

  return 0;

}



