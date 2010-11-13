#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include <math.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

int main(int argc, char* argv[]) 
{

  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  fwlite::TFileService fs = fwlite::TFileService("towers.root");
  TFileDirectory theDir = fs.mkdir("PlotCaloTowers");
  TH2F* energyEtaPhi_  = theDir.make<TH2F>("energyEtaPhi", "energyEtaPhi",    200,  -5.5, 5.5, 200, -3.14, 3.14 );

  std::vector<std::string> files;

  files.push_back( "dcap:///pnfs/cmsaf.mit.edu/t2bat/cms/store/user/vzhukova/HICorePhysics/FLOW_SKIMS_DATA2010_run150467/cff794fb2c241d2736cc3c6bcba83925/flow_skim_95_1_icC.root");

  fwlite::ChainEvent ev(files);

  unsigned int iEvent=0;
  for(ev.toBegin(); !ev.atEnd(); ++ev, ++iEvent){

    std::cout << "Processing event " << iEvent << "\n" ;

    edm::EventBase const & event = ev;

    edm::Handle<CaloTowerCollection>  towers;
    edm::InputTag towerLabel("towerMaker");
    event.getByLabel(towerLabel, towers);


    CaloTowerCollection::const_iterator tw = towers->begin();
    for ( ; tw != towers->end(); ++tw )
    {

      energyEtaPhi_->Fill( tw->eta(), tw->phi(), tw->energy() );

    }


   
  }


  return 0;

}

