//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
//
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>


class RpPbMCAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RpPbMCAnalyzer(const edm::ParameterSet&);
      ~RpPbMCAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void printWeirdParticle( const reco::GenParticle &);
      void printRoguePion( const reco::GenParticle &, double vdist);

      // ----------member data ---------------------------

      TH1F* chargedHadrons_;
      TH1F* pions_;
      TH1F* kaons_;
      TH1F* protons_;
      TH1F* antiprotons_;
      TH1F* chargedHadrons_finer_;
      TH1F* pions_finer_;
      TH1F* kaons_finer_;
      TH1F* protons_finer_;
      TH1F* antiprotons_finer_;

      TH1F* events_;

      TH1F* vtxNtrk_;
      TH1F* vtxN_;
      TH1F* distPV_;

      int nevt_;

      edm::InputTag genSrc_;
      edm::InputTag vertexSrc_;
      double etaMax_;
      double rapidityShift_;
      double vzMax_;
      bool printDebug_;

};

RpPbMCAnalyzer::RpPbMCAnalyzer(const edm::ParameterSet& iConfig):
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
printDebug_(iConfig.getParameter<bool>("printDebug")),
etaMax_(iConfig.getParameter<double>("etaMax")),
vzMax_(iConfig.getParameter<double>("vzMax")),
rapidityShift_(iConfig.getParameter<double>("rapidityShift"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;

   std::vector<Double_t> ptBins; std::vector<Double_t> ptBins_finer;
   const Double_t small = 1e-3;
   Double_t pt; Double_t pt_finer;

   for(pt =   0.2  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt);
   for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); //12 bins
   for(pt =   2.4; pt <   7.2-small; pt +=  0.4 ) ptBins.push_back(pt); //12 bins
   for(pt =   7.2; pt <  16.8-small; pt +=  1.6 ) ptBins.push_back(pt);//it was 3.2
   for(pt =  16.8; pt <  30.0-small; pt +=  6.6 ) ptBins.push_back(pt);
   ptBins.push_back(30.0);

   //finer bins
   for(pt_finer =   0.2  ; pt_finer <   1.2-small; pt_finer +=  0.005) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   1.2; pt_finer <   2.4-small; pt_finer +=  0.01 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   2.4; pt_finer <   7.2-small; pt_finer +=  0.04 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =   7.2; pt_finer <  16.8-small; pt_finer +=  0.16 ) ptBins_finer.push_back(pt_finer);
   for(pt_finer =  16.8; pt_finer <  30.0-small; pt_finer +=  0.66 ) ptBins_finer.push_back(pt_finer);
   ptBins_finer.push_back(30.0);

    events_ = fs->make<TH1F>("events","",10,0,10);

    vtxNtrk_ = fs->make<TH1F>("vtxNtrk","Tracks per vertex",200,0,200);
    vtxN_ = fs->make<TH1F>("vtxN","Vertices per event",20,0,20);

    distPV_ = fs->make<TH1F>("distPV_","3d distance from best primary vertex",100,0,20);

    chargedHadrons_ = fs->make<TH1F>(Form("chargedHadrons"),
      Form("Charged Hadron Spectrum, |#eta| < %f ",etaMax_),ptBins.size()-1, &ptBins[0]); 
    pions_ = fs->make<TH1F>(Form("pions"),
      Form("Charged Pion Spectrum, |#eta| < %f ",etaMax_),ptBins.size()-1, &ptBins[0]); 
    kaons_ = fs->make<TH1F>(Form("kaons"),
      Form("Charged Kaon Spectrum, |#eta| < %f ",etaMax_),ptBins.size()-1, &ptBins[0]); 
    protons_ = fs->make<TH1F>(Form("protons"),
      Form("Proton Spectrum, |#eta| < %f ",etaMax_),ptBins.size()-1, &ptBins[0]); 
    antiprotons_ = fs->make<TH1F>(Form("antiprotons"),
      Form("Antiproton Spectrum, |#eta| < %f ",etaMax_),ptBins.size()-1, &ptBins[0]); 

    chargedHadrons_finer_ = fs->make<TH1F>(Form("chargedHadrons_finer"),
      Form("Charged Hadron Spectrum, |#eta| < %f ",etaMax_),ptBins_finer.size()-1, &ptBins_finer[0]); 
    pions_finer_ = fs->make<TH1F>(Form("pions_finer"),
      Form("Charged Pion Spectrum, |#eta| < %f ",etaMax_),ptBins_finer.size()-1, &ptBins_finer[0]); 
    kaons_finer_ = fs->make<TH1F>(Form("kaons_finer"),
      Form("Charged Kaon Spectrum, |#eta| < %f ",etaMax_),ptBins_finer.size()-1, &ptBins_finer[0]); 
    protons_finer_ = fs->make<TH1F>(Form("protons_finer"),
      Form("Proton Spectrum, |#eta| < %f ",etaMax_),ptBins_finer.size()-1, &ptBins_finer[0]); 
    antiprotons_finer_ = fs->make<TH1F>(Form("antiprotons_finer"),
      Form("Antiproton Spectrum, |#eta| < %f ",etaMax_),ptBins_finer.size()-1, &ptBins_finer[0]); 

    nevt_ = 0;
}

RpPbMCAnalyzer::~RpPbMCAnalyzer()
{

  std::cout << "Total passed events: " << nevt_ << "\n\n";

}

void
RpPbMCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   Handle<reco::GenParticleCollection> gcol;
   iEvent.getByLabel(genSrc_, gcol);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   double vz =0.0, vx=0.0, vy=0.0;

   if(vertex->size()>0) {
     vz=vertex->begin()->z();
     vx=vertex->begin()->x();
     vy=vertex->begin()->y();
     } else { return; }

   unsigned int ntrk = 0;
   int vcount = 1;
   int vp = 1;

   for( const auto & vi : * vertex )
   {
     if( ntrk < vi.tracksSize() )
     { 
        ntrk = vi.tracksSize();
        vx = vi.x();
        vy = vi.y();
        vz = vi.z();
        vp = vcount;
     }
     vtxNtrk_->Fill(ntrk);
     if ( printDebug_ )
       std::cout << "reco::Vertex " << vcount << "  ntrk= " << vi.tracksSize()
                 << " x= " << vi.x()
                 << " y= " << vi.y()
                 << " z= " << vi.z() << std::endl;
     vcount++;
   }

   vtxN_->Fill(vertex->size() );

   if ( printDebug_ )
   {
     std::cout << "Total reconstructed vertices: " << vertex->size() << " using vertex " << vp << std::endl; 
     std::cout << "Primary vertex x= " << vx << " y= " << vy 
               << " z= " << vz << " ntrk= " << ntrk << std::endl;
   }

   // skip displaced vertices

   if( fabs(vz) > vzMax_ ) return;

   nevt_++;

   events_->Fill(0.5); 

   for( const auto & particle : * gcol )
   {
     if (printDebug_) printWeirdParticle( particle );

     if ( particle.status() == 1 && 
          particle.charge() != 0 &&
          fabs(particle.eta() - rapidityShift_ ) <= etaMax_ )
     {
       chargedHadrons_->Fill( particle.pt() );
       chargedHadrons_finer_->Fill( particle.pt() );

       double vdist = sqrt( (particle.vx() - vx)* (particle.vx() - vx)
                           + (particle.vy() - vy)* (particle.vy() - vy)
                           + (particle.vz() - vz)* (particle.vz() - vz) );
      
       distPV_->Fill(vdist);
 
       if ( abs(particle.pdgId()) == 211 )
       {
         pions_->Fill( particle.pt() );
         pions_finer_->Fill( particle.pt() );
         // find rogue pions
         if ( vdist > 1. && printDebug_) printRoguePion(particle,vdist);
       } 
       if ( abs(particle.pdgId()) == 321 )
         kaons_->Fill( particle.pt() );
         kaons_finer_->Fill( particle.pt() );
       if ( particle.pdgId() == 2212 )
         protons_->Fill( particle.pt() );
         protons_finer_->Fill( particle.pt() );
       if ( particle.pdgId() == -2212 )
         antiprotons_->Fill( particle.pt() );
         antiprotons_finer_->Fill( particle.pt() );
     }
   }
}

void 
RpPbMCAnalyzer::printWeirdParticle( const reco::GenParticle & particle)
{
   if (particle.status() > 3 ||
         ( particle.status() == 1 &&
           (  ( abs(particle.pdgId()) > 400 &&  abs(particle.pdgId()) < 600 )
             || ( abs(particle.pdgId()) > 4000 && abs(particle.pdgId()) < 6000 )
           )
         )
      )
   {
       std::cout << "  Stable B/D status= " << particle.status()
                 << " pdgId= " << particle.pdgId()
                 << " pt= " << particle.pt()
                 << " eta= " << particle.eta()
                 << " phi= " << particle.phi()
                 << std::endl;
   }
}

void 
RpPbMCAnalyzer::printRoguePion( const reco::GenParticle & particle, double vdist)
{
  std::cout << "  Rogue pion px= " << particle.vx()
            << " py= " << particle.vy()
            << " pz= " << particle.vz()
            << " vdist= " << vdist
            << " pt= " << particle.pt()
            << " eta= " << particle.eta()
            << " phi= " << particle.phi()
            << std::endl;
  if( particle.numberOfMothers() > 0 )
    std::cout << "  Mother pdgId= " << particle.mother()->pdgId() << std::endl;
}

void
RpPbMCAnalyzer::beginJob()
{
}

void
RpPbMCAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(RpPbMCAnalyzer);
