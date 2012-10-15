// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

//
// class declaration
//

class LowPtSpectrumAnalyzer : public edm::EDAnalyzer {
   public:
      explicit LowPtSpectrumAnalyzer(const edm::ParameterSet&);
      ~LowPtSpectrumAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag trackSrc_;
      edm::InputTag vertexSrc_;
      edm::InputTag centSrc_;
    
      double vertexCut_;
      double dzCut_;
      double chi2Cut_;

      CentralityProvider * centrality_;

      TH3F* tracks3D_;
      TH3F* tracks3D_c94_;
      TH3F* tracks3D_c97_;
      TH3F* tracks3D_c100_;


      TH1F* nevts_;
      TH1F* nevts_c94_;
      TH1F* nevts_c97_;
      TH1F* nevts_c100_;

      std::vector<double> etaBins, ptBins, centBins;

};

LowPtSpectrumAnalyzer::LowPtSpectrumAnalyzer(const edm::ParameterSet& iConfig):
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
centSrc_(iConfig.getParameter<edm::InputTag>("centSrc")),
vertexCut_(iConfig.getParameter<double>("vertexCut")),
dzCut_(iConfig.getParameter<double>("dzCut")),
chi2Cut_(iConfig.getParameter<double>("chi2Cut"))
{

  edm::Service<TFileService> fs;

  ptBins.clear();
  etaBins.clear();
  centBins.clear();

  const double small = 1e-3;
  double pt;
  for(pt =   0  ; pt <   23-small; pt +=  0.05) ptBins.push_back(pt); 

  for(double eta = -2.4; eta < 2.5; eta += 0.2)
    etaBins.push_back(eta);
  for(double c = 0.; c < 40.5; c += 1.0)
    centBins.push_back(c);
  


  tracks3D_ = fs->make<TH3F>("tracks3D", "Tracks - default;#eta;p_{T};centrality",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      centBins.size()-1, &centBins[0]); 
  tracks3D_c94_ = fs->make<TH3F>("tracks3D_c94", "Tracks - c94;#eta;p_{T};centrality",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      centBins.size()-1, &centBins[0]); 
  tracks3D_c97_ = fs->make<TH3F>("tracks3D_c97", "Tracks - c97;#eta;p_{T};centrality",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      centBins.size()-1, &centBins[0]); 
  tracks3D_c100_ = fs->make<TH3F>("tracks3D_c100", "Tracks - c100;#eta;p_{T};centrality",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      centBins.size()-1, &centBins[0]); 

  nevts_ = fs->make<TH1F>("nevts","Centrality Distribution - default", 40, 0, 40);
  nevts_c94_ = fs->make<TH1F>("nevts_c94","Centrality Distribution - 94", 40, 0, 40);
  nevts_c97_ = fs->make<TH1F>("nevts_c97","Centrality Distribution - 97", 40, 0, 40);
  nevts_c100_ = fs->make<TH1F>("nevts_c100","Centrality Distribution - 100", 40, 0, 40);

  // safety
  centrality_ = 0;
}


LowPtSpectrumAnalyzer::~LowPtSpectrumAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
LowPtSpectrumAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<std::vector<reco::Track> > tracks;
   iEvent.getByLabel(trackSrc_, tracks);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   Handle<reco::Centrality> cent;
   iEvent.getByLabel(centSrc_,cent);

   // Check that vertex abs(z) in range, and skip event if not
   if(vertex->size()<=0) return;
   if( fabs(vertex->begin()->z()) > vertexCut_ ) return;

   double vzErr = 0.;
   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   vtxPoint=vertex->begin()->position();
   vzErr=vertex->begin()->zError();

   // Get Centrality bin for the event

   if(!centrality_) centrality_ = new CentralityProvider(iSetup);

   centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event

   int bin = centrality_->getBin();
   nevts_->Fill(bin);

   double hfSum = cent->EtHFtowerSum();

   double c94low[40] = { 3621.35, 3253.01, 2934.71, 2650.77, 2395.67, 2167.02, 1955.71, 1762.4, 1585.1, 1422.41, 1273.68, 1137.43, 1012.88, 897.87, 793.592, 697.745, 611.215, 531.688, 460.817, 396.542, 339.401, 288.506, 243.63, 204.231, 170.035, 140.762, 115.854, 94.6972, 76.887, 62.2551, 50.0825, 40.2278, 32.0248, 25.4344, 19.9651, 15.4957, 11.5786, 6.09992, 6.03154, 0. };

   double c97low[40] = { 3632.33, 3271.89, 2960.67, 2681.84, 2430.14, 2204.26, 1995.9, 1805.25, 1628.74, 1466.07, 1318.75, 1182.32, 1057.25, 941.853, 837.035, 740.282, 651.636, 570.939, 497.955, 431.956, 372.691, 319.463, 272.161, 230.486, 193.83, 161.774, 134.629, 111.215, 91.2233, 74.5078, 60.6155, 49.0586, 39.6377, 31.7498, 25.3771, 20.0591, 15.6855, 11.8705, 6.03154, 0. };

   double c100low[40] = { 3643.67, 3290.47, 2986.31, 2712.91, 2465.38, 2243.22, 2037.14, 1848.25, 1673.34, 1512.3, 1364.51, 1228.63, 1103.12, 988.065, 881.488, 784.001, 694.358, 612.759, 537.393, 469.926, 408.04, 352.862, 303.186, 258.948, 220.051, 185.57, 155.669, 130.067, 107.895, 88.9554, 73.0311, 59.7535, 48.6524, 39.5531, 31.8866, 25.6635, 20.4462, 16.109, 12.3786, 0. };

   int bin94 =0;
   int bin97 =0;
   int bin100 =0;

   for( int i=0;i<39;i++) { if( hfSum < c94low[i] ) bin94 = i+1; }   
   for( int i=0;i<39;i++) { if( hfSum < c97low[i] ) bin97 = i+1; }   
   for( int i=0;i<39;i++) { if( hfSum < c100low[i] ) bin100 = i+1; }   

   nevts_c94_->Fill(0.5 + (double)bin94);
   nevts_c97_->Fill(0.5 + (double)bin97);
   nevts_c100_->Fill(0.5 + (double)bin100);

   

   // loop over tracks, fill histos

   std::vector<reco::Track>::const_iterator tk = tracks->begin();
   for ( ; tk != tracks->end(); ++tk )
   {
     //skip pixel tracks not within cuts
     bool isPixel = false;     
     if ( tk->numberOfValidHits() < 8 ) isPixel = true;
     bool isBad = false;
     double dz = tk->dz(vtxPoint); 
     double dzsigma = sqrt(tk->dzError()*tk->dzError()+vzErr*vzErr);
     if ( fabs(dz/dzsigma) > dzCut_ ) isBad = true;  
     if ( tk->normalizedChi2() > chi2Cut_ ) isBad = true;  

     if ( ! ( isBad && isPixel ) )
     {
       tracks3D_->Fill( tk->eta(), tk->pt(), 0.5+(double)bin );
       tracks3D_c94_->Fill( tk->eta(), tk->pt(), 0.5 + (double)bin94 );
       tracks3D_c97_->Fill( tk->eta(), tk->pt(), 0.5+(double)bin97 );
       tracks3D_c100_->Fill( tk->eta(), tk->pt(), 0.5+(double)bin100 );
     }
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
LowPtSpectrumAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LowPtSpectrumAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(LowPtSpectrumAnalyzer);

