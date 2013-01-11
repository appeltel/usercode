//
// Original Author:  Eric Appelt
//         Created:  Fri Nov 19 16:55:16 CST 2010
//
//

#include <memory>
#include <algorithm>

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
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>


class RpPbTrackingAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RpPbTrackingAnalyzer(const edm::ParameterSet&);
      ~RpPbTrackingAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> evtPerf_;
      std::map<std::string,TH2F*> evtPerf2D_;
      std::map<std::string,TH1F*> trkPerf_;
      std::map<std::string,TH2F*> trkPerf2D_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;

      TH1F* events_;
      TH1F* vertices_;
      TH1F* tracks_;
      TH1F* trackseta_;

      int nevt_;
      int ntrack_;
      int nvertex_;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      double etaMin_;
      double etaMax_;

};

RpPbTrackingAnalyzer::RpPbTrackingAnalyzer(const edm::ParameterSet& iConfig):
nevt_(0),
ntrack_(0),
nvertex_(0),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
}

RpPbTrackingAnalyzer::~RpPbTrackingAnalyzer()
{
}

void
RpPbTrackingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::TrackCollection> tcol;
   iEvent.getByLabel(trackSrc_, tcol);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);
  
   std::vector<reco::Vertex> vsorted = *vertex;
   // sort the vertcies by number of tracks in descending order
   std::sort( vsorted.begin(), vsorted.end(), 
              []( reco::Vertex a, reco::Vertex b) 
   {
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   });

  
   nevt_++;
   events_->Fill(0.5); 
   evtPerf_["Nvtx"]->Fill(vsorted.size());
   evtPerf_["Ntrk"]->Fill(tcol->size());

   int vcount = 0; 
   for( const auto & vi :  vsorted )
   {
     vtxPerf_["Ntrk"]->Fill(vi.tracksSize());
     vtxPerf_["x"]->Fill(vi.x());
     vtxPerf_["y"]->Fill(vi.y());
     vtxPerf_["z"]->Fill(vi.z());
     vtxPerf2D_["Ntrk2D"]->Fill(vcount,vi.tracksSize());
     vertices_->Fill(0.5);
     vcount++;
     nvertex_++;
   }

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
  
   // use vertex w most tracks as primary vertex
   if( vsorted.size() != 0 )
   {
     vtxPoint=vsorted.begin()->position();
     vzErr=vsorted.begin()->zError();
     vxErr=vsorted.begin()->xError();
     vyErr=vsorted.begin()->yError();
   }

   for( const auto & track : * tcol )
   {

     trkPerf_["eta"]->Fill( track.eta() );
     trkPerf2D_["etaphi"]->Fill( track.eta(), track.phi() );
     tracks_->Fill(0.5);

     if( track.eta() <= etaMax_ && track.eta() >= etaMin_ )
     {

       double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
       dxy = track.dxy(vtxPoint);
       dz = track.dz(vtxPoint);
       dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
       dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
       
       trackseta_->Fill(0.5);
       trkPerf_["Nhit"]->Fill(track.numberOfValidHits()); 
       trkPerf_["pt"]->Fill(track.pt()); 
       trkPerf_["phi"]->Fill(track.phi()); 
       trkPerf_["dxyErr"]->Fill(dxy/dxysigma); 
       trkPerf_["dzErr"]->Fill(dz/dzsigma); 
       trkPerf_["chi2"]->Fill(track.normalizedChi2());
       trkPerf_["pterr"]->Fill(track.ptError() / track.pt() );
       ntrack_++;
     }

   }
}


void
RpPbTrackingAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  events_ = fs->make<TH1F>("events","",1,0,1);
  tracks_ = fs->make<TH1F>("tracks","",1,0,1);
  trackseta_ = fs->make<TH1F>("trackseta","",1,0,1);
  vertices_ = fs->make<TH1F>("vertices","",1,0,1);

  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);

  vtxPerf_["Ntrk"] = fs->make<TH1F>("vtxNtrk","Tracks per vertex",50,0,200);
  vtxPerf_["x"] = fs->make<TH1F>("vtxX","Vertex x position",100,-1,1);
  vtxPerf_["y"] = fs->make<TH1F>("vtxY","Vertex y position",100,-1,1);
  vtxPerf_["z"] = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);

  vtxPerf2D_["Ntrk2D"] = fs->make<TH2F>("vtxNtrk2D","Tracks per vertex;vertex (sorted by Ntrk);Ntrk"
                                ,10,0,10,200,0,200);
  
  trkPerf_["Nhit"] = fs->make<TH1F>("trkNhit", "Tracks by Number of Valid Hits;N hits",    30,  0,30);
  trkPerf_["pt"] = fs->make<TH1F>("trkPt", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,10);
  trkPerf_["eta"] = fs->make<TH1F>("trkEta", "Track Pseudorapidity Distribution;#eta",50,-2.5,2.5);
  trkPerf_["phi"] = fs->make<TH1F>("trkPhi", "Track Azimuthal Distribution;#phi",100,-3.15,3.15);
  trkPerf_["chi2"] = fs->make<TH1F>("trkChi2", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",60,0,20);
  trkPerf_["pterr"] = fs->make<TH1F>("trkPterr", "Track p_{T} error;#delta p_{T} / p_{T}",50,0,0.5);
  trkPerf_["dxyErr"] = fs->make<TH1F>("trkDxyErr", "Transverse DCA Significance;dxy / #sigma_{dxy}",100,-20,20);
  trkPerf_["dzErr"] = fs->make<TH1F>("trkDzErr", "Longitudinal DCA Significance;dz / #sigma_{dz}",100,-20,20);

  trkPerf2D_["etaphi"] = fs->make<TH2F>("trkEtaPhi","Track Eta-Phi Map;#eta;#phi",50,-2.5,2.5,100,-3.15,3.15);

}




void
RpPbTrackingAnalyzer::beginJob()
{
}

void
RpPbTrackingAnalyzer::endJob()
{
  // normalize histograms
  for( auto elem : trkPerf_ )
  {
    auto & histo = elem.second;
    histo->Scale(1./(double)histo->Integral());
  }  
  for( auto elem : trkPerf2D_ )
  {
    auto & histo = elem.second;
    histo->Scale(1./(double)histo->Integral());
  }  
  for( auto elem : vtxPerf_ )
  {
    auto & histo = elem.second;
    histo->Scale(1./(double)histo->Integral());
  }  
  for( auto elem : vtxPerf2D_ )
  {
    auto & histo = elem.second;
    histo->Scale(1./(double)histo->Integral());
  }  
  for( auto elem : evtPerf_ )
  {
    auto & histo = elem.second;
    histo->Scale(1./(double)histo->Integral());
  }  



}

DEFINE_FWK_MODULE(RpPbTrackingAnalyzer);
