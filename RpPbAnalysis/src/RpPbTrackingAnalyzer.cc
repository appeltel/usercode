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
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>

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
      std::map<std::string,TH3F*> trkPerf3D_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;
      
      TH3F* trkSpectrum_; 

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
      double ptMin_;
 
      CentralityProvider * centrality_;

      bool applyCuts_;
      std::string qualityString_;
      double dxyErrMax_;
      double dzErrMax_;
      double ptErrMax_;

      std::vector<double> ptBins_;
      std::vector<double> etaBins_;
      std::vector<double> occBins_;

      bool occByCentrality_;
      bool occByNPixelTrk_;

};

RpPbTrackingAnalyzer::RpPbTrackingAnalyzer(const edm::ParameterSet& iConfig):
nevt_(0),
ntrack_(0),
nvertex_(0),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
ptMin_(iConfig.getParameter<double>("ptMin")),
applyCuts_(iConfig.getParameter<bool>("applyCuts")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
dxyErrMax_(iConfig.getParameter<double>("dzErrMax")),
dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
occBins_(iConfig.getParameter<std::vector<double> >("occBins")),
occByCentrality_(iConfig.getParameter<bool>("occByCentrality")),
occByNPixelTrk_(iConfig.getParameter<bool>("occByNPixelTrk"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   centrality_ = 0;
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
   // use chi2 as tiebreaker
   std::sort( vsorted.begin(), vsorted.end(), 
              []( reco::Vertex a, reco::Vertex b) 
   {
      if( a.tracksSize() != b.tracksSize() )
        return  a.tracksSize() > b.tracksSize() ? true : false ;
      else
        return  a.chi2() < b.chi2() ? true : false ;  
   });

  
   nevt_++;
   events_->Fill(0.5); 
   evtPerf_["Nvtx"]->Fill(vsorted.size());
   evtPerf_["Ntrk"]->Fill(tcol->size());

   int lumi = iEvent.getLuminosityBlock().luminosityBlock();
   evtPerf_["Lumi"]->Fill(lumi);
   evtPerf_["NvtxLumi"]->Fill(lumi,vsorted.size());

   if (!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); 
   evtPerf_["centrality"]->Fill(centrality_->getBin());

   // determine occupancy variable for event
   double occ = 0.;  
   if( occByCentrality_) occ = centrality_->centralityValue();   
   if( occByNPixelTrk_) occ = centrality_->raw()->NpixelTracks();   

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

   for (unsigned int i =1; i<vsorted.size(); i++)
   {
     double dz = fabs( vsorted[i].z() - vsorted[0].z() );
     double dx = fabs( vsorted[i].x() - vsorted[0].x() );
     double dy = fabs( vsorted[i].y() - vsorted[0].y() );
     double dxy  = sqrt ( dx*dx + dy*dy );
     vtxPerf_["assocVtxDz"]->Fill(dz);
     vtxPerf2D_["assocVtxDzNtrk"]->Fill(dz,vsorted[i].tracksSize() );
     vtxPerf_["assocVtxDxy"]->Fill(dxy);
     vtxPerf2D_["assocVtxDxyNtrk"]->Fill(dxy,vsorted[i].tracksSize() );
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

     tracks_->Fill(0.5);

     double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
     dxy = track.dxy(vtxPoint);
     dz = track.dz(vtxPoint);
     dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
     dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
       
     if ( applyCuts_)
     {
       if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1) 
        continue;
       if(fabs(dxy/dxysigma) > dxyErrMax_) continue;
       if(fabs(dz/dzsigma) > dzErrMax_) continue;
       if(track.ptError() / track.pt() > ptErrMax_) continue;
     }

     trkSpectrum_->Fill( track.eta(), track.pt(), occ);
         
     if( track.eta() <= etaMax_ && track.eta() >= etaMin_ && track.pt() > ptMin_)
     {

       trackseta_->Fill(0.5);
       trkPerf_["Nhit"]->Fill(track.numberOfValidHits()); 
       trkPerf_["pt"]->Fill(track.pt()); 
       trkPerf_["eta"]->Fill( track.eta() );
       trkPerf2D_["etaphi"]->Fill( track.eta(), track.phi() );
       trkPerf_["ptHigh"]->Fill(track.pt()); 
       trkPerf_["phi"]->Fill(track.phi()); 
       trkPerf_["dxyErr"]->Fill(dxy/dxysigma); 
       trkPerf_["dzErr"]->Fill(dz/dzsigma); 
       trkPerf_["chi2"]->Fill(track.normalizedChi2());
       trkPerf_["pterr"]->Fill(track.ptError() / track.pt() );
       trkPerf2D_["etavz"]->Fill( vtxPoint.z(), track.eta() );
       trkPerf3D_["Nhit3D"]->Fill(track.eta(), track.pt(), track.numberOfValidHits());
       trkPerf3D_["phi3D"]->Fill(track.eta(), track.pt(), track.phi());
       trkPerf3D_["dxyErr3D"]->Fill(track.eta(), track.pt(), dxy/dxysigma);
       trkPerf3D_["dzErr3D"]->Fill(track.eta(), track.pt(), dz/dzsigma);
       trkPerf3D_["chi23D"]->Fill(track.eta(), track.pt(), track.normalizedChi2());
       trkPerf3D_["pterr3D"]->Fill(track.eta(), track.pt(), track.ptError() / track.pt() );
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

  trkSpectrum_ = fs->make<TH3F>("trkSpectrum",";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 

  std::vector<double> dumBins;
  dumBins.clear(); 
  for( double i = 0.; i<36.; i += 1.) dumBins.push_back(i);
  trkPerf3D_["Nhit3D"] = fs->make<TH3F>("trkNhit", "Tracks by Number of Valid Hits;N hits",    
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();  
  for( double i = -3.15; i<3.151; i += 6.30/100.) dumBins.push_back(i);
  trkPerf3D_["phi3D"] = fs->make<TH3F>("trkPhi", "Track Azimuthal Distribution;#phi",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();    
  for( double i = 0.; i<6.01; i += 6./60.) dumBins.push_back(i);
  trkPerf3D_["chi23D"] = fs->make<TH3F>("trkChi2", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = 0.0; i<0.201; i += 0.2/50.) dumBins.push_back(i);
  trkPerf3D_["pterr3D"] = fs->make<TH3F>("trkPterr", "Track p_{T} error;#delta p_{T} / p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = -8.; i<8.01; i += 16./100.) dumBins.push_back(i);
  trkPerf3D_["dxyErr3D"] = fs->make<TH3F>("trkDxyErr", "Transverse DCA Significance;dxy / #sigma_{dxy}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  trkPerf3D_["dzErr3D"] = fs->make<TH3F>("trkDzErr", "Longitudinal DCA Significance;dz / #sigma_{dz}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);


  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);
  evtPerf_["centrality"] = fs->make<TH1F>("centrality","Event centrality bin",100,0,100);

  evtPerf_["NvtxLumi"] = fs->make<TH1F>("evtNvtxLumi","Primary Vertices by Lumi",200,0,2000);
  evtPerf_["Lumi"] = fs->make<TH1F>("evtLumi","Events by Lumi",200,0,2000);

  vtxPerf_["Ntrk"] = fs->make<TH1F>("vtxNtrk","Tracks per vertex",50,0,200);
  vtxPerf_["x"] = fs->make<TH1F>("vtxX","Vertex x position",1000,-1,1);
  vtxPerf_["y"] = fs->make<TH1F>("vtxY","Vertex y position",1000,-1,1);
  vtxPerf_["z"] = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);
  vtxPerf_["assocVtxDz"] = fs->make<TH1F>("assocVtxDz","Z Distance from first PV; dz (cm)",200,0,50);
  vtxPerf_["assocVtxDxy"] = fs->make<TH1F>("assocVtxDxy","Rho Distance from first PV; dxy (cm)",200,0,4);

  vtxPerf2D_["Ntrk2D"] = fs->make<TH2F>("vtxNtrk2D","Tracks per vertex;vertex (sorted by Ntrk);Ntrk"
                                ,10,0,10,200,0,200);
  vtxPerf2D_["assocVtxDzNtrk"] = fs->make<TH2F>("assocVtxDzNtrk",
                                 "Z Distance from first PV vs Ntrk of assoc; dz (cm); Ntrk",
                                 200,0,50,50,0,200);
  vtxPerf2D_["assocVtxDxyNtrk"] = fs->make<TH2F>("assocVtxDxyNtrk",
                                 "Rho Distance from first PV vs Ntrk of assoc; dxy (cm); Ntrk",
                                 200,0,4,50,0,200);
  
  trkPerf_["Nhit"] = fs->make<TH1F>("trkNhit", "Tracks by Number of Valid Hits;N hits",    35,  0,35);
  trkPerf_["pt"] = fs->make<TH1F>("trkPt", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,6);
  trkPerf_["ptHigh"] = fs->make<TH1F>("trkPtHigh", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,200);
  trkPerf_["eta"] = fs->make<TH1F>("trkEta", "Track Pseudorapidity Distribution;#eta",50,-2.5,2.5);
  trkPerf_["phi"] = fs->make<TH1F>("trkPhi", "Track Azimuthal Distribution;#phi",100,-3.15,3.15);
  trkPerf_["chi2"] = fs->make<TH1F>("trkChi2", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",60,0,6);
  trkPerf_["pterr"] = fs->make<TH1F>("trkPterr", "Track p_{T} error;#delta p_{T} / p_{T}",50,0,0.2);
  trkPerf_["dxyErr"] = fs->make<TH1F>("trkDxyErr", "Transverse DCA Significance;dxy / #sigma_{dxy}",100,-8,8);
  trkPerf_["dzErr"] = fs->make<TH1F>("trkDzErr", "Longitudinal DCA Significance;dz / #sigma_{dz}",100,-8,8);

  trkPerf2D_["etaphi"] = fs->make<TH2F>("trkEtaPhi","Track Eta-Phi Map;#eta;#phi",50,-2.5,2.5,100,-3.15,3.15);
  trkPerf2D_["etavz"] = fs->make<TH2F>("trkEtaVz","Track Eta vs Vertex z;Vertex z (cm);#eta",
                                       100,-30,30,100,-3.0,3.0);

}




void
RpPbTrackingAnalyzer::beginJob()
{
}

void
RpPbTrackingAnalyzer::endJob()
{
}

DEFINE_FWK_MODULE(RpPbTrackingAnalyzer);
