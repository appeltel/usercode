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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

class RpPbTrackingCorrections : public edm::EDAnalyzer {
   public:
      explicit RpPbTrackingCorrections(const edm::ParameterSet&);
      ~RpPbTrackingCorrections();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);
      bool passesCuts(const reco::Track & track, const reco::Vertex & vertex);

      // ----------member data ---------------------------

      std::map<std::string,TH3F*> trkCorr3D_;
      std::map<std::string,TH2F*> trkCorr2D_;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      edm::InputTag tpFakSrc_;
      edm::InputTag tpEffSrc_;
      edm::InputTag associatorMap_;

      std::vector<double> ptBins_;
      std::vector<double> etaBins_;
      std::vector<double> occBins_;

      bool occByCentrality_;
      bool occByNPixelTrk_;
 
      CentralityProvider * centrality_;

      bool applyCuts_;
      std::string qualityString_;
      double dxyErrMax_;
      double dzErrMax_;
      double ptErrMax_;

};

RpPbTrackingCorrections::RpPbTrackingCorrections(const edm::ParameterSet& iConfig):
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
tpEffSrc_(iConfig.getParameter<edm::InputTag>("tpEffSrc")),
tpFakSrc_(iConfig.getParameter<edm::InputTag>("tpFakSrc")),
associatorMap_(iConfig.getParameter<edm::InputTag>("associatorMap")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
occBins_(iConfig.getParameter<std::vector<double> >("occBins")),
occByCentrality_(iConfig.getParameter<bool>("occByCentrality")),
occByNPixelTrk_(iConfig.getParameter<bool>("occByNPixelTrk")),
applyCuts_(iConfig.getParameter<bool>("applyCuts")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
dxyErrMax_(iConfig.getParameter<double>("dzErrMax")),
dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
ptErrMax_(iConfig.getParameter<double>("ptErrMax"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   centrality_ = 0;
}

RpPbTrackingCorrections::~RpPbTrackingCorrections()
{
}

void
RpPbTrackingCorrections::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // collections of simulated particles 
   edm::Handle<TrackingParticleCollection>  TPCollectionHeff, TPCollectionHfake;
   iEvent.getByLabel(tpEffSrc_,TPCollectionHeff);
   iEvent.getByLabel(tpFakSrc_,TPCollectionHfake);

   // association map between tracks and simulated particles
   reco::RecoToSimCollection recSimColl;
   reco::SimToRecoCollection simRecColl;
   edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
   edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
   iEvent.getByLabel(associatorMap_,simtorecoCollectionH);
   simRecColl= *(simtorecoCollectionH.product());
   iEvent.getByLabel(associatorMap_,recotosimCollectionH);
   recSimColl= *(recotosimCollectionH.product());

   // reconstructed tracks
   Handle<edm::View<reco::Track> > tcol;
   iEvent.getByLabel(trackSrc_, tcol);

   // primary vertices
   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);
  
   // sort the vertcies by number of tracks in descending order
   std::vector<reco::Vertex> vsorted = *vertex;
   std::sort( vsorted.begin(), vsorted.end(), 
              []( reco::Vertex a, reco::Vertex b) 
   {
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   });

   // obtain event centrality
   if (!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); 

   // skip events with no PV, this should not happen
   if( vsorted.size() == 0) return;

   // determine occupancy variable for event
   double occ = 0.;  
   if( occByCentrality_) occ = centrality_->centralityValue();   
   if( occByNPixelTrk_) occ = centrality_->raw()->NpixelTracks();   


   // ---------------------
   // loop through reco tracks to fill fake, reco, and secondary histograms
   // ---------------------

   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
    
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());
     // skip tracks that fail cuts, using vertex with most tracks as PV       
     if( !passesCuts(*tr, vsorted[0]) ) continue;

     trkCorr3D_["hrec3D"]->Fill(tr->eta(), tr->pt(), occ);
     trkCorr2D_["hrec"]->Fill(tr->eta(), tr->pt());

     // look for match to simulated particle, use first match if it exists
     std::vector<std::pair<TrackingParticleRef, double> > tp;
     const TrackingParticle *mtp=0;
     if(recSimColl.find(track) != recSimColl.end())
     {
       tp = recSimColl[track];
       mtp = tp.begin()->first.get();  
       if( mtp->status() < 0) trkCorr3D_["hsec3D"]->Fill(tr->eta(), tr->pt(), occ);     
       if( mtp->status() < 0) trkCorr2D_["hsec"]->Fill(tr->eta(), tr->pt());     
     }
     else
     {
       trkCorr3D_["hfak3D"]->Fill(tr->eta(), tr->pt(), occ);
       trkCorr2D_["hfak"]->Fill(tr->eta(), tr->pt());
     }

   }

   // ---------------------
   // loop through sim particles to fill matched, multiple,  and sim histograms 
   // ---------------------

   for(TrackingParticleCollection::size_type i=0; i<TPCollectionHeff->size(); i++) 
   {      
     TrackingParticleRef tpr(TPCollectionHeff, i);
     TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
         
     if(tp->status() < 0 || tp->charge()==0) continue; //only charged primaries
     
     trkCorr3D_["hsim3D"]->Fill(tp->eta(),tp->pt(), occ);
     trkCorr2D_["hsim"]->Fill(tp->eta(),tp->pt());

     // find number of matched reco tracks that pass cuts
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     size_t nrec=0;
     if(simRecColl.find(tpr) != simRecColl.end())
     {
       rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
       for (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit = rt.begin(); rtit != rt.end(); ++rtit)
       {
         const reco::Track* tmtr = rtit->first.get();
         if( ! passesCuts(*tmtr, vsorted[0]) ) continue;
         nrec++;
       }
     }
     
     if(nrec>0) trkCorr3D_["heff3D"]->Fill(tp->eta(),tp->pt(), occ);
     if(nrec>1) trkCorr3D_["hmul3D"]->Fill(tp->eta(),tp->pt(), occ);
     if(nrec>0) trkCorr2D_["heff"]->Fill(tp->eta(),tp->pt());
     if(nrec>1) trkCorr2D_["hmul"]->Fill(tp->eta(),tp->pt());

   }


}

bool
RpPbTrackingCorrections::passesCuts(const reco::Track & track, const reco::Vertex & vertex)
{
   if ( ! applyCuts_ ) return true;

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   vtxPoint=vertex.position();
   vzErr=vertex.zError();
   vxErr=vertex.xError();
   vyErr=vertex.yError();

   double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
   dxy = track.dxy(vtxPoint);
   dz = track.dz(vtxPoint);
   dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
   dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
 
   if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1)
       return false;
   if(fabs(dxy/dxysigma) > dxyErrMax_) return false;
   if(fabs(dz/dzsigma) > dzErrMax_) return false;
   if(track.ptError() / track.pt() > ptErrMax_) return false;

   return true;
}


void
RpPbTrackingCorrections::initHistos(const edm::Service<TFileService> & fs)
{

  std::vector<std::string> hNames3D = { "hsim3D", "hrec3D", "hmul3D", "hfak3D",
                                        "heff3D", "hsec3D" };

  for( auto name : hNames3D )
  {
     trkCorr3D_[name] = fs->make<TH3F>(name.c_str(),";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 
  }

  std::vector<std::string> hNames2D = { "hsim", "hrec", "hmul", "hfak",
                                        "heff", "hsec" };

  for( auto name : hNames2D )
  {
     trkCorr2D_[name] = fs->make<TH2F>(name.c_str(),";#eta;p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0]);
  }


}




void
RpPbTrackingCorrections::beginJob()
{
}

void
RpPbTrackingCorrections::endJob()
{
}

DEFINE_FWK_MODULE(RpPbTrackingCorrections);
