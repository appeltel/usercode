//
// Original Author:  Eric Appelt
//         Created:  Tue Nov 26, 2013
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
#include <TF1.h>

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

class RpPbSimpleTrackingCorrections : public edm::EDAnalyzer {
   public:
      explicit RpPbSimpleTrackingCorrections(const edm::ParameterSet&);
      ~RpPbSimpleTrackingCorrections();

      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);
      bool passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex);

      // ----------member data ---------------------------


      std::map<std::string,TH2F*> trkCorr2D_;
      TH3F * momRes_;
      TH1F * vtxZ_;
      TF1 * vtxWeightFunc_;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      edm::InputTag tpFakSrc_;
      edm::InputTag tpEffSrc_;
      edm::InputTag associatorMap_;

      std::vector<double> ptBins_;
      std::vector<double> etaBins_;

      
      std::vector<double> vtxWeightParameters_;
      bool doVtxReweighting_;

      bool applyVertexZCut_;
      double vertexZMax_;

      bool applyTrackCuts_;
      std::string qualityString_;
      double dxyErrMax_;
      double dzErrMax_;
      double ptErrMax_;

      bool doMomRes_;

};

RpPbSimpleTrackingCorrections::RpPbSimpleTrackingCorrections(const edm::ParameterSet& iConfig):
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
tpFakSrc_(iConfig.getParameter<edm::InputTag>("tpFakSrc")),
tpEffSrc_(iConfig.getParameter<edm::InputTag>("tpEffSrc")),
associatorMap_(iConfig.getParameter<edm::InputTag>("associatorMap")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
vtxWeightParameters_(iConfig.getParameter<std::vector<double> >("vtxWeightParameters")),
doVtxReweighting_(iConfig.getParameter<bool>("doVtxReweighting")),
applyVertexZCut_(iConfig.getParameter<bool>("applyVertexZCut")),
vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
applyTrackCuts_(iConfig.getParameter<bool>("applyTrackCuts")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
dxyErrMax_(iConfig.getParameter<double>("dzErrMax")),
dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
doMomRes_(iConfig.getParameter<bool>("doMomRes"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);

   vtxWeightFunc_ = new TF1("vtxWeight","gaus(0)/gaus(3)",-50.,50.);
   // vtxWeightParameters should have size 6,
   // one really should throw an error if not
   if( (int)vtxWeightParameters_.size() == 6 )
   {
     for( unsigned int i=0;i<vtxWeightParameters_.size(); i++)
       vtxWeightFunc_->SetParameter(i,vtxWeightParameters_[i]);
   }

}

RpPbSimpleTrackingCorrections::~RpPbSimpleTrackingCorrections()
{
   delete vtxWeightFunc_;
}

void
RpPbSimpleTrackingCorrections::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   // obtain collections of simulated particles 
   edm::Handle<TrackingParticleCollection>  TPCollectionHeff, TPCollectionHfake;
   iEvent.getByLabel(tpEffSrc_,TPCollectionHeff);
   iEvent.getByLabel(tpFakSrc_,TPCollectionHfake);

   // obtain association map between tracks and simulated particles
   reco::RecoToSimCollection recSimColl;
   reco::SimToRecoCollection simRecColl;
   edm::Handle<reco::SimToRecoCollection > simtorecoCollectionH;
   edm::Handle<reco::RecoToSimCollection > recotosimCollectionH;
   iEvent.getByLabel(associatorMap_,simtorecoCollectionH);
   simRecColl= *(simtorecoCollectionH.product());
   iEvent.getByLabel(associatorMap_,recotosimCollectionH);
   recSimColl= *(recotosimCollectionH.product());

   // obtain reconstructed tracks
   Handle<edm::View<reco::Track> > tcol;
   iEvent.getByLabel(trackSrc_, tcol);

   // obtain primary vertices
   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);
  
   // sort the vertcies by number of tracks in descending order
   std::vector<reco::Vertex> vsorted = *vertex;
   std::sort( vsorted.begin(), vsorted.end(), RpPbSimpleTrackingCorrections::vtxSort );

   // skip events with no PV, this should not happen
   if( vsorted.size() == 0) return;

   // skip events failing vertex cut
   if( applyVertexZCut_)
   {
     if( fabs(vsorted[0].z()) > vertexZMax_ ) return;
   }

   // determine vertex reweighting factor
   double w = 1.0;
   if ( doVtxReweighting_ )
     w *= vtxWeightFunc_->Eval(vsorted[0].z());

   vtxZ_->Fill(vsorted[0].z(),w);

   // determine event multipliticy 
   int multiplicity =0;
   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());
     if( passesTrackCuts(*tr, vsorted[0]) ) 
       multiplicity++;
   }


   // ---------------------
   // loop through reco tracks to fill fake, reco, and secondary histograms
   // ---------------------

   for(edm::View<reco::Track>::size_type i=0; i<tcol->size(); ++i){
    
     edm::RefToBase<reco::Track> track(tcol, i);
     reco::Track* tr=const_cast<reco::Track*>(track.get());
     // skip tracks that fail cuts, using vertex with most tracks as PV       
     if( ! passesTrackCuts(*tr, vsorted[0]) ) continue;

     trkCorr2D_["hrec"]->Fill(tr->eta(), tr->pt(), w);

     // look for match to simulated particle, use first match if it exists
     std::vector<std::pair<TrackingParticleRef, double> > tp;
     const TrackingParticle *mtp=0;
     if(recSimColl.find(track) != recSimColl.end())
     {
       tp = recSimColl[track];
       mtp = tp.begin()->first.get();  
       if( mtp->status() < 0 ) 
       {
         trkCorr2D_["hsec"]->Fill(tr->eta(), tr->pt(), w);     
       }
     }
     else
     {
       trkCorr2D_["hfak"]->Fill(tr->eta(), tr->pt(), w);
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

     trkCorr2D_["hsim"]->Fill(tp->eta(),tp->pt(), w);

     // find number of matched reco tracks that pass cuts
     std::vector<std::pair<edm::RefToBase<reco::Track>, double> > rt;
     size_t nrec=0;
     if(simRecColl.find(tpr) != simRecColl.end())
     {
       rt = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[tpr];
       std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator rtit;
       for (rtit = rt.begin(); rtit != rt.end(); ++rtit)
       {
         const reco::Track* tmtr = rtit->first.get();
         if( ! passesTrackCuts(*tmtr, vsorted[0]) ) continue;
         nrec++;
         if( doMomRes_ ) momRes_->Fill( tp->eta(), tp->pt(), tmtr->pt(), w);
       }
     }
     if(nrec>0) trkCorr2D_["heff"]->Fill(tp->eta(),tp->pt(), w);
     if(nrec>1) trkCorr2D_["hmul"]->Fill(tp->eta(),tp->pt(), w);
   }
}

bool
RpPbSimpleTrackingCorrections::passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex)
{
   if ( ! applyTrackCuts_ ) return true;

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
RpPbSimpleTrackingCorrections::initHistos(const edm::Service<TFileService> & fs)
{


  std::vector<std::string> hNames2D = { "hsim", "hrec", "hmul", "hfak",
                                        "heff", "hsec" };

  for( auto name : hNames2D )
  {
     trkCorr2D_[name] = fs->make<TH2F>(name.c_str(),";#eta;p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0]);
  }


  vtxZ_ = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);

  std::vector<double> ptBinsFine;
  for( unsigned int bin = 0; bin<ptBins_.size()-1; bin++)
  {
    double bStart = ptBins_[bin];
    double bWid = ptBins_[bin+1] - ptBins_[bin];
    for( int i=0;i<5;i++)
      ptBinsFine.push_back( bStart + (double)i * bWid / 5. );
  }
  ptBinsFine.push_back(ptBins_[ptBins_.size()-1]);

  momRes_ = fs->make<TH3F>("momRes","momentum resolution sim vs reco",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBinsFine.size()-1, &ptBinsFine[0],
                           ptBinsFine.size()-1, &ptBinsFine[0]);
}

bool
RpPbSimpleTrackingCorrections::vtxSort( const reco::Vertex &  a, const reco::Vertex & b )
{
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;  
}

void
RpPbSimpleTrackingCorrections::beginJob()
{
}

void
RpPbSimpleTrackingCorrections::endJob()
{
}

DEFINE_FWK_MODULE(RpPbSimpleTrackingCorrections);
