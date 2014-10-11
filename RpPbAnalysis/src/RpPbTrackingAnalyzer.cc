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
#include <TF1.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "HepPDT/ParticleID.hh"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "HepMC/HeavyIon.h"

class RpPbTrackingAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RpPbTrackingAnalyzer(const edm::ParameterSet&);
      ~RpPbTrackingAnalyzer();
      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);
      bool passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex);

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> evtPerf_;
      std::map<std::string,TH2F*> evtPerf2D_;
      std::map<std::string,TH1F*> trkPerf_;
      std::map<std::string,TH2F*> trkPerf2D_;
      std::map<std::string,TH3F*> trkPerf3D_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;
      
      TH3F* trkSpectrum_; 
      TH3F* trkSpectrumDS_; 
      TH3F* genSpectrum_; 

      TH1F* events_;

      TF1 * vtxWeightFunc_;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      edm::InputTag jetSrc_;
      double vertexZMax_;
 
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
      bool occByMultiplicity_;
      bool occByNcoll_;
      bool occByLeadingJetEt_;
      double jetEtaMax_;

      bool cutByLeadingTrackPt_;
      double leadingTrackPtMin_;
      double leadingTrackPtMax_;
      double leadingTrackEtaMin_;
      double leadingTrackEtaMax_;

      bool doPerformance_;

      bool doMC_;
      bool isHIMC_;
      edm::InputTag genSrc_;
      bool doMCbyTP_;
      edm::InputTag tpSrc_;

      bool doTrigEffCorrection_;
      double zeroMultFraction_;
      std::vector<double> trigEffByMult_;
      std::vector<double> trigContaminationByMult_;

      std::vector<double> vtxWeightParameters_;
      bool doVtxReweighting_;     

};

RpPbTrackingAnalyzer::RpPbTrackingAnalyzer(const edm::ParameterSet& iConfig):
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc")),
vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
applyCuts_(iConfig.getParameter<bool>("applyCuts")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
dxyErrMax_(iConfig.getParameter<double>("dzErrMax")),
dzErrMax_(iConfig.getParameter<double>("dzErrMax")),
ptErrMax_(iConfig.getParameter<double>("ptErrMax")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
occBins_(iConfig.getParameter<std::vector<double> >("occBins")),
occByCentrality_(iConfig.getParameter<bool>("occByCentrality")),
occByNPixelTrk_(iConfig.getParameter<bool>("occByNPixelTrk")),
occByMultiplicity_(iConfig.getParameter<bool>("occByMultiplicity")),
occByNcoll_(iConfig.getParameter<bool>("occByNcoll")),
occByLeadingJetEt_(iConfig.getParameter<bool>("occByLeadingJetEt")),
jetEtaMax_(iConfig.getParameter<double>("jetEtaMax")),
cutByLeadingTrackPt_(iConfig.getParameter<bool>("cutByLeadingTrackPt")),
leadingTrackPtMin_(iConfig.getParameter<double>("leadingTrackPtMin")),
leadingTrackPtMax_(iConfig.getParameter<double>("leadingTrackPtMax")),
leadingTrackEtaMin_(iConfig.getParameter<double>("leadingTrackEtaMin")),
leadingTrackEtaMax_(iConfig.getParameter<double>("leadingTrackEtaMax")),
doPerformance_(iConfig.getParameter<bool>("doPerformance")),
doMC_(iConfig.getParameter<bool>("doMC")),
isHIMC_(iConfig.getParameter<bool>("isHIMC")),
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
doMCbyTP_(iConfig.getParameter<bool>("doMCbyTP")),
tpSrc_(iConfig.getParameter<edm::InputTag>("tpSrc")),
doTrigEffCorrection_(iConfig.getParameter<bool>("doTrigEffCorrection")),
zeroMultFraction_(iConfig.getParameter<double>("zeroMultFraction")),
trigEffByMult_(iConfig.getParameter<std::vector<double> >("trigEffByMult")),
trigContaminationByMult_(iConfig.getParameter<std::vector<double> >("trigContaminationByMult")),
vtxWeightParameters_(iConfig.getParameter<std::vector<double> >("vtxWeightParameters")),
doVtxReweighting_(iConfig.getParameter<bool>("doVtxReweighting"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   centrality_ = 0;
   vtxWeightFunc_ = new TF1("vtxWeight","gaus(0)/gaus(3)",-50.,50.);
   // vtxWeightParameters should have size 6,
   // one really should throw an error if not
   if( (int)vtxWeightParameters_.size() == 6 )
   {
     for( unsigned int i=0;i<vtxWeightParameters_.size(); i++)
       vtxWeightFunc_->SetParameter(i,vtxWeightParameters_[i]);
   }
}

RpPbTrackingAnalyzer::~RpPbTrackingAnalyzer()
{
   delete vtxWeightFunc_;
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
   std::sort( vsorted.begin(), vsorted.end(), RpPbTrackingAnalyzer::vtxSort );
 
   // skip events with no PV, this should not happen
   if( vsorted.size() == 0) return; 
   // skip events failing vertex cut
   if( fabs(vsorted[0].z()) > vertexZMax_ ) return; 

   evtPerf_["evtRaw"]->Fill(0.5);

   // determine event multipliticy and leading track pt (selected tracks)
   int multiplicity =0;
   double leadingTrackPt = 0.0;
   double leadingTrackPtGlobal = 0.0;
   // force track cuts for selection
   bool realApplyCuts = applyCuts_;
   applyCuts_ = true;
   reco::TrackCollection::const_iterator track;
   for(  track = tcol->begin(); track != tcol->end() ; ++track )
   { 
     if( passesTrackCuts(*track, vsorted[0]) ) 
     {
       multiplicity++;
       if( track->pt() > leadingTrackPt ) leadingTrackPtGlobal =  track->pt();
       if( track->pt() > leadingTrackPt && track->eta() < leadingTrackEtaMax_
           && track->eta() > leadingTrackEtaMin_  ) 
         leadingTrackPt =  track->pt();
     }
   }
   applyCuts_ = realApplyCuts;

   // determine weighting factor based on multiplicity
   // if weighting is to be performed
   double w = 1.0;
   if( doTrigEffCorrection_ )
   {
     if( multiplicity > 0 && multiplicity < (int) trigEffByMult_.size())
       w = w / trigEffByMult_[multiplicity];
     if( multiplicity > 0 && multiplicity < (int) trigContaminationByMult_.size())
       w = w * (1.0 - trigContaminationByMult_[multiplicity]);
   }

   // do vertex reweighting of each event if desired
   if ( doVtxReweighting_ )
     w *= vtxWeightFunc_->Eval(vsorted[0].z());


   evtPerf_["eventsLeadTrack"]->Fill( leadingTrackPtGlobal, w); 
   evtPerf_["eventsLeadTrackEta"]->Fill( leadingTrackPt, w); 

   // skip events  by leading track pt
   if (cutByLeadingTrackPt_ )
   {
      if( leadingTrackPt < leadingTrackPtMin_ ) return;
      if( leadingTrackPt > leadingTrackPtMax_ ) return;
   }

   // obtain jets, if we are considering them
   Handle<std::vector<pat::Jet> > jets;
   std::vector<const pat::Jet *> sortedJets;
   if( occByLeadingJetEt_) 
   {
      iEvent.getByLabel(jetSrc_, jets);
      // take only accepted jets and sort them by Et
      for(unsigned it=0; it<jets->size(); ++it)
      {
        const pat::Jet* jet = &((*jets)[it]);
        if(fabs(jet->eta())>=jetEtaMax_) continue; 
        sortedJets.push_back(jet);
      }
      sortByPtRef(&sortedJets); 
   }

   // fill event based histograms
   events_->Fill(0.5,w); 
   evtPerf_["Nvtx"]->Fill(vsorted.size(),w);
   evtPerf_["Ntrk"]->Fill(tcol->size(),w);
   int lumi = iEvent.getLuminosityBlock().luminosityBlock();
   evtPerf_["Lumi"]->Fill(lumi,w);
   evtPerf_["NvtxLumi"]->Fill(lumi,vsorted.size()*w);

   if (!centrality_) centrality_ = new CentralityProvider(iSetup);
   centrality_->newEvent(iEvent,iSetup); 
   evtPerf_["centrality"]->Fill(centrality_->getBin(),w);

   // determine occupancy variable for event
   double occ = 0.;  
   if( occByCentrality_) occ = (double) centrality_->getBin();   
   if( occByNPixelTrk_) occ = centrality_->raw()->NpixelTracks();   
   if( occByMultiplicity_ ) occ = multiplicity;

   double leadJetEt = 0.;
   if ( ! sortedJets.empty() ) leadJetEt = sortedJets[0]->pt();
   if( occByLeadingJetEt_ ) occ = leadJetEt;

   // is this a Double-Sided Event?
   bool isDS = false;

   // find event Ncoll if MC, fill gen Histos
   // also determine if this is a DS event
   if ( doMC_)
   {
     Handle<edm::HepMCProduct> hepmc;
     iEvent.getByLabel("generator",hepmc);
     int ncoll = 2;
     float b = 0;
     if( isHIMC_ )
     {
       ncoll = hepmc->GetEvent()->heavy_ion()->Ncoll();
       b = hepmc->GetEvent()->heavy_ion()->impact_parameter();
     }
     evtPerf_["ncoll"]->Fill( ncoll ,w);
     evtPerf_["b"]->Fill( b,w );
     evtPerf2D_["ncollCent"]->Fill( ncoll, centrality_->getBin(),w);
     evtPerf2D_["bCent"]->Fill( b, centrality_->getBin(),w);
     evtPerf2D_["bncoll"]->Fill( b, ncoll,w);
     if( occByNcoll_) occ = ncoll;

     bool posDS = false; bool negDS = false;
     edm::ESHandle<ParticleDataTable> particleDataTable_;
     iSetup.getData(particleDataTable_);

     Handle<reco::GenParticleCollection> gcol;
     iEvent.getByLabel(genSrc_, gcol);
     for( const auto & gen : * gcol )
     {
       // see if genpartice counts for DS
       HepPDT::ParticleID particleID(gen.pdgId());
       if (particleID.isValid())
       {
         ParticleData const * particleData = particleDataTable_->particle(particleID);
         if (particleData)
         { 
           double ctau =  particleDataTable_->particle(particleID)->lifetime();
           if ( ctau  > 1e-18 || ctau == 0.0 )
           {
             if( gen.energy() > 3.0 && gen.eta() > 3.0 && gen.eta() < 5.0 ) posDS = true;
             if( gen.energy() > 3.0 && gen.eta() < -3.0 && gen.eta() > -5.0 ) negDS = true;
           }
         }
       }
       // fill gen spectrum with genParticles
       if( gen.status() == 1 && gen.charge() != 0 && !doMCbyTP_ )
         genSpectrum_->Fill( gen.eta(), gen.pt(), occ,w );
     }
     if( posDS && negDS ) isDS = true;
   }

   // fill gen spectrum using TP 
   if( doMCbyTP_ )
   {
     edm::Handle<TrackingParticleCollection>  TPCollection;
     iEvent.getByLabel(tpSrc_,TPCollection);
     for(TrackingParticleCollection::size_type i=0; i<TPCollection->size(); i++) 
     {      
       TrackingParticleRef tpr(TPCollection, i);
       TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
       if( tp->status() < 0 || tp->charge()==0) continue; 
       genSpectrum_->Fill( tp->eta(), tp->pt(), occ,w);
     }
   }

   // fill event histos relating to DS 
   evtPerf_["multiplicity"]->Fill(multiplicity,w);
   evtPerf2D_["hevtEM"]->Fill(multiplicity,occ,w);
   if( isDS) evtPerf_["multDS"]->Fill(multiplicity,w);
   if( !isDS) evtPerf_["multNonDS"]->Fill(multiplicity,w);
   if( isDS) evtPerf_["eventsLeadTrackDS"]->Fill( leadingTrackPtGlobal, w);

   // fill vertex performance histograms
   int vcount = 0; 
   for( const auto & vi :  vsorted )
   {
     vtxPerf_["Ntrk"]->Fill(vi.tracksSize(),w);
     vtxPerf_["x"]->Fill(vi.x(),w);
     vtxPerf_["y"]->Fill(vi.y(),w);
     vtxPerf_["z"]->Fill(vi.z(),w);
     vtxPerf2D_["Ntrk2D"]->Fill(vcount,vi.tracksSize(),w);
     vcount++;
   }
   // comparisons between first and additional primary vertices
   for (unsigned int i =1; i<vsorted.size(); i++)
   {
     double dz = fabs( vsorted[i].z() - vsorted[0].z() );
     double dx = fabs( vsorted[i].x() - vsorted[0].x() );
     double dy = fabs( vsorted[i].y() - vsorted[0].y() );
     double dxy  = sqrt ( dx*dx + dy*dy );
     vtxPerf_["assocVtxDz"]->Fill(dz,w);
     vtxPerf2D_["assocVtxDzNtrk"]->Fill(dz,vsorted[i].tracksSize(),w );
     vtxPerf_["assocVtxDxy"]->Fill(dxy);
     vtxPerf2D_["assocVtxDxyNtrk"]->Fill(dxy,vsorted[i].tracksSize(),w );
   }
  
   // use vertex with most tracks as primary vertex
   // determine position and error
   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;
   if( vsorted.size() != 0 )
   {
     vtxPoint=vsorted.begin()->position();
     vzErr=vsorted.begin()->zError();
     vxErr=vsorted.begin()->xError();
     vyErr=vsorted.begin()->yError();
   }
  
   // loop through reconstructed tracks and fill
   // track spectrum and performance histograms
   for( const auto & track : * tcol )
   {
     double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
     dxy = track.dxy(vtxPoint);
     dz = track.dz(vtxPoint);
     dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
     dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);

     if( !passesTrackCuts(track, vsorted[0]) ) continue;

     trkSpectrum_->Fill( track.eta(), track.pt(), occ,w);
     if( isDS  ) trkSpectrumDS_->Fill( track.eta(), track.pt(), occ,w);
     if( doPerformance_)
     {
       trkPerf_["Nhit"]->Fill(track.numberOfValidHits(),w); 
       trkPerf_["pt"]->Fill(track.pt(),w); 
       trkPerf_["eta"]->Fill( track.eta(),w );
       trkPerf2D_["etaphi"]->Fill( track.eta(), track.phi(),w );
       trkPerf_["ptHigh"]->Fill(track.pt(),w); 
       trkPerf_["phi"]->Fill(track.phi(),w); 
       trkPerf_["dxyErr"]->Fill(dxy/dxysigma,w); 
       trkPerf_["dzErr"]->Fill(dz/dzsigma,w); 
       trkPerf_["chi2"]->Fill(track.normalizedChi2(),w);
       trkPerf_["pterr"]->Fill(track.ptError() / track.pt(),w );
       trkPerf2D_["etavz"]->Fill( vtxPoint.z(), track.eta(),w );
       trkPerf3D_["Nhit3D"]->Fill(track.eta(), track.pt(), track.numberOfValidHits(),w);
       trkPerf3D_["phi3D"]->Fill(track.eta(), track.pt(), track.phi(),w);
       trkPerf3D_["dxyErr3D"]->Fill(track.eta(), track.pt(), dxy/dxysigma,w);
       trkPerf3D_["dzErr3D"]->Fill(track.eta(), track.pt(), dz/dzsigma,w);
       trkPerf3D_["chi23D"]->Fill(track.eta(), track.pt(), track.normalizedChi2(),w);
       trkPerf3D_["pterr3D"]->Fill(track.eta(), track.pt(), track.ptError() / track.pt(),w );
     }
   }
}


void
RpPbTrackingAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  events_ = fs->make<TH1F>("events","",1,0,1);

  trkSpectrum_ = fs->make<TH3F>("trkSpectrum",";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 

  trkSpectrumDS_ = fs->make<TH3F>("trkSpectrumDS",";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 

  genSpectrum_ = fs->make<TH3F>("genSpectrum",";#eta;p_{T};occ var",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           occBins_.size()-1, &occBins_[0]); 

  std::vector<double> dumBins;
  dumBins.clear(); 
  for( double i = 0.; i<36.; i += 1.) dumBins.push_back(i);
  trkPerf3D_["Nhit3D"] = fs->make<TH3F>("trkNhit3D", "Tracks by Number of Valid Hits;N hits",    
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();  
  for( double i = -3.15; i<3.151; i += 6.30/50.) dumBins.push_back(i);
  trkPerf3D_["phi3D"] = fs->make<TH3F>("trkPhi3D", "Track Azimuthal Distribution;#phi",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();    
  for( double i = 0.; i<20.01; i += 20./80.) dumBins.push_back(i);
  trkPerf3D_["chi23D"] = fs->make<TH3F>("trkChi23D", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = 0.0; i<0.201; i += 0.2/50.) dumBins.push_back(i);
  trkPerf3D_["pterr3D"] = fs->make<TH3F>("trkPterr3D", "Track p_{T} error;#delta p_{T} / p_{T}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = -20.; i<20.01; i += 40./100.) dumBins.push_back(i);
  trkPerf3D_["dxyErr3D"] = fs->make<TH3F>("trkDxyErr3D", "Transverse DCA Significance;dxy / #sigma_{dxy}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);
  trkPerf3D_["dzErr3D"] = fs->make<TH3F>("trkDzErr3D", "Longitudinal DCA Significance;dz / #sigma_{dz}",
                           etaBins_.size()-1, &etaBins_[0],
                           ptBins_.size()-1, &ptBins_[0],
                           dumBins.size()-1, &dumBins[0]);


  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["evtRaw"] = fs->make<TH1F>("evtRaw","Raw number of events",1,0,1);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);
  evtPerf_["centrality"] = fs->make<TH1F>("centrality","Event centrality bin",100,0,100);
  evtPerf_["ncoll"] = fs->make<TH1F>("ncoll","Event Ncoll from Generator",50,0,50);
  evtPerf_["b"] = fs->make<TH1F>("b","Impact Parameter from Generator",100,0,20);

  evtPerf_["NvtxLumi"] = fs->make<TH1F>("evtNvtxLumi","Primary Vertices by Lumi",200,0,2000);
  evtPerf_["Lumi"] = fs->make<TH1F>("evtLumi","Events by Lumi",200,0,2000);

  evtPerf_["multiplicity"] = fs->make<TH1F>("multiplicity","Event Multiplicity (selected tracks)",500,0,500);
  evtPerf_["multDS"] = fs->make<TH1F>("multDS","DS Event Multiplicity (selected tracks)",500,0,500);
  evtPerf_["multNonDS"] = fs->make<TH1F>("multNonDS","Non-DS Event Multiplicity (selected tracks)",500,0,500);

  evtPerf_["eventsLeadTrack"] = fs->make<TH1F>("eventsLeadTrack", "Number of events by leading track pt",
                                               100,0,200);
  evtPerf_["eventsLeadTrackEta"] = fs->make<TH1F>("eventsLeadTrackEta", "Number of events by leading track pt",
                                               100,0,200);
  evtPerf_["eventsLeadTrackDS"] = fs->make<TH1F>("eventsLeadTrackDS", "Number of events by leading track pt",
                                               100,0,200);

  evtPerf2D_["ncollCent"] = fs->make<TH2F>("ncollCent","Ncoll versus Centrality",
                                50,0,50,100,0,100);
  evtPerf2D_["bCent"] = fs->make<TH2F>("bCent","Impact Parameter versus Centrality",
                                100,0,20,100,0,100);
  evtPerf2D_["bncoll"] = fs->make<TH2F>("bncoll","Impact Parameter versus Ncoll",
                                100,0,20,50,0,50);

  std::vector<double> multBins;
  multBins.clear();
  for( double i = 0.; i<600.; i += 10.) multBins.push_back(i);
     evtPerf2D_["hevtEM"] = fs->make<TH2F>("hevtEM",";mult;occ var",
                           multBins.size()-1, &multBins[0],
                           occBins_.size()-1, &occBins_[0]);

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


bool
RpPbTrackingAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b )
{
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;  
}

void
RpPbTrackingAnalyzer::beginJob()
{
}

void
RpPbTrackingAnalyzer::endJob()
{
   // scale number of events by missing zero multiplicity DS events
   if( doTrigEffCorrection_ )
     events_->Scale(1./(1.-zeroMultFraction_));
}

bool
RpPbTrackingAnalyzer::passesTrackCuts(const reco::Track & track, const reco::Vertex & vertex)
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

DEFINE_FWK_MODULE(RpPbTrackingAnalyzer);
