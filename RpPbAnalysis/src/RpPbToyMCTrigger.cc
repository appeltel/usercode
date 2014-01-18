//
// Original Author:  Eric Appelt
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
#include <TRandom3.h>
#include <TFile.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>

class RpPbToyMCTrigger : public edm::EDAnalyzer {
   public:
      explicit RpPbToyMCTrigger(const edm::ParameterSet&);
      ~RpPbToyMCTrigger();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);
      int findPtBin( const reco::GenParticle & gen );
      double getOfflineEff( double pt, double eta );
      double getHLTEff( double pt, double eta );
      double getJointEff( double pt, double eta );

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> hist1D_;
      std::map<std::string,TH2F*> hist2D_;
      
      std::vector<double> ptBins_;

      double etaMin_;
      double etaMax_;
      double etaShift_;

      double mbPrescale_;
      double trk12Prescale_;
      double trk20Prescale_;

      edm::InputTag genSrc_; 
 
      TRandom3 rand_;

      TFile correctionFile_;
      std::string correctionDir_;

      TH2D * hEffOffline_;
      TH2D * hEffHLT_;
      TH2D * hEffJoint_;

};

RpPbToyMCTrigger::RpPbToyMCTrigger(const edm::ParameterSet& iConfig):
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
etaShift_(iConfig.getParameter<double>("etaShift")),
mbPrescale_(iConfig.getParameter<double>("mbPrescale")),
trk12Prescale_(iConfig.getParameter<double>("trk12Prescale")),
trk20Prescale_(iConfig.getParameter<double>("trk20Prescale")),
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
rand_(),
correctionFile_(iConfig.getParameter<std::string>("correctionFile").c_str()),
correctionDir_(iConfig.getParameter<std::string>("correctionDir"))
{
   edm::Service<TFileService> fs;
   initHistos(fs);
   hEffOffline_ = (TH2D*) correctionFile_.Get(Form("%s/heff",correctionDir_.c_str()));
   hEffHLT_ = (TH2D*) correctionFile_.Get(Form("%s/heffHLT",correctionDir_.c_str()));
   hEffJoint_ = (TH2D*) correctionFile_.Get(Form("%s/heffHLTj",correctionDir_.c_str()));
  
   TH2D * hSim = (TH2D*) correctionFile_.Get(Form("%s/hsim",correctionDir_.c_str()));
   hEffOffline_->Divide(hEffOffline_,hSim,1,1,"B");
   hEffHLT_->Divide(hEffHLT_,hSim,1,1,"B");
   hEffJoint_->Divide(hEffJoint_,hSim,1,1,"B");
   delete hSim;
}

RpPbToyMCTrigger::~RpPbToyMCTrigger()
{
   delete hEffOffline_;
   delete hEffHLT_;
   delete hEffJoint_;
}

void
RpPbToyMCTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::GenParticleCollection> gcol;
   iEvent.getByLabel(genSrc_, gcol);

   // triggers  
   bool isMB = false;
   bool isTrk12 = false;
   bool isTrk20 = false;
   bool isTrk30 = false;
   bool isTrk12HLT = false;
   bool isTrk20HLT = false;
   bool isTrk30HLT = false;

   // make temporary collection of offline  
   // and HLT "reconstructed" particles.
   // and fill gen particle histos
   std::vector<reco::GenParticle> tracks;
   std::vector<reco::GenParticle> tracksHLT;
   double leadGenPt = 0.0;
   for( const auto & gen : * gcol )
   {
      // only consider charged primaries in the 
      // acceptance range of the tracker (|eta|<2.5) 
      if( gen.status() != 1 ||
          gen.charge() == 0 ||
          fabs(gen.eta()+etaShift_) > 2.5 ) continue;

      if( gen.eta()+etaShift_ >= etaMin_ && gen.eta()+etaShift_ <= etaMax_ )      
      {
        hist1D_["genSpectrum"]->Fill(gen.pt());
        if( leadGenPt < gen.pt() ) leadGenPt = gen.pt(); 
      }     
      // reconstruct randomly by efficiency
      // first the offline, and then conditionally the HLT
      double offlineEff = getOfflineEff( gen.pt(), gen.eta()+etaShift_ );
      double jointEff = getJointEff( gen.pt(), gen.eta()+etaShift_ );
      if( rand_.Uniform(0.0,1.0) < offlineEff ) 
      {
         tracks.push_back(gen);
         if( rand_.Uniform(0.0,1.0) < jointEff / offlineEff)
           tracksHLT.push_back(gen);
      } 
      else
      {
         if( rand_.Uniform(0.0,1.0) < jointEff / (1. - offlineEff))
           tracksHLT.push_back(gen);
      }
   }  
   hist1D_["genLeads"]->Fill(leadGenPt);

   // see if the HLT trigger paths fired
   for( const auto & gen : tracksHLT )
   {
     if( gen.pt() >= 12.0 ) isTrk12HLT = true;
     if( gen.pt() >= 20.0 ) isTrk20HLT = true;
     if( gen.pt() >= 30.0 ) isTrk30HLT = true;
   }
   // apply prescales 
   if( rand_.Uniform(0.0,mbPrescale_) < 1.0  ) isMB = true;
   if( rand_.Uniform(0.0,trk12Prescale_) >= 1.0  ) isTrk12HLT = false;
   if( rand_.Uniform(0.0,trk20Prescale_) >= 1.0  ) isTrk20HLT = false;


   // analyze offline "reconstructed" particles for event selection
   double leadTrackPt = 0;
   for( const auto & gen : tracks )
   {
      // only consider offline tracks in selected eta region
      if( gen.eta()+etaShift_ < etaMin_ || gen.eta()+etaShift_ > etaMax_ )      
        continue;

      if( gen.pt() > leadTrackPt ) leadTrackPt = gen.pt();

      if( gen.pt() >= 14.0 && isTrk12HLT ) isTrk12 = true;
      if( gen.pt() >= 22.0 && isTrk20HLT ) isTrk20 = true;
      if( gen.pt() >= 32.0 && isTrk30HLT ) isTrk30 = true;
   }


   // fill histrograms for reconstructed tracks
   // according to fired triggers
   for( const auto & gen : tracks )
   {
      // only consider offline tracks in selected eta region
      if( gen.eta()+etaShift_ < etaMin_ || gen.eta()+etaShift_ > etaMax_ )      
        continue;

      if(isMB) hist1D_["mbSpectrum"]->Fill(gen.pt());   
      if(isTrk12) hist1D_["trk12Spectrum"]->Fill(gen.pt());   
      if(isTrk20) hist1D_["trk20Spectrum"]->Fill(gen.pt());   
      if(isTrk30) hist1D_["trk30Spectrum"]->Fill(gen.pt());   
      if(isMB && leadTrackPt < 14.0 ) 
        hist1D_["mbSpectrumSel"]->Fill(gen.pt());   
      if(isTrk12 && leadTrackPt < 22.0 ) 
        hist1D_["trk12SpectrumSel"]->Fill(gen.pt());   
      if(isTrk20 && leadTrackPt < 32.0 ) 
        hist1D_["trk20SpectrumSel"]->Fill(gen.pt());   
      if(isTrk30) 
        hist1D_["trk30SpectrumSel"]->Fill(gen.pt());   
   }
   if(isMB) hist1D_["mbLeads"]->Fill(leadTrackPt);
   if(isTrk12) hist1D_["trk12Leads"]->Fill(leadTrackPt);
   if(isTrk20) hist1D_["trk20Leads"]->Fill(leadTrackPt);
   if(isTrk30) hist1D_["trk30Leads"]->Fill(leadTrackPt);
}

int 
RpPbToyMCTrigger::findPtBin( const reco::GenParticle & gen )
{
  int ptbin = 0;
  for( unsigned int i=1; i<ptBins_.size(); i++)
  {
    if( gen.pt() < ptBins_[i] )
    {
      ptbin = i-1;
      break;
    }
  }
  return ptbin;
}

void
RpPbToyMCTrigger::initHistos(const edm::Service<TFileService> & fs)
{

  hist1D_["genSpectrum"] = fs->make<TH1F>("genSpectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["genLeads"] = fs->make<TH1F>("genLeads",";p_{T};yield",
                           100,0,200);

  hist1D_["mbSpectrum"] = fs->make<TH1F>("mbSpectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["mbLeads"] = fs->make<TH1F>("mbLeads",";p_{T};yield",
                           100,0,200);
  hist1D_["mbSpectrumSel"] = fs->make<TH1F>("mbSpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk12Spectrum"] = fs->make<TH1F>("trk12Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk12Leads"] = fs->make<TH1F>("trk12Leads",";p_{T};yield",
                           100,0,200);
  hist1D_["trk12SpectrumSel"] = fs->make<TH1F>("trk12SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk20Spectrum"] = fs->make<TH1F>("trk20Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk20Leads"] = fs->make<TH1F>("trk20Leads",";p_{T};yield",
                           100,0,200);
  hist1D_["trk20SpectrumSel"] = fs->make<TH1F>("trk20SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk30Spectrum"] = fs->make<TH1F>("trk30Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk30Leads"] = fs->make<TH1F>("trk30Leads",";p_{T};yield",
                           100,0,200);
  hist1D_["trk30SpectrumSel"] = fs->make<TH1F>("trk30SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

}

double
RpPbToyMCTrigger::getOfflineEff( double pt, double eta )
{
  double eff = hEffOffline_->GetBinContent(
                  hEffOffline_->GetXaxis()->FindBin(eta),
                  hEffOffline_->GetYaxis()->FindBin(pt) );
  if( eff >= 0.9999 || eff <= 0.0001) eff = 1;

  return eff;
}

double
RpPbToyMCTrigger::getHLTEff( double pt, double eta )
{
  double eff = hEffHLT_->GetBinContent(
                  hEffHLT_->GetXaxis()->FindBin(eta),
                  hEffHLT_->GetYaxis()->FindBin(pt) );
  if( eff >= 0.9999 || eff <= 0.0001) eff = 1;

  return eff;
}

double
RpPbToyMCTrigger::getJointEff( double pt, double eta )
{
  double eff = hEffJoint_->GetBinContent(
                  hEffJoint_->GetXaxis()->FindBin(eta),
                  hEffJoint_->GetYaxis()->FindBin(pt) );
  if( eff >= 0.9999 || eff <= 0.0001) eff = 1;

  return eff;
}

void
RpPbToyMCTrigger::beginJob()
{
}

void
RpPbToyMCTrigger::endJob()
{
}


DEFINE_FWK_MODULE(RpPbToyMCTrigger);
