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

      // ----------member data ---------------------------

      std::map<std::string,TH1F*> hist1D_;
      std::map<std::string,TH2F*> hist2D_;
      
      std::vector<double> ptBins_;
      std::vector<double> hltTrk12Eff_;
      std::vector<double> hltTrk20Eff_;
      std::vector<double> hltTrk30Eff_;
      std::vector<double> offlineEff_;

      double mbPrescale_;
      double trk12Prescale_;
      double trk20Prescale_;

      edm::InputTag genSrc_; 

 
      TRandom3 rand_;

};

RpPbToyMCTrigger::RpPbToyMCTrigger(const edm::ParameterSet& iConfig):
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
hltTrk12Eff_(iConfig.getParameter<std::vector<double> >("hltTrk12Eff")),
hltTrk20Eff_(iConfig.getParameter<std::vector<double> >("hltTrk20Eff")),
hltTrk30Eff_(iConfig.getParameter<std::vector<double> >("hltTrk30Eff")),
offlineEff_(iConfig.getParameter<std::vector<double> >("offlineEff")),
mbPrescale_(iConfig.getParameter<double>("mbPrescale")),
trk12Prescale_(iConfig.getParameter<double>("trk12Prescale")),
trk20Prescale_(iConfig.getParameter<double>("trk20Prescale")),
genSrc_(iConfig.getParameter<edm::InputTag>("genSrc")),
rand_()
{
   edm::Service<TFileService> fs;
   initHistos(fs);
}

RpPbToyMCTrigger::~RpPbToyMCTrigger()
{
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

   // make a temporary collection of offline "reconstructed" particles.
   // and fill gen particle histos
   std::vector<reco::GenParticle> tracks;
   double leadGenPt = 0.0;
   for( const auto & gen : * gcol )
   {
      // only consider charged primaries in the 
      // acceptance range of interest (|eta| < 1.0) 
      if( gen.status() != 1 ||
          gen.charge() == 0 ||
          fabs(gen.eta()) > 2.4 ) continue;

      hist1D_["genSpectrum"]->Fill(gen.pt());
      if( leadGenPt < gen.pt() ) leadGenPt = gen.pt(); 
     
      // find pt bin of the particle
      int ptbin = findPtBin(gen);

      // reconstruct randomly by efficiency
      if( rand_.Uniform(0.0,1.0) < offlineEff_[ptbin] ) tracks.push_back(gen);
   }  
   hist1D_["genLeads"]->Fill(leadGenPt);

   // analyze offline "reconstructed" particles for triggering
   double leadTrackPt = 0;
   for( const auto & gen : tracks )
   {
      if( gen.pt() > leadTrackPt ) leadTrackPt = gen.pt();

      // find pt bin of the particle
      int ptbin = findPtBin(gen);

      // don't do this, all events assumed to fire MB trigger 
      // check for over 400 MeV/c track (MB trigger criterion)
      //if( gen.pt() > 0.4 ) has400MeVTrk = true;
 
      // fire triggers based on effiency
      if( rand_.Uniform(0.0,1.0) < hltTrk12Eff_[ptbin] ) isTrk12 = true;
      if( rand_.Uniform(0.0,1.0) < hltTrk20Eff_[ptbin] ) isTrk20 = true;
      if( rand_.Uniform(0.0,1.0) < hltTrk30Eff_[ptbin] ) isTrk30 = true;
   }

   // apply prescales 
   if( rand_.Uniform(0.0,mbPrescale_) < 1.0  ) isMB = true;
   if( rand_.Uniform(0.0,trk12Prescale_) >= 1.0  ) isTrk12 = false;
   if( rand_.Uniform(0.0,trk20Prescale_) >= 1.0  ) isTrk20 = false;


   // fill histrograms for reconstructed tracks
   // according to fired triggers
   for( const auto & gen : tracks )
   {
      if(isMB) hist1D_["mbSpectrum"]->Fill(gen.pt());   
      if(isTrk12) hist1D_["trk12Spectrum"]->Fill(gen.pt());   
      if(isTrk20) hist1D_["trk20Spectrum"]->Fill(gen.pt());   
      if(isTrk30) hist1D_["trk30Spectrum"]->Fill(gen.pt());   
      if(isMB && leadTrackPt < 12.0) 
        hist1D_["mbSpectrumSel"]->Fill(gen.pt());   
      if(isTrk12 && leadTrackPt < 19.2 && leadTrackPt >= 12.0 ) 
        hist1D_["trk12SpectrumSel"]->Fill(gen.pt());   
      if(isTrk20 && leadTrackPt < 28.8 && leadTrackPt >= 19.2 ) 
        hist1D_["trk20SpectrumSel"]->Fill(gen.pt());   
      if(isTrk30 && leadTrackPt >= 28.8 ) 
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
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["mbSpectrum"] = fs->make<TH1F>("mbSpectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["mbLeads"] = fs->make<TH1F>("mbLeads",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["mbSpectrumSel"] = fs->make<TH1F>("mbSpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk12Spectrum"] = fs->make<TH1F>("trk12Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk12Leads"] = fs->make<TH1F>("trk12Leads",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk12SpectrumSel"] = fs->make<TH1F>("trk12SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk20Spectrum"] = fs->make<TH1F>("trk20Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk20Leads"] = fs->make<TH1F>("trk20Leads",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk20SpectrumSel"] = fs->make<TH1F>("trk20SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

  hist1D_["trk30Spectrum"] = fs->make<TH1F>("trk30Spectrum",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk30Leads"] = fs->make<TH1F>("trk30Leads",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);
  hist1D_["trk30SpectrumSel"] = fs->make<TH1F>("trk30SpectrumSel",";p_{T};yield",
                           ptBins_.size()-1, &ptBins_[0]);

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
