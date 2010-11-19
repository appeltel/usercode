#include "Appeltel/PixelTrackAnalysis/interface/HiPixelTrkEffHistograms.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"

#include <iostream>
#include <cmath>
using namespace std;

HiPixelTrkEffHistograms::HiPixelTrkEffHistograms(const edm::ParameterSet& pset)
{
  fillHistograms         = pset.getParameter<bool>("fillHistograms");
  fillNtuples            = pset.getParameter<bool>("fillNtuples");
  constPtBins            = pset.getParameter<bool>("constPtBins");
  lowPtMode              = pset.getParameter<bool>("lowPtMode");
}


HiPixelTrkEffHistograms::~HiPixelTrkEffHistograms()
{
}

void 
HiPixelTrkEffHistograms::declareHistograms()
{
  if(fillNtuples) {
    
    TString leafStr;
    
    trackTrees.push_back(f->make<TTree>("simTrackTree","simTrackTree"));
    leafStr = "ids/I:etas/F:pts/F:hits/I:status/I:acc/I:nrec/I:ptr/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I";
    trackTrees[0]->Branch("simTrackValues", &simTrackValues, leafStr.Data());
    
    trackTrees.push_back(f->make<TTree>("recTrackTree","recTrackTree"));
    leafStr = "charge/I:etar/F:ptr/F:phir/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I:nsim/I:status/I:ids/I:parids/I:etas/F:pts/F";
    trackTrees[1]->Branch("recTrackValues", &recTrackValues, leafStr.Data());
    
  }
  
  if(fillHistograms) {

    // pt bins
    if(!constPtBins){
       const double small = 1e-3;
       double pt;

       // simple rebinning possible with a rebinning factor n = 2, 3, 4 !
       for(pt =   0  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 24 bins
       for(pt =   1.2; pt <   2.4-small; pt +=  0.1 ) ptBins.push_back(pt); // 12 bins
       for(pt =   2.4; pt <   7.2-small; pt +=  0.2 ) ptBins.push_back(pt); // 24 bins
       for(pt =   7.2; pt <  13.2-small; pt +=  0.5 ) ptBins.push_back(pt); // 12 bins
       for(pt =  13.2; pt <  25.2-small; pt +=  1.0 ) ptBins.push_back(pt); // 12 bins
       for(pt =  25.2; pt <  61.2-small; pt +=  3.0 ) ptBins.push_back(pt); // 12 bins
       for(pt =  61.2; pt < 121.2-small; pt +=  5.0 ) ptBins.push_back(pt); // 12 bins
       ptBins.push_back(121.2);

    }else if(lowPtMode){

       static float ptMin   =  0.0;
       static float ptMax   =  2.0;
       static float ptWidth =  0.02;

       for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }else{

       static float ptMin   =  0.0;
       static float ptMax   =  200.0;
       static float ptWidth =  0.2;

       for(double pt = ptMin; pt < ptMax + ptWidth/2; pt += ptWidth)
	  ptBins.push_back(pt);

    }
    
    // eta bins
    static float etaMin   = -3.0;
    static float etaMax   =  3.0;
    static float etaWidth =  0.2;

    for(double eta = etaMin; eta < etaMax + etaWidth/2; eta += etaWidth)
      etaBins.push_back(eta);


    // jet et (centrality) bins
    static float jetMin = 0.0;
    static float jetMax = 40.0;  
    static float jetWidth = 2.0;  // 5% centrality bins

    for(double jet = jetMin; jet < jetMax + jetWidth/2; jet += jetWidth)
       jetBins.push_back(jet);


    // simulated
    hsim = f->make<TH2F>("hsim","Sim Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // accepted
    hacc = f->make<TH2F>("hacc","Accepted Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // efficiency
    heff = f->make<TH2F>("heff","Effic Rec Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // multiply reconstructed
    hmul = f->make<TH2F>("hmul","Mult Rec Tracks;#eta;p_{T} (GeV/c)",
			 etaBins.size()-1, &etaBins[0],
			 ptBins.size()-1, &ptBins[0]);

    // reconstructed
    hrec = f->make<TH2F>("hrec","Rec Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // fakes
    hfak = f->make<TH2F>("hfak","Fake Tracks;#eta;p_{T} (GeV/c)",
		    etaBins.size()-1, &etaBins[0],
		    ptBins.size()-1, &ptBins[0]);

    // secondary
    hsec = f->make<TH2F>("hsec","Secondary Tracks;#eta;p_{T} (GeV/c)",
			 etaBins.size()-1, &etaBins[0],
			 ptBins.size()-1, &ptBins[0]);
    

    // simulated 3D 
    hsim3D = f->make<TH3F>("hsim3D","Sim Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // efficiency  3D 
    heff3D = f->make<TH3F>("heff3D","Effic Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // multiply reconstructed 3D 
    hmul3D = f->make<TH3F>("hmul3D","Mult Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                           etaBins.size()-1, &etaBins[0],
                           ptBins.size()-1, &ptBins[0],
                           jetBins.size()-1, &jetBins[0]);


    // reconstructed 3D 
    hrec3D = f->make<TH3F>("hrec3D","Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // fakes 3D 
    hfak3D = f->make<TH3F>("hfak3D","Fake Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // secondary
    hsec3D = f->make<TH3F>("hsec3D","Secondary Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
			   etaBins.size()-1, &etaBins[0],
			   ptBins.size()-1, &ptBins[0],
			   jetBins.size()-1, &jetBins[0]);

    // mom resolution (Sim to Rec) 
    hresStoR3D = f->make<TH3F>("hresStoR3D","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               ptBins.size()-1, &ptBins[0]);

    /*
    hresStoR3D_etaS = f->make<TH3F>("hresStoR3D_etaS","Momentum resolution (sim to rec);#jet E_{T} (GeV/c);sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
				    jetBins.size()-1, &jetBins[0],
				    ptBins.size()-1, &ptBins[0],
				    ptBins.size()-1, &ptBins[0]);

    hresStoR3D_etaL = f->make<TH3F>("hresStoR3D_etaL","Momentum resolution (sim to rec);#jet E_{T} (GeV/c);sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
                                    jetBins.size()-1, &jetBins[0],
                                    ptBins.size()-1, &ptBins[0],
                                    ptBins.size()-1, &ptBins[0]);
    */

    // mom resolution (Sim to Rec) v2      
    //hresStoR3D_v2 = f->make<TH3F>("hresStoR3D_v2","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c);rec p_{T} (GeV/c)",
    //etaBins.size()-1, &etaBins[0],
    //ptBins.size()-1, &ptBins[0],
    //ptBins.size()-1, &ptBins[0]);


  }

}

void 
HiPixelTrkEffHistograms::fillSimHistograms(const SimTrack_t & s)
{

  if(fillNtuples && s.status>0){
    simTrackValues = s;
    trackTrees[0]->Fill();
  }

  if(fillHistograms && s.status>0) {
    hsim->Fill(s.etas, s.pts);
    hsim3D->Fill(s.etas, s.pts, s.jetr);
    if(s.acc)    hacc->Fill(s.etas, s.pts);
    if(s.nrec==1) {
       hresStoR3D->Fill(s.etas, s.pts, s.ptr);
       //if(fabs(s.etas)<1.0) hresStoR3D_etaS->Fill(s.jetr, s.pts, s.ptr);
       //if(fabs(s.etas)<2.4) hresStoR3D_etaL->Fill(s.jetr, s.pts, s.ptr);
    }
    if(s.nrec>0) heff->Fill(s.etas, s.pts), heff3D->Fill(s.etas, s.pts, s.jetr);
    if(s.nrec>1) hmul->Fill(s.etas, s.pts), hmul3D->Fill(s.etas, s.pts, s.jetr);
  }

}

void 
HiPixelTrkEffHistograms::fillRecHistograms(const RecTrack_t & r)
{

  if(fillNtuples){
    recTrackValues = r;
    trackTrees[1]->Fill();
  }

  if(fillHistograms) {
    hrec->Fill(r.etar, r.ptr);
    hrec3D->Fill(r.etar, r.ptr, r.jetr);
    if(!r.nsim) hfak->Fill(r.etar, r.ptr), hfak3D->Fill(r.etar, r.ptr, r.jetr);
    if(r.nsim>0 && r.status<0) hsec->Fill(r.etar, r.ptr), hsec3D->Fill(r.etar, r.ptr, r.jetr); // nsim>0 redudant?
  }

}

void 
HiPixelTrkEffHistograms::writeHistograms()
{

}



