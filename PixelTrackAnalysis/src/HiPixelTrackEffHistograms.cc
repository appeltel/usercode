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
  fillHistograms  = pset.getParameter<bool>("fillHistograms");
  fillNtuples     = pset.getParameter<bool>("fillNtuples");
  constPtBins     = pset.getParameter<bool>("constPtBins");
  lowPtMode       = pset.getParameter<bool>("lowPtMode");
  flow2010Mode    = pset.exists("flow2010Mode") ? pset.getParameter<bool>("flow2010Mode") : false;
  jet2010Mode    = pset.exists("jet2010Mode") ? pset.getParameter<bool>("jet2010Mode") : false;
  flow2010ShengquanMode    = pset.exists("flow2010ShengquanMode") ? pset.getParameter<bool>("flow2010ShengquanMode") : false;
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
    leafStr = "ids/I:etas/F:pts/F:phis/F:hits/I:status/I:acc/I:nrec/I:ptr/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I";
    trackTrees[0]->Branch("simTrackValues", &simTrackValues, leafStr.Data());
    
    trackTrees.push_back(f->make<TTree>("recTrackTree","recTrackTree"));
    leafStr = "charge/I:etar/F:ptr/F:phir/F:dz/F:d0/F:pterr/F:d0err/F:dzerr/F:hitr/I:algo/I:nsim/I:status/I:ids/I:parids/I:etas/F:pts/F";
    trackTrees[1]->Branch("recTrackValues", &recTrackValues, leafStr.Data());
    
  }
  

  if(fillHistograms) {

    for( double res = -0.2; res < 0.2000001; res += 0.002) resBins.push_back(res);


    for( double cosb = -1.0; cosb < 1.00001; cosb += 0.02) cosBins.push_back(cosb);

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

    if (flow2010Mode)
    {

       ptBins.clear();
       etaBins.clear();

      ptBins.push_back(0.2);
      ptBins.push_back(0.3);
      ptBins.push_back(0.4);
      ptBins.push_back(0.5);
      ptBins.push_back(0.6);
      ptBins.push_back(0.8);
      ptBins.push_back(1.0);
      ptBins.push_back(1.2);
      ptBins.push_back(1.6);
      ptBins.push_back(2.0);
      ptBins.push_back(2.5);
      ptBins.push_back(3.0);
      ptBins.push_back(4.0);
      ptBins.push_back(6.0);
      ptBins.push_back(8.0);
      ptBins.push_back(12.0);
      for(double eta = -2.4; eta < 2.5; eta += 0.2)
      etaBins.push_back(eta);
    }

    if (flow2010ShengquanMode)
    {


       const double small = 1e-3;
       double pt;
       ptBins.clear();
       etaBins.clear();


       for(pt =   0  ; pt <   1.2-small; pt +=  0.05) ptBins.push_back(pt); // 24 bins
       for(pt =   1.2; pt <   2.0-small; pt +=  0.1 ) ptBins.push_back(pt); // 12 bins
      ptBins.push_back(2.0);
      ptBins.push_back(2.5);
      ptBins.push_back(3.0);
      ptBins.push_back(4.0);
      ptBins.push_back(6.0);
      ptBins.push_back(8.0);
      ptBins.push_back(12.0);

       for(double eta = -2.4; eta < 2.5; eta += 0.2)
       etaBins.push_back(eta);
    }

if (jet2010Mode)
{

   ptBins.clear();
   etaBins.clear();

  ptBins.push_back(0.2);
  ptBins.push_back(0.3);
  ptBins.push_back(0.4);
  ptBins.push_back(0.5);
  ptBins.push_back(0.6);
  ptBins.push_back(0.7);
  ptBins.push_back(0.8);
  ptBins.push_back(0.9);
  ptBins.push_back(1.0);
  ptBins.push_back(1.1);
  ptBins.push_back(1.2);
  ptBins.push_back(1.3);
  ptBins.push_back(1.4);
  ptBins.push_back(1.5);
  ptBins.push_back(1.6);
  ptBins.push_back(1.7);
  ptBins.push_back(1.8);
  ptBins.push_back(1.9);
  ptBins.push_back(2.0);
  ptBins.push_back(2.5);
  ptBins.push_back(3.0);
  ptBins.push_back(4.0);
  ptBins.push_back(5.0);
  ptBins.push_back(6.0);
  ptBins.push_back(8.0);
  ptBins.push_back(12.0);
  ptBins.push_back(16.0);
  ptBins.push_back(20.0);
  ptBins.push_back(30.0);
  ptBins.push_back(50.0);
  ptBins.push_back(80.0);
  ptBins.push_back(120.0);
  etaBins.push_back(0.0);
  etaBins.push_back(0.4);
  etaBins.push_back(0.8);
  etaBins.push_back(1.2);
  etaBins.push_back(1.6);
  etaBins.push_back(2.0);
  etaBins.push_back(2.4);
}



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
    
    // accepted 3D 
    hacc3D = f->make<TH3F>("hacc3D","Accepted Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


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

    hreff3D = f->make<TH3F>("hreff3D","Effic Rec Tracks;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
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

    // simulated 3D v2
    hsim3D_v2 = f->make<TH3F>("hsim3D_v2","Sim Tracks v2;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // reconstructed  3D v2
    hrec3D_v2 = f->make<TH3F>("hrec3D_v2","Rec Tracks v2;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


    // fakes 3D v2
    hfak3D_v2 = f->make<TH3F>("hfak3D_v2","Fake Tracks v2;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // matched 3D v2
    heff3D_v2 = f->make<TH3F>("heff3D_v2","Matched Tracks v2;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    hreff3D_v2 = f->make<TH3F>("hreff3D_v2","Matched Tracks v2;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // simulated 3D v42
    hsim3D_v42 = f->make<TH3F>("hsim3D_v42","Sim Tracks v42;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // reconstructed  3D v42
    hrec3D_v42 = f->make<TH3F>("hrec3D_v42","Rec Tracks v42;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


    // fakes 3D v42
    hfak3D_v42 = f->make<TH3F>("hfak3D_v42","Fake Tracks v42;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // matched 3D v42
    heff3D_v42 = f->make<TH3F>("heff3D_v42","Matched Tracks v42;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    hreff3D_v42 = f->make<TH3F>("hreff3D_v42","Matched Tracks v42;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // simulated 3D v62
    hsim3D_v62 = f->make<TH3F>("hsim3D_v62","Sim Tracks v62;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // reconstructed  3D v62
    hrec3D_v62 = f->make<TH3F>("hrec3D_v62","Rec Tracks v62;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


    // fakes 3D v62
    hfak3D_v62 = f->make<TH3F>("hfak3D_v62","Fake Tracks v62;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // matched 3D v62
    heff3D_v62 = f->make<TH3F>("heff3D_v62","Matched Tracks v62;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    hreff3D_v62 = f->make<TH3F>("hreff3D_v62","Matched Tracks v62;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

     // RXN Plane vn

    for( int i = 0; i<5; i++)
    {
      hsim3D_vn[i] = f->make<TH3F>(Form("hsim3D_v%dRP",i+2),
                        Form("Sim Tracks v%d RP;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",i+2),
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0],
                        jetBins.size()-1, &jetBins[0]);

      hrec3D_vn[i] = f->make<TH3F>(Form("hrec3D_v%dRP",i+2),
                        Form("Reco Tracks v%d RP;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",i+2),
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0],
                        jetBins.size()-1, &jetBins[0]);

      hfak3D_vn[i] = f->make<TH3F>(Form("hfak3D_v%dRP",i+2),
                        Form("Fake Tracks v%d RP;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",i+2),
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0],
                        jetBins.size()-1, &jetBins[0]);

      heff3D_vn[i] = f->make<TH3F>(Form("heff3D_v%dRP",i+2),
                        Form("Matched Tracks v%d RP;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",i+2),
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0],
                        jetBins.size()-1, &jetBins[0]);

      hreff3D_vn[i] = f->make<TH3F>(Form("hreff3D_v%dRP",i+2),
                        Form("Matched Tracks v%d RP;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",i+2),
                        etaBins.size()-1, &etaBins[0],
                        ptBins.size()-1, &ptBins[0],
                        jetBins.size()-1, &jetBins[0]);

      hepplus[i] = f->make<TH1F>(Form("epplus_%d",i),Form("Order %d EP+ angle;#phi;events",i+2),100,-3.2,3.2);
      hepminus[i] = f->make<TH1F>(Form("epminus_%d",i),Form("Order %d EP- angle;#phi;events",i+2),100,-3.2,3.2);
      hepcomp[i] = f->make<TH2F>(Form("epcomp_%d",i),Form("Order %d EP+ vs EP- angle;#phi EP+;#phi EP-",i+2),100,-3.2,3.2,100,-3.2,3.2);   
      hepres[i] = f->make<TH2F>(Form("epres_%d",i),Form("Order %d EP cos( n(#Psi_{EP+} - #Psi_{EP-})); centrality; cosine",i+2),40,0.,40.,100,-1.0,1.0);   
 
    }


    // simulated 3D v42e
    hsim3D_v42e = f->make<TH3F>("hsim3D_v42e","Sim Tracks v42e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // reconstructed  3D v42e
    hrec3D_v42e = f->make<TH3F>("hrec3D_v42e","Rec Tracks v42e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


    // fakes 3D v42e
    hfak3D_v42e = f->make<TH3F>("hfak3D_v42e","Fake Tracks v42e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // matched 3D v42e
    heff3D_v42e = f->make<TH3F>("heff3D_v42e","Matched Tracks v42e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    hreff3D_v42e = f->make<TH3F>("hreff3D_v42e","Matched Tracks v42e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // simulated 3D v62e
    hsim3D_v62e = f->make<TH3F>("hsim3D_v62e","Sim Tracks v62e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // reconstructed  3D v62e
    hrec3D_v62e = f->make<TH3F>("hrec3D_v62e","Rec Tracks v62e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);


    // fakes 3D v62e
    hfak3D_v62e = f->make<TH3F>("hfak3D_v62e","Fake Tracks v62e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // matched 3D v62e
    heff3D_v62e = f->make<TH3F>("heff3D_v62e","Matched Tracks v62e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    hreff3D_v62e = f->make<TH3F>("hreff3D_v62e","Matched Tracks v62e;#eta;p_{T} (GeV/c);jet E_{T} (GeV/c)",
                      etaBins.size()-1, &etaBins[0],
                      ptBins.size()-1, &ptBins[0],
                      jetBins.size()-1, &jetBins[0]);

    // Individual centrality histograms storing cosine distributions.

    for( int c=0; c<40; c++)
    {


	    // simulated  v2
	    chsim_v2[c] = f->make<TH3F>(Form("chsim_v2_%d",c),"Sim Tracks v2;#eta;p_{T} (GeV/c);cosine",
			    etaBins.size()-1, &etaBins[0],
			    ptBins.size()-1, &ptBins[0],
			    cosBins.size()-1, &cosBins[0]);

	    // reconstructed   v2
	    chrec_v2[c] = f->make<TH3F>(Form("chrec_v2_%d",c),"Rec Tracks v2;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);


	    // fakes  v2
	    chfak_v2[c] = f->make<TH3F>(Form("chfak_v2_%d",c),"Fake Tracks v2;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // matcched  v2
	    cheff_v2[c] = f->make<TH3F>(Form("cheff_v2_%d",c),"Matcched Tracks v2;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    chreff_v2[c] = f->make<TH3F>(Form("chreff_v2_%d",c),"Matcched Tracks v2;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // simulated  v42
	    chsim_v42[c] = f->make<TH3F>(Form("chsim_v42_%d",c),"Sim Tracks v42;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // reconstructed   v42
	    chrec_v42[c] = f->make<TH3F>(Form("chrec_v42_%d",c),"Rec Tracks v42;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);


	    // fakes  v42
	    chfak_v42[c] = f->make<TH3F>(Form("chfak_v42_%d",c),"Fake Tracks v42;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // matcched  v42
	    cheff_v42[c] = f->make<TH3F>(Form("cheff_v42_%d",c),"Matcched Tracks v42;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    chreff_v42[c] = f->make<TH3F>(Form("chreff_v42_%d",c),"Matcched Tracks v42;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // simulated  v62
	    chsim_v62[c] = f->make<TH3F>(Form("chsim_v62_%d",c),"Sim Tracks v62;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // reconstructed   v62
	    chrec_v62[c] = f->make<TH3F>(Form("chrec_v62_%d",c),"Rec Tracks v62;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);


	    // fakes  v62
	    chfak_v62[c] = f->make<TH3F>(Form("chfak_v62_%d",c),"Fake Tracks v62;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // matcched  v62
	    cheff_v62[c] = f->make<TH3F>(Form("cheff_v62_%d",c),"Matcched Tracks v62;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    chreff_v62[c] = f->make<TH3F>(Form("chreff_v62_%d",c),"Matcched Tracks v62;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	     // RXN Plane vn

	    for( int i = 0; i<5; i++)
	    {
	      chsim_vn[i][c] = f->make<TH3F>(Form("chsim_v%dRP_%d",i+2,c),
				Form("Sim Tracks v%d RP;#eta;p_{T} (GeV/c);cosine",i+2),
				etaBins.size()-1, &etaBins[0],
				ptBins.size()-1, &ptBins[0],
				cosBins.size()-1, &cosBins[0]);

	      chrec_vn[i][c] = f->make<TH3F>(Form("chrec_v%dRP_%d",i+2,c),
				Form("Reco Tracks v%d RP;#eta;p_{T} (GeV/c);cosine",i+2),
				etaBins.size()-1, &etaBins[0],
				ptBins.size()-1, &ptBins[0],
				cosBins.size()-1, &cosBins[0]);

	      chfak_vn[i][c] = f->make<TH3F>(Form("chfak_v%dRP_%d",i+2,c),
				Form("Fake Tracks v%d RP;#eta;p_{T} (GeV/c);cosine",i+2),
				etaBins.size()-1, &etaBins[0],
				ptBins.size()-1, &ptBins[0],
				cosBins.size()-1, &cosBins[0]);

	      cheff_vn[i][c] = f->make<TH3F>(Form("cheff_v%dRP_%d",i+2,c),
				Form("Matcched Tracks v%d RP;#eta;p_{T} (GeV/c);cosine",i+2),
				etaBins.size()-1, &etaBins[0],
				ptBins.size()-1, &ptBins[0],
				cosBins.size()-1, &cosBins[0]);

	      chreff_vn[i][c] = f->make<TH3F>(Form("chreff_v%dRP_%d",i+2,c),
				Form("Matcched Tracks v%d RP;#eta;p_{T} (GeV/c);cosine",i+2),
				etaBins.size()-1, &etaBins[0],
				ptBins.size()-1, &ptBins[0],
				cosBins.size()-1, &cosBins[0]);
	 
	    }


	    // simulated  v42e
	    chsim_v42e[c] = f->make<TH3F>(Form("chsim_v42e_%d",c),"Sim Tracks v42e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // reconstructed   v42e
	    chrec_v42e[c] = f->make<TH3F>(Form("chrec_v42e_%d",c),"Rec Tracks v42e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);


	    // fakes  v42e
	    chfak_v42e[c] = f->make<TH3F>(Form("chfak_v42e_%d",c),"Fake Tracks v42e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // matcched  v42e
	    cheff_v42e[c] = f->make<TH3F>(Form("cheff_v42e_%d",c),"Matcched Tracks v42e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    chreff_v42e[c] = f->make<TH3F>(Form("chreff_v42e_%d",c),"Matcched Tracks v42e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // simulated  v62e
	    chsim_v62e[c] = f->make<TH3F>(Form("chsim_v62e_%d",c),"Sim Tracks v62e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // reconstructed   v62e
	    chrec_v62e[c] = f->make<TH3F>(Form("chrec_v62e_%d",c),"Rec Tracks v62e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);


	    // fakes  v62e
	    chfak_v62e[c] = f->make<TH3F>(Form("chfak_v62e_%d",c),"Fake Tracks v62e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    // matcched  v62e
	    cheff_v62e[c] = f->make<TH3F>(Form("cheff_v62e_%d",c),"Matcched Tracks v62e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
			      cosBins.size()-1, &cosBins[0]);

	    chreff_v62e[c] = f->make<TH3F>(Form("chreff_v62e_%d",c),"Matcched Tracks v62e;#eta;p_{T} (GeV/c);cosine",
			      etaBins.size()-1, &etaBins[0],
			      ptBins.size()-1, &ptBins[0],
                              cosBins.size()-1, &cosBins[0]);


    }


    // mom resolution (Sim to Rec) 
    hresStoR3D = f->make<TH3F>("hresStoR3D","Momentum resolution (sim to rec);#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_cent0 = f->make<TH3F>("hresStoR3D_cent0","Momentum resolution (sim to rec) 0-10 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c0 = f->make<TH3F>("hresStoR3D_c0","Momentum resolution (sim to rec) 0-5 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c1 = f->make<TH3F>("hresStoR3D_c1","Momentum resolution (sim to rec) 10-15 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c2 = f->make<TH3F>("hresStoR3D_c2","Momentum resolution (sim to rec) 20-30 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c3 = f->make<TH3F>("hresStoR3D_c3","Momentum resolution (sim to rec) 40-50 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c4 = f->make<TH3F>("hresStoR3D_c4","Momentum resolution (sim to rec) 70-80 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c5 = f->make<TH3F>("hresStoR3D_c5","Momentum resolution (sim to rec) 50-60 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);

    hresStoR3D_c6 = f->make<TH3F>("hresStoR3D_c6","Momentum resolution (sim to rec) 60-70 Centrality;#eta;sim p_{T} (GeV/c); #frac{p_{T,rec} - p_{T,sim}}{p_{T,rec}}  (GeV/c)",
                               etaBins.size()-1, &etaBins[0],
                               ptBins.size()-1, &ptBins[0],
                               resBins.size()-1, &resBins[0]);
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
HiPixelTrkEffHistograms::fillEPHistograms( float ep[2][5], float cent )
{
  for( int i=0;i<5;i++)
  {
    hepplus[i]->Fill(ep[0][i]);
    hepminus[i]->Fill(ep[1][i]);
    hepcomp[i]->Fill(ep[0][i],ep[1][i]);
    hepres[i]->Fill( cent, cos( (float (i+2)) * ( ep[0][i]-ep[1][i] ) ) );
  }
}

void 
HiPixelTrkEffHistograms::fillSimHistograms(const SimTrack_t & s)
{


  if(fillNtuples && s.status>0){
    simTrackValues = s;
    trackTrees[0]->Fill();
  }

  float eta = jet2010Mode ? fabs(s.etas) : s.etas ; 

  if(fillHistograms && s.status>0) {
    int c = (int)(floor(s.jetr));

    hsim->Fill(eta, s.pts);
    hsim3D->Fill(eta, s.pts, s.jetr);
    hsim3D_v2->Fill(eta, s.pts, s.jetr, cos(2*(s.phis - s.rxn)));
    hsim3D_v42->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rxn)));
    hsim3D_v62->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rxn)));
    chsim_v2[c]->Fill(eta, s.pts, cos(2*(s.phis - s.rxn)));
    chsim_v42[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rxn)));
    chsim_v62[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rxn)));
    if ( eta <= 0. )
    {
      hsim3D_v42e->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rp[0][0])));
      hsim3D_v62e->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rp[0][0])));
      for(int i=0;i<5;i++) hsim3D_vn[i]->Fill(eta, s.pts, s.jetr, cos((i+2)*(s.phis -s.rp[0][i])));
      chsim_v42e[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rp[0][0])));
      chsim_v62e[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rp[0][0])));
      for(int i=0;i<5;i++) chsim_vn[i][c]->Fill(eta, s.pts,  cos((i+2)*(s.phis -s.rp[0][i])));
    }
    if ( eta > 0. )
    {
      hsim3D_v42e->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rp[1][0])));
      hsim3D_v62e->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rp[1][0])));
      for(int i=0;i<5;i++) hsim3D_vn[i]->Fill(eta, s.pts, s.jetr, cos((i+2)*(s.phis -s.rp[1][i])));
      chsim_v42e[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rp[1][0])));
      chsim_v62e[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rp[1][0])));
      for(int i=0;i<5;i++) chsim_vn[i][c]->Fill(eta, s.pts,  cos((i+2)*(s.phis -s.rp[1][i])));
    }
    if(s.acc) 
    {
       hacc->Fill(eta, s.pts);
       hacc3D->Fill(eta,s.pts,s.jetr);
    }
    if(s.nrec==1) {
       hresStoR3D->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr);
       if ( s.jetr < 4.0 ) hresStoR3D_cent0->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr < 2.0 ) hresStoR3D_c0->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 4.0 && s.jetr < 6.0 ) hresStoR3D_c1->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 8.0 && s.jetr < 12.0 ) hresStoR3D_c2->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 16.0 && s.jetr < 20.0 ) hresStoR3D_c3->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 28.0 && s.jetr < 32.0 ) hresStoR3D_c4->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 20.0 && s.jetr < 24.0 ) hresStoR3D_c5->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       if ( s.jetr > 24.0 && s.jetr < 28.0 ) hresStoR3D_c6->Fill(eta, s.pts, (s.ptr-s.pts)/s.ptr );
       //if(fabs(s.etas)<1.0) hresStoR3D_etaS->Fill(s.jetr, s.pts, s.ptr);
       //if(fabs(s.etas)<2.4) hresStoR3D_etaL->Fill(s.jetr, s.pts, s.ptr);
    }
    if(s.nrec>0) 
    {
      heff->Fill(eta, s.pts); 
      heff3D->Fill(eta, s.pts, s.jetr);
      heff3D_v2->Fill(eta, s.pts, s.jetr, cos(2*(s.phis - s.rxn)));
      heff3D_v42->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rxn)));
      heff3D_v62->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rxn)));
      cheff_v2[c]->Fill(eta, s.pts, cos(2*(s.phis - s.rxn)));
      cheff_v42[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rxn)));
      cheff_v62[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rxn)));
      if ( eta <= 0. )
      {
        heff3D_v42e->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rp[0][0])));
        heff3D_v62e->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rp[0][0])));
        for(int i=0;i<5;i++) heff3D_vn[i]->Fill(eta, s.pts, s.jetr, cos((i+2)*(s.phis -s.rp[0][i])));
        cheff_v42e[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rp[0][0])));
        cheff_v62e[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rp[0][0])));
        for(int i=0;i<5;i++) cheff_vn[i][c]->Fill(eta, s.pts, cos((i+2)*(s.phis -s.rp[0][i])));
      }
      if ( eta > 0. )
      {
        heff3D_v42e->Fill(eta, s.pts, s.jetr, cos(4*(s.phis - s.rp[1][0])));
        heff3D_v62e->Fill(eta, s.pts, s.jetr, cos(6*(s.phis - s.rp[1][0])));
        for(int i=0;i<5;i++) heff3D_vn[i]->Fill(eta, s.pts, s.jetr, cos((i+2)*(s.phis -s.rp[1][i])));
        cheff_v42e[c]->Fill(eta, s.pts, cos(4*(s.phis - s.rp[1][0])));
        cheff_v62e[c]->Fill(eta, s.pts, cos(6*(s.phis - s.rp[1][0])));
        for(int i=0;i<5;i++) cheff_vn[i][c]->Fill(eta, s.pts, cos((i+2)*(s.phis -s.rp[1][i])));
      }

    }
    if(s.nrec>1) hmul->Fill(eta, s.pts), hmul3D->Fill(eta, s.pts, s.jetr);
  }


}

void 
HiPixelTrkEffHistograms::fillRecHistograms(const RecTrack_t & r)
{

  if(fillNtuples){
    recTrackValues = r;
    trackTrees[1]->Fill();
  }

  float eta = jet2010Mode ? fabs(r.etar) : r.etar ; 

  if(fillHistograms) {
    int c = (int)(floor(r.jetr));

    hrec->Fill(eta, r.ptr);
    hrec3D->Fill(eta, r.ptr, r.jetr);
    hrec3D_v2->Fill(eta, r.ptr, r.jetr, cos( 2*(r.phir - r.rxn)) );
    hrec3D_v42->Fill(eta, r.ptr, r.jetr, cos( 4*(r.phir - r.rxn)) );
    hrec3D_v62->Fill(eta, r.ptr, r.jetr, cos( 6*(r.phir - r.rxn)) );
    chrec_v2[c]->Fill(eta, r.ptr, cos( 2*(r.phir - r.rxn)) );
    chrec_v42[c]->Fill(eta, r.ptr, cos( 4*(r.phir - r.rxn)) );
    chrec_v62[c]->Fill(eta, r.ptr, cos( 6*(r.phir - r.rxn)) );
    if ( eta <= 0. )
    {
      hrec3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[0][0])));
      hrec3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[0][0])));
      for(int i=0;i<5;i++) hrec3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[0][i])));
      chrec_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[0][0])));
      chrec_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[0][0])));
      for(int i=0;i<5;i++) chrec_vn[i][c]->Fill(eta, r.ptr, cos((i+2)*(r.phir -r.rp[0][i])));
    }
    if ( eta > 0. )
    {
      hrec3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[1][0])));
      hrec3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[1][0])));
      for(int i=0;i<5;i++) hrec3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[1][i])));
      chrec_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[1][0])));
      chrec_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[1][0])));
      for(int i=0;i<5;i++) chrec_vn[i][c]->Fill(eta, r.ptr, cos((i+2)*(r.phir -r.rp[1][i])));
    }
    if(!r.nsim) 
    { 
      hfak->Fill(eta, r.ptr); 
      hfak3D->Fill(eta, r.ptr, r.jetr);
      hfak3D_v2->Fill( eta, r.ptr, r.jetr, cos( 2*(r.phir - r.rxn )) );
      hfak3D_v42->Fill( eta, r.ptr, r.jetr, cos( 4*(r.phir - r.rxn )) );
      hfak3D_v62->Fill( eta, r.ptr, r.jetr, cos( 6*(r.phir - r.rxn )) );
      chfak_v2[c]->Fill( eta, r.ptr,  cos( 2*(r.phir - r.rxn )) );
      chfak_v42[c]->Fill( eta, r.ptr, cos( 4*(r.phir - r.rxn )) );
      chfak_v62[c]->Fill( eta, r.ptr, cos( 6*(r.phir - r.rxn )) );
      if ( eta <= 0. )
      {
        hfak3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[0][0])));
        hfak3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[0][0])));
        for(int i=0;i<5;i++) hfak3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[0][i])));
        chfak_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[0][0])));
        chfak_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[0][0])));
        for(int i=0;i<5;i++) chfak_vn[i][c]->Fill(eta, r.ptr,  cos((i+2)*(r.phir -r.rp[0][i])));
      } 
      if ( eta > 0. )
      { 
        hfak3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[1][0])));
        hfak3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[1][0])));
        for(int i=0;i<5;i++) hfak3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[1][i])));
        chfak_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[1][0])));
        chfak_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[1][0])));
        for(int i=0;i<5;i++) chfak_vn[i][c]->Fill(eta, r.ptr,  cos((i+2)*(r.phir -r.rp[1][i])));
      }
    }
    if(r.nsim>0) 
    {
       hreff3D_v2->Fill( eta, r.ptr, r.jetr, cos( 2*(r.phir - r.rxn )) );
       hreff3D_v42->Fill( eta, r.ptr, r.jetr, cos( 4*(r.phir - r.rxn )) );
       hreff3D_v62->Fill( eta, r.ptr, r.jetr, cos( 6*(r.phir - r.rxn )) );
       chreff_v2[c]->Fill( eta, r.ptr,  cos( 2*(r.phir - r.rxn )) );
       chreff_v42[c]->Fill( eta, r.ptr,  cos( 4*(r.phir - r.rxn )) );
       chreff_v62[c]->Fill( eta, r.ptr,  cos( 6*(r.phir - r.rxn )) );
       hreff3D->Fill( eta, r.ptr, r.jetr ); 
      if ( eta <= 0. )
      {
        hreff3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[0][0])));
        hreff3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[0][0])));
        for(int i=0;i<5;i++) hreff3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[0][i])));
        chreff_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[0][0])));
        chreff_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[0][0])));
        for(int i=0;i<5;i++) chreff_vn[i][c]->Fill(eta, r.ptr,  cos((i+2)*(r.phir -r.rp[0][i])));
      }
      if ( eta > 0. )
      {
        hreff3D_v42e->Fill(eta, r.ptr, r.jetr, cos(4*(r.phir - r.rp[1][0])));
        hreff3D_v62e->Fill(eta, r.ptr, r.jetr, cos(6*(r.phir - r.rp[1][0])));
        for(int i=0;i<5;i++) hreff3D_vn[i]->Fill(eta, r.ptr, r.jetr, cos((i+2)*(r.phir -r.rp[1][i])));
        chreff_v42e[c]->Fill(eta, r.ptr,  cos(4*(r.phir - r.rp[1][0])));
        chreff_v62e[c]->Fill(eta, r.ptr,  cos(6*(r.phir - r.rp[1][0])));
        for(int i=0;i<5;i++) chreff_vn[i][c]->Fill(eta, r.ptr,  cos((i+2)*(r.phir -r.rp[1][i])));
      }

    }
    if(r.nsim>0 && r.status<0) hsec->Fill(eta, r.ptr), hsec3D->Fill(eta, r.ptr, r.jetr); // nsim>0 redudant?
  }

}

void 
HiPixelTrkEffHistograms::writeHistograms()
{

}



