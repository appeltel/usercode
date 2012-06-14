#ifndef PixelTrackAnalysis_HiTPixelTrackEffHistograms
#define PixelTrackAnalysis_HiTPixelTrackEffHistograms

#include <vector>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

namespace edm { class ParameterSet ; }
class TFile;
class TTree;
class TH1F;
class TH2F;
class TH3F;

// define matched track structs

typedef struct
{
  int ids;
  float etas;
  float pts;
  float phis;
  int hits;
  int status;
  int acc;
  int nrec;
  float ptr;
  float dz;
  float d0;
  float pterr;
  float d0err;
  float dzerr;
  int hitr;
  int algo;
  float jetr;
  float rxn;
  float rp[2][5];
} SimTrack_t;

typedef struct
{
  int charge;
  float etar;
  float ptr;
  float phir;
  float dz;
  float d0;
  float pterr;
  float d0err;
  float dzerr;
  int hitr; 
  int algo;
  int nsim;
  int status;
  int ids;
  int parids;
  float etas;
  float pts;
  float jetr;
  float rxn;
  float rp[2][5];
} RecTrack_t;

// define class to hold histograms and track trees

class HiPixelTrkEffHistograms
{

 public:
   HiPixelTrkEffHistograms(const edm::ParameterSet& pset);
   ~HiPixelTrkEffHistograms();

  void declareHistograms();
  void fillSimHistograms(const SimTrack_t & s);
  void fillRecHistograms(const RecTrack_t & r);
  void fillEPHistograms( float  ep[2][5], float cent );
  void writeHistograms();

 private: 
   std::vector<TTree*> trackTrees;
   SimTrack_t simTrackValues;
   RecTrack_t recTrackValues;
   bool fillHistograms;
   bool fillNtuples;
   bool constPtBins;
   bool lowPtMode;
   bool flow2010Mode;
   bool jet2010Mode;
   bool flow2010ShengquanMode;
   edm::Service<TFileService> f;

   std::vector<double> etaBins, ptBins, jetBins, resBins, cosBins;

   // SimTrack
   TH2F* hsim;
   TH2F* hacc;
   TH2F* heff;
   TH2F* hmul;

   TH3F* hacc3D;
   TH3F* hsim3D;
   TH3F* heff3D;
   TH3F* hreff3D;
   TH3F* hmul3D;
   TH3F* hresStoR3D;
   TH3F* hresStoR3D_cent0;
   TH3F* hresStoR3D_c0;
   TH3F* hresStoR3D_c1;
   TH3F* hresStoR3D_c2;
   TH3F* hresStoR3D_c3;
   TH3F* hresStoR3D_c4;
   TH3F* hresStoR3D_c5;
   TH3F* hresStoR3D_c6;
   TH3F* hresStoR3D_etaS;
   TH3F* hresStoR3D_etaL;
   TH3F* hresStoR3D_v2;

   // uses the generator reaction plane

   TH3F* hsim3D_v2;
   TH3F* heff3D_v2;
   TH3F* hreff3D_v2;

   TH3F* hsim3D_v42;
   TH3F* heff3D_v42;
   TH3F* hreff3D_v42;
   TH3F* hsim3D_v62;
   TH3F* heff3D_v62;
   TH3F* hreff3D_v62;

   TH3F* chsim_v2[40];
   TH3F* cheff_v2[40];
   TH3F* chreff_v2[40];

   TH3F* chsim_v42[40];
   TH3F* cheff_v42[40];
   TH3F* chreff_v42[40];
   TH3F* chsim_v62[40];
   TH3F* cheff_v62[40];
   TH3F* chreff_v62[40];

   // calculates the reaction plane based on genparticles

   TH3F* hsim3D_vn[5];
   TH3F* heff3D_vn[5];
   TH3F* hreff3D_vn[5];

   TH3F* hsim3D_v42e;
   TH3F* heff3D_v42e;
   TH3F* hreff3D_v42e;

   TH3F* hsim3D_v62e;
   TH3F* heff3D_v62e;
   TH3F* hreff3D_v62e;

   TH3F* chsim_vn[5][40];
   TH3F* cheff_vn[5][40];
   TH3F* chreff_vn[5][40];

   TH3F* chsim_v42e[40];
   TH3F* cheff_v42e[40];
   TH3F* chreff_v42e[40];

   TH3F* chsim_v62e[40];
   TH3F* cheff_v62e[40];
   TH3F* chreff_v62e[40];


   // RecTrack
   TH2F* hrec;
   TH2F* hfak;
   TH2F* hsec;

   TH3F* hrec3D;
   TH3F* hfak3D;
   TH3F* hsec3D;

   TH3F* hrec3D_v2;
   TH3F* hfak3D_v2;
   TH3F* hrec3D_v42;
   TH3F* hfak3D_v42;
   TH3F* hrec3D_v42e;
   TH3F* hfak3D_v42e;
   TH3F* hrec3D_v62;
   TH3F* hfak3D_v62;
   TH3F* hrec3D_v62e;
   TH3F* hfak3D_v62e;

   TH3F* hrec3D_vn[5];
   TH3F* hfak3D_vn[5];

   TH3F* chrec_v2[40];
   TH3F* chfak_v2[40];
   TH3F* chrec_v42[40];
   TH3F* chfak_v42[40];
   TH3F* chrec_v42e[40];
   TH3F* chfak_v42e[40];
   TH3F* chrec_v62[40];
   TH3F* chfak_v62[40];
   TH3F* chrec_v62e[40];
   TH3F* chfak_v62e[40];

   TH3F* chrec_vn[5][40];
   TH3F* chfak_vn[5][40];


   TH1F* hepplus[5];
   TH1F* hepminus[5];
   TH2F* hepcomp[5];
   TH2F* hepres[5];

};

#endif

