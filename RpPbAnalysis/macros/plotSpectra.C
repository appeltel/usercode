#include "TrackCorrector.C"
#include "trackingPlotUtilities.C"

plotSpectra()
{

  

  double etaMin = 1.3 - 0.46;
  double etaMax = 1.8 - 0.46;
    
  double etaMinRev = -etaMax;
  double etaMaxRev = -etaMin;
    
    
  TrackCorrector t("../corrections/trackCorrections_HIN12017v4_HijingCombined.root");
  t.load("trkCorr_HIN12017");
  t.setOption1(true);
  t.setOption2(true);

  TrackCorrector tr("../corrections/trackCorrections_HIN12017v4_HijingCombined.root");
  tr.load("trkCorr_HIN12017_Pbp");
  tr.setOption1(true);
  tr.setOption2(true);

  TFile * fmb = new TFile("trackAnalysis_HIN12017_pPb_MB_spectra.root");
  TFile * fmbrev = new TFile("trackAnalysis_HIN12017_Pbp_MB_spectra.root");

  TFile * ftrk = new TFile("trackAnalysis_HIN12017_pPb_HighPt_spectra.root");
  TFile * ftrkrev = new TFile("trackAnalysis_HIN12017_Pbp_HighPt_spectra.root");

  // array goes as 0=MB, 1=Trk12, 2=Trk20, 3=Trk30

  TH3F *hTrk[4];
  TH3F *hTrkRev[4];
 
  hTrk[0] = (TH3F*) fmb->Get("trkAna_HIN12017/trkSpectrum");
  hTrkRev[0] = (TH3F*) fmbrev->Get("trkAna_HIN12017_Pbp/trkSpectrum");
  hTrk[1] = (TH3F*) ftrk->Get("trkAna_HIN12017_trk12/trkSpectrum");
  hTrkRev[1] = (TH3F*) ftrkrev->Get("trkAna_HIN12017_Pbp_trk12/trkSpectrum");
  hTrk[2] = (TH3F*) ftrk->Get("trkAna_HIN12017_trk20/trkSpectrum");
  hTrkRev[2] = (TH3F*) ftrkrev->Get("trkAna_HIN12017_Pbp_trk20/trkSpectrum");
  hTrk[3] = (TH3F*) ftrk->Get("trkAna_HIN12017_trk30/trkSpectrum");
  hTrkRev[3] = (TH3F*) ftrkrev->Get("trkAna_HIN12017_Pbp_trk30/trkSpectrum");

  TH1F * hnevt[4];
  TH1F * hnevtrev[4];
 
  hnevt[0] = (TH1F *) fmb->Get("trkAna_HIN12017/eventsLeadTrack");
  hnevt[1] = (TH1F *) ftrk->Get("trkAna_HIN12017_trk12/eventsLeadTrack");
  hnevt[2] = (TH1F *) ftrk->Get("trkAna_HIN12017_trk20/eventsLeadTrack");
  hnevt[3] = (TH1F *) ftrk->Get("trkAna_HIN12017_trk30/eventsLeadTrack");
  hnevtrev[0] = (TH1F *) fmbrev->Get("trkAna_HIN12017_Pbp/eventsLeadTrack");
  hnevtrev[1] = (TH1F *) ftrkrev->Get("trkAna_HIN12017_Pbp_trk12/eventsLeadTrack");
  hnevtrev[2] = (TH1F *) ftrkrev->Get("trkAna_HIN12017_Pbp_trk20/eventsLeadTrack");
  hnevtrev[3] = (TH1F *) ftrkrev->Get("trkAna_HIN12017_Pbp_trk30/eventsLeadTrack");
  
  // Calculate Nevt for each sample

  double nevt[4];
  double nevtrev[4];
    
  // Adding the nevent histogram
    
  for( int i=0;i<4;i++)
  {
    hnevt[i]->Add(hnevtrev[i],1.);
  }

  nevt[0] = hnevt[0]->Integral( hnevt[0]->GetXaxis()->FindBin(0.1),
                               hnevt[0]->GetXaxis()->FindBin(13.9) );
    
  nevt[1] = nevt[0] * hnevt[1]->Integral( hnevt[1]->GetXaxis()->FindBin(14.1),
                                         hnevt[1]->GetXaxis()->FindBin(21.9) )
            / hnevt[0]->Integral( hnevt[0]->GetXaxis()->FindBin(14.1),
                                  hnevt[0]->GetXaxis()->FindBin(21.9) );
  
  nevt[2] = nevt[1] * hnevt[2]->Integral( hnevt[2]->GetXaxis()->FindBin(22.1),
                                         hnevt[2]->GetXaxis()->FindBin(31.9) )
            / hnevt[1]->Integral( hnevt[1]->GetXaxis()->FindBin(22.1),
                                  hnevt[1]->GetXaxis()->FindBin(31.9) );
  nevt[3] = nevt[2] * hnevt[3]->Integral( hnevt[3]->GetXaxis()->FindBin(32.1),
                                          hnevt[3]->GetXaxis()->FindBin(999.9) )
            / hnevt[2]->Integral( hnevt[2]->GetXaxis()->FindBin(32.1),
                                  hnevt[2]->GetXaxis()->FindBin(999.9) );
    
/////////////////// Applying efficiency corrections first ///////////////////
    
  for( int i=0;i<4;i++)
  {
     applyTrackingCorrections( hTrk[i], t );
     applyTrackingCorrections( hTrkRev[i], tr );
  }

/////////////////////////////////////////////////////////////////////////////
    
  // Select Jet-Et Range


  // create 1D histos
  TH1D * hTrkPt[4];
  TH1D * hTrkPtRev[4];

  for( int i = 0; i<4; i++) 
  {
     hTrkPt[i] = getPtProjectedTH1D(Form("trkPt_%d",i),hTrk[i],etaMin,etaMax,0.0,999.0);
     hTrkPtRev[i] = getPtProjectedTH1D(Form("trkPtRev_%d",i),hTrkRev[i],etaMinRev,etaMaxRev,0.0,999.0);
  }

  double etaWid = etaMax - etaMin;

  TFile * fbcorr = new TFile("BinningAndResolutionCorrection_TrackTrigger.root");
  TH1F * hbincorr = (TH1F *) fbcorr->Get("hPt_pseudo2_copy1"); 

  for( int i=0; i<4; i++)
  {
    normalizeBinAndEtaWidth( hTrkPt[i], etaMax-etaMin );
    normalizeBinAndEtaWidth( hTrkPtRev[i], etaMaxRev-etaMinRev );
    applyCorrectionFactorRange( hTrkPt[i], hbincorr,1.0,999.0);
    applyCorrectionFactorRange( hTrkPtRev[i], hbincorr,1.0,999.0);
  }
    
    ///// Merging the pPb and Pbp samples together /////////////
  for( int i=0;i<4;i++)
  {
    hTrkPt[i]->Add(hTrkPtRev[i],1.);
    hTrkPt[i]->Scale(1./nevt[i]);
  }
     
     // add together spectra
     
  for( int i=1;i<4;i++)
  {
     hTrkPt[0]->Add(hTrkPt[i],1.);
  }
    
  ////////////////////////////////////////////////////////////
    
  hTrkPt[0]->SetMarkerStyle(20);

   TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   //c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.06);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

   //c1->Divide(1,2,0,0);
   c1->cd();
   c1->SetLogy();
   c1->SetLogx();
   c1->SetLogx();
  // c1->SetLogy();
  // c1->SetLogx();

   TH1F * hDum = new TH1F("hhdum","",10,0.4,110.);
   hDum->SetMinimum(1e-10);
   hDum->SetMaximum(100);
   hDum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum->GetXaxis()->SetTitleSize(0.035);
   hDum->GetXaxis()->SetTitleOffset(1.5);
   hDum->GetXaxis()->SetLabelSize(0.035);
   hDum->GetYaxis()->SetTitle("1/N_{evt} d^{2}N_{ch}/d_{T}d#eta");
   hDum->GetYaxis()->SetTitleSize(0.035);
   hDum->GetYaxis()->SetLabelSize(0.035);
   hDum->GetYaxis()->SetTitleOffset(1.5);

  hDum->Draw();
  hTrkPt[0]->Draw("same");

  TLatex * tex = new TLatex(0.6,1.0e-7,Form("%5.3f < #eta_{COM} < %5.3f",etaMin+0.46,etaMax+0.46));
  //TLatex * tex = new TLatex(2.6,1.05,Form(|#eta_{COM}| < 1.0");
  tex->SetTextSize(0.04);
  tex->Draw();


  TLegend * leg = new TLegend(0.5,0.63,0.9,0.9);
  leg->AddEntry(hTrkPt[0],"Corrected Tracks pPb Eta+","lp");
 // leg->AddEntry(hTrkPtRev[0],"Corrected Tracks Pbp","lp");
  leg->SetFillColor(kWhite);
  //leg->Draw();



}
