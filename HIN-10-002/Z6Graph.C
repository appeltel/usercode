void Z6Graph()
{
double temp;

double dndeta[4]= {1612,731.6, 266.5, 36.3 };

const int NCentBin=4;
const int NBinLow=20;
double xypt[NCentBin][NBinLow];
double xyY[NCentBin][NBinLow];
double xyYerr[NCentBin][NBinLow];
double xyYsys[NCentBin][NBinLow];

  gROOT->Reset();
  gROOT->SetStyle("MyStyle");
  gROOT->LoadMacro("MyDef.C");
  gStyle->SetTextFont(43);
  gStyle->SetLabelFont(43,"x");
  gStyle->SetLabelFont(43,"y");
  gStyle->SetLabelFont(43,"z");
  gStyle->SetTitleFont(43,"x");
  gStyle->SetTitleFont(43,"y");
  gStyle->SetTitleFont(43,"z");
  gStyle->SetEndErrorSize(0);
  int textsize = 23;

  int i=-1;
  cfn[++i]="./finaldata/Spectra_00_04_0.txt"   ;cft[i]="Yield";
  cfn[++i]="./finaldata/Spectra_00_04_4.txt"   ;cft[i]="";
  cfn[++i]="./finaldata/Spectra_00_04_8.txt"   ;cft[i]="";
  cfn[++i]="./finaldata/Spectra_00_04_11.txt"   ;cft[i]="";

  htit[0]="0-5%";
  htit[1]="20-25% (x2)";
  htit[2]="40-50% (x5)";
  htit[3]="70-80% (x30)";

  for (int ip=0;ip<Mpads;ip++)
  {
    hxmin[ip]=-0.2; hxmax[ip]=4.999;
    hymin[ip]=0.0001; hymax[ip]=2599;
    lxmin[ip]=0.4; lxmax[ip]=0.6;
    lymin[ip]=0.4; lymax[ip]=0.70;
    ptit[ip][0]=2.4; ptit[ip][1]=0.269*0.88*1500;
    atit[ip][0]="p_{T} (GeV/c)"; atit[ip][1]="1/N_{evt}d^{2}N/dp_{T}d#eta";
  }
  //
  // Read data and create vector files
  //
  for (int ifile=0;ifile<Mfiles;ifile++)
  {
    cout << ifile << "'th data" << endl;
    int np=0;
    ifstream infile;
    infile.open(cfn[ifile]);
    if (infile.is_open())
    {
      cout << "reading " << cfn[ifile] << endl;
      std::string line;
      while(std::getline(infile,line,'\n'))
      {
        double v1,v2,v3,v4,v5,v6,v7;
        int res = sscanf(line.c_str(),"%lg %lg %lg",&v1,&v2,&v3);
        if (res!=3) continue;
        if (v1 < 0.3 ) continue;
      //  cout <<"v1="<<v1<<" v2="<<v2<<" v3="<<v3 << endl;
        xa[ifile][np]=v1;
        xe[ifile][np]=0.0;
        ya[ifile][np]=v2;
        ye[ifile][np]=v3;
        // add systematic error
        //ye[ifile][np] = sqrt( v3*v3 + 0.06*0.06*v2*v2);

        xesys[ifile][np]=0.2;
        yesys[ifile][np]=0.;

        if( ifile == 1)
        {
          ya[ifile][np]=ya[ifile][np]*2.;
          ye[ifile][np]=ye[ifile][np]*2.;
          yesys[ifile][np]=yesys[ifile][np]*2.;
        }
        if( ifile == 2)
        {
          ya[ifile][np]=ya[ifile][np]*5.;
          ye[ifile][np]=ye[ifile][np]*5.;
          yesys[ifile][np]=yesys[ifile][np]*5.;
        }
        if( ifile == 3)
        {
          ya[ifile][np]=ya[ifile][np]*30.;
          ye[ifile][np]=ye[ifile][np]*30.;
          yesys[ifile][np]=yesys[ifile][np]*30.;
        }


        ye[ifile][np]=sqrt(ye[ifile][np]*ye[ifile][np]+yesys[ifile][np]*yesys[ifile][np]);

        np++;
      }
      infile.close();
    }
    else
    {
      //cout << ifile << "'th data not exist!: " << cfn[ifile] << endl;
    }
    ndp[ifile]=np;
  }
  //
  //
  //
ofstream xmpt;
ofstream xpara;
xmpt.open("ympt.txt");
xpara.open("zpara.txt");
ifstream inpart;
inpart.open("npart.txt");  double npartx[12];
for(int i=0;i<4;i++) inpart>>npartx[i];

  const int nw=4;
  TCanvas *can=new TCanvas("can","can",820,700);
  can->SetLeftMargin(0.2);
  can->Divide(2,2,0,0);
  TH1D *h1[nw];
  for (int iw=0;iw<nw;iw++)
  {
    can->cd(iw+1);
    char ch1[8];
    sprintf(ch1,"h1_%d",iw);
    h1[iw] = new TH1D(ch1,"",500,hxmin[iw],hxmax[iw]);
    h1[iw]->SetMinimum(hymin[iw]);
    h1[iw]->SetMaximum(hymax[iw]);
    h1[iw]->SetXTitle(atit[iw][0]);
    h1[iw]->SetYTitle(atit[iw][1]);
    h1[iw]->GetXaxis()->CenterTitle(1);
    h1[iw]->GetYaxis()->CenterTitle(1);
    h1[iw]->Draw();
  }
  cout << "Now Draw!" << endl;
  TGraphErrors *ge;
  for (int iw=0;iw<nw;iw++)
  {
    can->cd(iw+1);
  //gPad->SetLogy();
    //
    if (iw==0)
    {
      tex=new TLatex(0.7,ptit[iw][1]*0.86,"CMS Preliminary");
      tex->SetTextSize(20);
      tex->SetTextFont(43);
      //tex->Draw();
      tex1=new TLatex(0.7,ptit[iw][1]*0.74,"PbPb  #sqrt{s_{NN}} = 2.76TeV");
      tex1->SetTextFont(43);
      tex1->SetTextSize(20);
      //tex1->Draw();
      tex1=new TLatex(0.7,ptit[iw][1]*0.62,"|#eta| < 0.8");
      tex1->SetTextFont(43);
      tex1->SetTextSize(20);
      //tex1->Draw();

    }
    tex=new TLatex(ptit[iw][0],ptit[iw][1],htit[iw]);
    tex->SetTextFont(43);
    //tex->SetTextSize((iw>6)?0.07:0.09);

    tex->SetTextSize(20);

      h1[iw]->GetYaxis()->SetLabelSize(textsize);
      h1[iw]->GetYaxis()->SetTitleOffset(3.0);
      h1[iw]->GetYaxis()->SetTitleSize(textsize);



      h1[iw]->GetXaxis()->SetTitleOffset(2.5);
      h1[iw]->GetXaxis()->SetTitleSize(textsize);
      h1[iw]->GetXaxis()->SetLabelSize(textsize);
  //    h1[iw]->GetXaxis()->SetLabelOffset(0.001);


    tex->Draw();
    //
    // Legend
    //
    TLegend *leg = new TLegend(lxmin[iw],lymin[iw],lxmax[iw],lymax[iw]);
    leg->SetTextSize(18);
    leg->SetTextFont(43);
    //
    // Read data
    //
    for (int im=0;im<1;im++)
    {
      int j=im*4+iw;
      int col=mcol[im]+4;
      if (ndp[j]==0) continue;

      if(iw==3||iw==7||iw==11) gPad->SetRightMargin(0.005);

      TGraphErrors *ge = new TGraphErrors(ndp[j],&xa[j][0],&ya[j][0],&xe[j][0],&ye[j][0]);
      TGraphErrors *ge_sys = new TGraphErrors(ndp[j],&xa[j][0],&ya[j][0],&xesys[j][0],&yesys[j][0]);

      //sys err
      //ge_sys->SetMarkerStyle(20); 
      ge_sys->SetMarkerStyle(kFullCircle);
      ge_sys->SetMarkerColor(kBlack);
      ge_sys->SetMarkerSize(1.1);
      ge_sys->SetFillColor(kRed-9);
      ge_sys->Draw("PE2");

      //ge->SetMarkerStyle(20);
      ge->SetMarkerStyle(kFullCircle);
      ge->SetMarkerSize(1.1);
      ge->SetMarkerColor(kBlack);
      //ge->SetLineWidth(1.2);
      ge->SetLineColor(kBlack);
///
        //TF1 *fitptm = new TF1("fitptm","([0]-1.2)*([0]-2)/2/3.14159265359/[1]/[1]*[2]*pow((1+x/[1]),-[0])",0.,5.1);

  // Power law fit

        TF1 *fitptm1 = new TF1("fitptm1","2*3.14159265359*x*[2]*pow((1+x/[1]),-[0])",0.3,5.1);
        fitptm1->SetParameters(10.8327,2.221289,1201);

  // mT fit
        TF1 * fitptm2 = new TF1("fitptm2",ExpFuncScaled,0.3,0.7,3);
        fitptm2->SetParameters(0.2,0.,dndeta[iw]);
        //fitptm->FixParameter(2,dndeta[iw]);

  // Tsallis fit
  TF1 * fitptm3 = new TF1("fitptm3","6.28318531 * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.3,3.1);
  fitptm3->SetParameters(10.8327,2.221289,dndeta[iw]);

  TF1 * fitptmx1 = new TF1("fitptmx1","2*3.14159265359*x*x*[2]*pow((1+x/[1]),-[0])",0.0,3.1);
  TF1 * fitptmx2 = new TF1("fitptmx2",ExpFuncME,0.0,3.1,3);
  TF1 * fitptmx3 = new TF1("fitptmx3","6.28318531 * x * x * [2]*pow((1+sqrt(0.01947978 + x*x)-0.13957) /[1],-[0])",0.0,3.1);

        fitptm1->SetLineColor(kRed);
        fitptm1->SetLineWidth(2);
        fitptm1->SetLineStyle(1);
        ge->Fit("fitptm1","R0Q");
        fitptm2->SetLineColor(kBlue);
        fitptm2->SetLineWidth(2);
        fitptm2->SetLineStyle(1);
        ge->Fit("fitptm2","R0Q");
        fitptm3->SetLineColor(kGreen+2);
        fitptm3->SetLineWidth(2);
        fitptm3->SetLineStyle(1);
        ge->Fit("fitptm3","R0Q");
    fitptmx1->SetParameters(fitptm1->GetParameter(0),fitptm1->GetParameter(1), fitptm1->GetParameter(2) );
    fitptmx2->SetParameters(fitptm2->GetParameter(0),fitptm2->GetParameter(1), fitptm2->GetParameter(2) );
    fitptmx3->SetParameters(fitptm3->GetParameter(0),fitptm3->GetParameter(1), fitptm3->GetParameter(2) );
       
        float ndf=fitptm1->GetNDF();
        double chi2=fitptm1->GetChisquare();
        double chi2_1=fitptm1->GetChisquare() / fitptm1->GetNDF();
        double chi2_2=fitptm2->GetChisquare() / fitptm2->GetNDF();
        double chi2_3=fitptm3->GetChisquare() / fitptm3->GetNDF();
        cout << "fitptm3 chi2 = " << chi2_3 << endl;
        if(ndf==0) ndf=1;
       // cout<<"iw="<<iw<<", ndf= "<<ndf<<" #chi^2/ndf= "<<chi2/ndf<<endl;
       // cout<<fitptm->GetParameter(0)<<"       "<<fitptm->GetParameter(1)<<"       "<<fitptm->GetParameter(2)<<"   #chi^2/ndf=   "<<fitptm->GetChisquare()/fitptm->GetNDF()<<endl;
        //cout<<fitptm->GetParError(0)<<"       "<<fitptm->GetParError(1)<<"       "<<fitptm->GetParError(2)<<"   #chi^2/ndf=   "<<fitptm->GetChisquare()/fitptm->GetNDF()<<endl;
      //  xpara<<fitptm->GetParameter(0)<<"       "<<fitptm->GetParameter(1)<<"       "<<fitptm->GetParameter(2)<<"     "<<fitptm->GetChisquare()/fitptm->GetNDF()<<endl;

double xn=fitptm1->GetParameter(0);
double xp0=fitptm1->GetParameter(1);
double xne=fitptm1->GetParError(0);
double xp0e=fitptm1->GetParError(1);
//double salices=fitptm->Integral(0.3,5.0);

    double numerator, denominator, Mpt1, Mpt2, Mpt3;

    numerator = fitptmx1->Integral(0.0,0.3);
    denominator = fitptm1->Integral(0.0,0.3);
    for( int p = 0; p<ge->GetN(); p++)
    {
       Double_t  ppx, ppy;
       ge->GetPoint(p,ppx,ppy);
       if ( 0.3 < ppx && ppx < 0.6 ) numerator += ppx * ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) numerator += ppx * ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) numerator += ppx * ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) numerator += ppx * ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) numerator += ppx * ppy * 1.0;
       if ( 0.3 < ppx && ppx < 0.6 ) denominator += ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) denominator += ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) denominator += ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) denominator += ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) denominator += ppy * 1.0;
    }
    Mpt1 = numerator / denominator; 

    numerator = fitptmx2->Integral(0.0,0.3);
    denominator = fitptm2->Integral(0.0,0.3);
    for( int p = 0; p<ge->GetN(); p++)
    {
       Double_t  ppx, ppy;
       ge->GetPoint(p,ppx,ppy);
       if ( 0.3 < ppx && ppx < 0.6 ) numerator += ppx * ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) numerator += ppx * ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) numerator += ppx * ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) numerator += ppx * ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) numerator += ppx * ppy * 1.0;
       if ( 0.3 < ppx && ppx < 0.6 ) denominator += ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) denominator += ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) denominator += ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) denominator += ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) denominator += ppy * 1.0;
    }
    Mpt2 = numerator / denominator; 

    numerator = fitptmx3->Integral(0.0,0.3);
    denominator = fitptm3->Integral(0.0,0.3);
    for( int p = 0; p<ge->GetN(); p++)
    {
       Double_t  ppx, ppy;
       ge->GetPoint(p,ppx,ppy);
       if ( 0.3 < ppx && ppx < 0.6 ) numerator += ppx * ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) numerator += ppx * ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) numerator += ppx * ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) numerator += ppx * ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) numerator += ppx * ppy * 1.0;
       if ( 0.3 < ppx && ppx < 0.6 ) denominator += ppy * 0.1;
       if ( 0.6 < ppx && ppx < 1.2 ) denominator += ppy * 0.2;
       if ( 1.2 < ppx && ppx < 2.0 ) denominator += ppy * 0.4;
       if ( 2.0 < ppx && ppx < 4.0 ) denominator += ppy * 0.5;
       if ( 4.0 < ppx && ppx < 6.0 ) denominator += ppy * 1.0;
    }
    Mpt3 = numerator / denominator; 
      cout << denominator << endl;
      //ge->Draw("P");

      tex=new TLatex(ptit[iw][0]+0.5,ptit[iw][1]+1500,Form("<pT> = %1.3f",Mpt1));
      tex->SetTextFont(43);
      tex->SetTextColor(kRed);
      tex->SetTextSize(18);
      tex->Draw();
      tex=new TLatex(ptit[iw][0]+0.5,ptit[iw][1]+1200,Form("<pT> = %1.3f",Mpt2));
      tex->SetTextFont(43);
      tex->SetTextColor(kBlue);
      tex->SetTextSize(18);
      tex->Draw();
      tex=new TLatex(ptit[iw][0]+0.5,ptit[iw][1]+900,Form("<pT> = %1.3f",Mpt3));
      tex->SetTextFont(43);
      tex->SetTextColor(kGreen+2);
      tex->SetTextSize(18);
      tex->Draw();

      tex=new TLatex(ptit[iw][0],ptit[iw][1]+1500,Form("Chi2/ndof = %.1e",chi2_1));
      tex->SetTextFont(43);
      tex->SetTextColor(kRed);
      tex->SetTextSize(18);
      //tex->Draw();
      tex=new TLatex(ptit[iw][0],ptit[iw][1]+1300,Form("Chi2/ndof = %.1e",chi2_2));
      tex->SetTextFont(43);
      tex->SetTextColor(kBlue);
      tex->SetTextSize(18);
      //tex->Draw();
      tex=new TLatex(ptit[iw][0],ptit[iw][1]+1100,Form("Chi2/ndof = %.1e",chi2_3));
      tex->SetTextFont(43);
      tex->SetTextColor(kGreen+2);
      tex->SetTextSize(18);
      //tex->Draw();

      fitptm1->SetRange(0.0,5.1);
      fitptm1->Draw("SAME");
      fitptm2->SetRange(0.0,5.1);
      fitptm2->Draw("SAME");
      fitptm3->SetRange(0.0,5.1);
      fitptm3->Draw("SAME");

      ge->Draw("P");
      leg->AddEntry(ge,"Data","p");
      leg->AddEntry(fitptm1,"Power","l");
      leg->AddEntry(fitptm2,"mT-Exp","l");
      leg->AddEntry(fitptm3,"Tsallis","l");
    }
    if (iw==0) leg->Draw();
  }
  cout << "end of process" << endl;
  can->cd();


}


Double_t ExpFuncScaled(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t T = par[0];
  Double_t m = par[1];
  Double_t dNdy = par[2];

  Double_t mt = TMath::Sqrt( m*m + x[0]*x[0] );
  f = 1./(2.*TMath::Pi()*T*(T+m));

  f = f*dNdy;

  f = f*TMath::Exp(-(mt-m)/T);
   f*=2.*TMath::Pi()*x[0];

  return f;

}

Double_t ExpFuncME(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t T = par[0];
  Double_t m = par[1];
  Double_t dNdy = par[2];

  Double_t mt = TMath::Sqrt( m*m + x[0]*x[0] );
  f = 1./(2.*TMath::Pi()*T*(T+m));

  f = f*dNdy;

  f = f*TMath::Exp(-(mt-m)/T);
   f*=2.*TMath::Pi()*x[0]*x[0];

  return f;

}

