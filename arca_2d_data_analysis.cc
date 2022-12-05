#include "TVirtualFitter.h"
#include "TFile.h"
#include "TLegend.h"
#include <TMarker.h>
#include "TApplication.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TColor.h"
#include <chrono> 

using namespace std;
using namespace std::chrono;

// 2D Gauss fit function
Double_t gaus2(Double_t *x, Double_t *par) {

   Double_t bkg_x = 1.+par[8]*x[0] + par[9]*x[0]*x[0] + par[10]*x[0]*x[0]*x[0];
   Double_t bkg_y = 1.+par[5]*x[1] + par[6]*x[1]*x[1]+par[7]*x[1]*x[1]*x[1];
   Double_t sig = par[3]; // Angular resolution
   Double_t amp = par[4];
   Double_t r1 = ((x[0]-par[1])/sig);
   Double_t r2 = ((x[1]-par[2])/sig);
   return par[0]*bkg_x*bkg_y*(1. - (amp*TMath::Pi()*0.26*0.26)/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
   //return par[0]*bkg_x*bkg_y*(1. - amp/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
}

// 2D Gauss fit function at [0,0]
Double_t gaus20(Double_t *x, Double_t *par) {

   Double_t bkg_x = 1.+par[6]*x[0] + par[7]*x[0]*x[0] + par[8]*x[0]*x[0]*x[0];
   Double_t bkg_y = 1.+par[3]*x[1] + par[4]*x[1]*x[1]+par[5]*x[1]*x[1]*x[1];
   //Double_t SizeMS = TMath::Pi()*0.26*0.26; //angular size of moon or sun
   Double_t sig = par[1]; // Angular resolution 
   //Double_t amp = SizeMS/(2.*TMath::Pi()*sig*sig);
   Double_t amp = par[2];
   //amp *= par[2];
   Double_t r1 = ((x[0])/sig);
   Double_t r2 = ((x[1])/sig);
   return par[0]*bkg_x*bkg_y*(1. - (amp*TMath::Pi()*0.26*0.26)/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
   //return par[0]*bkg_x*bkg_y*(1. - amp/(2.*TMath::Pi()*sig*sig)*TMath::Exp(-0.5*(r1*r1+r2*r2)));
}

// 2D Background fit function
Double_t bkg(Double_t *x, Double_t *par) {

  Double_t bkg_x = 1.+par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  Double_t bkg_y = 1.+par[1]*x[1] + par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1];
  return par[0]*bkg_x*bkg_y;
  //return par[0]+par[1]*x[1]+par[2]*x[1]*x[1]+par[3]*x[1]*x[1]*x[1]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0];
}



int main(int argc, char** argv) {
  
  if (argc < 5) {
    cout << "5 arguments are required: <input_file> <output_plots_name> <beta_cut> <lambda_map_max_degree> <lambda_map_binning> "<<endl;
    exit(1);
  }

  
  auto start = high_resolution_clock::now();


  //TApplication app("app",NULL,NULL); 
  TFile *file=new TFile("data_sun.root", "RECREATE");
  
  string filename, daz, dalt, fitinf0, sundist, MC_sun_dist, lik, len;
  double chi;
  int nHitFit;
 
 double sig = 0.5; //0.25 for ARCA8/ARCA 0.26 for ARCA6
 //double beta_cut = 0.32;//<----------------------------------
 double beta_cut = std::stof(argv[3]);
 //double beta_cut = 0.32; //<---------- 0.31 for comb. ARCA6, 0.32 for ARCA8 and ARCA
 //cout<<"sig = "<<sig<<endl;
 //cout<<"beta_cut = "<<beta_cut<<endl;
 float size = 0.05; // label and title size for various axes

// ---------- define 2D histogram of event number distribution

 double MaxGrad = std::stof(argv[4]);
 double Bin = 0.1;
 int Nbin2 = (2.*MaxGrad + Bin/2. ) / Bin;
 TH2D *h2_noshad = new TH2D("h2_noshad","",Nbin2,-MaxGrad,MaxGrad,Nbin2,-MaxGrad,MaxGrad);
 h2_noshad->SetBinErrorOption(TH2::kPoisson);
 h2_noshad->SetStats(kFALSE);
 h2_noshad->SetXTitle("Azimuth [deg]");
 h2_noshad->SetYTitle("Altitude [deg]");
 h2_noshad->GetXaxis()->SetTitleSize(size);
 h2_noshad->GetXaxis()->SetTitleOffset(0.9);  
 h2_noshad->GetXaxis()->SetTitleFont(62);  
 h2_noshad->GetXaxis()->SetLabelFont(62);  
 h2_noshad->GetXaxis()->SetLabelSize(size);
 h2_noshad->GetYaxis()->SetTitleSize(size);
 h2_noshad->GetYaxis()->SetTitleOffset(0.9);  
 h2_noshad->GetYaxis()->SetTitleFont(62);  
 h2_noshad->GetYaxis()->SetLabelFont(62);  
 h2_noshad->GetYaxis()->SetLabelSize(size);

 TH2D *h2_shad = new TH2D("h2_shad","",Nbin2,-MaxGrad,MaxGrad,Nbin2,-MaxGrad,MaxGrad);
 h2_shad->SetBinErrorOption(TH2::kPoisson);
 h2_shad->SetStats(kFALSE);
 h2_shad->SetXTitle("Azimuth [deg]");
 h2_shad->SetYTitle("Altitude [deg]");
 h2_shad->GetXaxis()->SetTitleSize(size);
 h2_shad->GetXaxis()->SetTitleOffset(0.9);
 h2_shad->GetXaxis()->SetTitleFont(62);
 h2_shad->GetXaxis()->SetLabelFont(62);
 h2_shad->GetXaxis()->SetLabelSize(size);
 h2_shad->GetYaxis()->SetTitleSize(size);
 h2_shad->GetYaxis()->SetTitleOffset(0.9);
 h2_shad->GetYaxis()->SetTitleFont(62);
 h2_shad->GetYaxis()->SetLabelFont(62);
 h2_shad->GetYaxis()->SetLabelSize(size);


// ---------- define 1D projection of event distribution to obtain BG Altitude fit

 TH1D *hAlt = new TH1D("hAlt","",Nbin2,-MaxGrad,MaxGrad);
 //hAlt->SetStats(kFALSE);
 hAlt->SetXTitle("Altitude [deg]");
 hAlt->SetYTitle("Events");
 hAlt->GetXaxis()->SetTitleSize(size);
 hAlt->GetXaxis()->SetTitleOffset(0.9);  
 hAlt->GetXaxis()->SetTitleFont(62);  
 hAlt->GetXaxis()->SetLabelFont(62);  
 hAlt->GetXaxis()->SetLabelSize(size);
 hAlt->GetYaxis()->SetTitleSize(size);
 hAlt->GetYaxis()->SetTitleOffset(0.9);  
 hAlt->GetYaxis()->SetTitleFont(62);  
 hAlt->GetYaxis()->SetLabelFont(62);  
 hAlt->GetYaxis()->SetLabelSize(size);

// ---------- define 1D projection of event distribution to obtain BG Azimuth fit

 TH1D *hAzi = new TH1D("hAzi","",Nbin2,-MaxGrad,MaxGrad);
 //hAzi->SetStats(kFALSE);
 hAzi->SetXTitle("Azimuth [deg]");
 hAzi->SetYTitle("Events");
 hAzi->GetXaxis()->SetTitleSize(size);
 hAzi->GetXaxis()->SetTitleOffset(0.9);  
 hAzi->GetXaxis()->SetTitleFont(62);  
 hAzi->GetXaxis()->SetLabelFont(62);  
 hAzi->GetXaxis()->SetLabelSize(size);
 hAzi->GetYaxis()->SetTitleSize(size);
 hAzi->GetYaxis()->SetTitleOffset(0.9);  
 hAzi->GetYaxis()->SetTitleFont(62);  
 hAzi->GetYaxis()->SetLabelFont(62);  
 hAzi->GetYaxis()->SetLabelSize(size);
 
 // ---------- loop over n_beta  --------------------------

 cout << "Analyzing: "<<argv[1]<<endl;
 cout << "2D plot will be saved as: "<<argv[2]<<endl; 
 cout << "Beta0 cut = "<<argv[3]<<endl;
 cout << "Maxdeg = "<<argv[4]<<endl;
 cout << "Bin size = "<<argv[5]<<endl;

 int n_h1 =0;
 double p_nomoon;
 double chi2Fit0;
 // ---------- loop over ntuple muon events --------------------------
 double dtr = TMath::DegToRad();
 ifstream fstream;
   
   //fstream.open("csv/arca_data_9635-11707_pre+anoise+combined_shadow.txt");
   //fstream.open("csv/arca19/data/arca19_data_pre+anoise+sun.txt");
   //fstream.open("csv/arca19/data/arca19_data_pre+anoise+moon_fake_-8.txt");
   fstream.open(argv[1]);
   //fstream.open("csv/arca6/arca6_data_9635-10286_pre+anoise+combined_shadow.txt");

 {
      while ( ! fstream.eof() )
    {
      fstream >> daz >> dalt >> fitinf0 >> sundist >> lik >> len;
     double DeltaAzi = std::stof(daz);
     double DeltaAlt = std::stof(dalt);
     double beta0 = std::stof(fitinf0)*TMath::RadToDeg();
     double sun_dist = std::stof(sundist);
     double likelihood = std::stof(lik);
     double track_length = std::stof(len);
     
     if ( beta0 < beta_cut ) // <------------------------------------------------------ CUT ----------------
       { 
	 if (DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad && DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad ) {
	     hAzi->Fill(DeltaAzi,1./Nbin2);
	 }
	 if (DeltaAzi < MaxGrad &&  DeltaAzi > -MaxGrad && DeltaAlt < MaxGrad && DeltaAlt > -MaxGrad ) {
	     hAlt->Fill(DeltaAlt,1./Nbin2);
	 }
	 h2_shad->Fill(DeltaAzi,DeltaAlt); 
	 n_h1++;
	 
       }

    }
 }// end event loop
 fstream.close();
  
  // --------------- event loop finished ------------ make fits !! 
 cout<<"N. of events in S1 = "<<h2_shad->GetEntries()<<endl;
 cout<<"N. of events in x distribution = "<<hAzi->GetEntries()<<endl;
 cout<<"N. of events in y distribution = "<<hAlt->GetEntries()<<endl;


   // fit projection of 2D histogram to correct for changing event rate
 TF1 *pol0 = new TF1("pol0", "pol0",-MaxGrad,MaxGrad);
 TF1 *pol2 = new TF1("pol2", "pol2",-MaxGrad,MaxGrad);
 TF1 *pol3 = new TF1("pol3", "pol3",-MaxGrad,MaxGrad);

 cout<<"Fit minchi2 Azimuth: "<<endl;
 hAzi->Fit("pol0","N");
 double a0 = pol0->GetParameter(0);
 double x1 = 0;
 double x2 = 0;
 double x3 = 0;
 //double x1 = pol3->GetParameter(1)/a0;
 //double x2 = pol3->GetParameter(2)/a0;
 //double x3 = pol3->GetParameter(3)/a0;
 
 double chi2azi = pol0->GetChisquare(); //<------------ no need of multiplying by 2 as the minimum fcn value is already the pearson chi2 value!! 
 


 hAlt->Fit("pol2","N");
 a0 = pol2->GetParameter(0);
 double y1 = pol2->GetParameter(1)/a0;
 double y2 = pol2->GetParameter(2)/a0;
 double y3 = 0; 
//double y3 = pol3->GetParameter(3)/a0;
 double chi2alt = pol2->GetChisquare();
  
 cout << "--------- 1D FITS --------  "<<endl;
 cout << "--------------------------  "<<endl;
 cout << "chi2azi_ = " << chi2azi<<endl;
 cout << "chi2alt = " << chi2alt<<endl;
 cout << "p(azi) = " << TMath::Prob(chi2azi,Nbin2-(4))<<endl;
 cout << "p(alt) = " << TMath::Prob(chi2alt,Nbin2-(4))<<endl;



// define 2d fit function, initialize parameters

 //gROOT->SetStyle("Plain");
 gStyle->SetOptFit(1111);
 gStyle->SetOptStat("nemr");


 //--------------PLOT MARGINAL DISTRIBUTIONS-------------------
 TCanvas *c0 = new TCanvas("c0", "Azimuth", 200, 10, 600, 400);
 //c0->cd();
 //hAzi->Draw("H");
 //hAzi->Draw("E, P, SAME");
 //c0->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/data_azimuth_distr.pdf");
 

 TCanvas *c1 = new TCanvas("c1", "Altitude", 200, 10, 600, 400);
 //c1->cd();
 //hAlt->Draw("E, P, SAME");
 //hAlt->Draw("H");
 /*
 gPad->Update();
 auto stat = dynamic_cast<TPaveStats*>(hAlt->FindObject("stats"));
 if (stat) {
   //cout << " X1NDC: " << stat->GetX1NDC() << " X2NDC: " << stat->GetX2NDC() << endl;
   stat->SetX1NDC(0.1); stat->SetX2NDC(0.3);
   stat->Draw();
 } else {
   cout << "No stats box found!\n";
 }
 */
 //c1->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/data_altitude_distr.pdf");

 //--------------PLOT 2D DISTRIBUTION--------------------

 TCanvas *c2 = new TCanvas("c3", "shad_distribution", 200, 10, 600, 400);
 //c2->cd();
 //h2_shad->Draw("colz");
 //c2->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/data_2D_distr.pdf");


 //TH1D *proj_x = h2_shad->ProjectionX();
  
 
 const Int_t npar = 11;
 double xs =  0.01;
 double ys = 0.01;
 double amp = 1.;
 double SigFit = sig; // <-----------------------------------------------------------
    
 Double_t f2params[npar] = {a0,xs,ys,SigFit,amp,y1,y2,y3,x1,x2,x3};
    TF2 *f2 = new TF2("f2",gaus2,-MaxGrad,MaxGrad,-MaxGrad,MaxGrad,npar);
    f2->SetParameters(f2params);
    
    //f2->SetParameter(3,SigFit);
    //f2->SetParameters(100,0.1,0.1,0.1,1.,0.001,0.001,0.001,0.001,0.001,0.001);
    f2->SetParNames("a0","XS","YS","sigma","amp","y1","y2","y3","x1","x2","x3");    
    f2->SetParLimits(4,0.,10.);
    f2->SetParLimits(3,0.,1.00);
    f2->FixParameter(3,SigFit);    
    f2->FixParameter(7,0);
    f2->FixParameter(8,0);
    f2->FixParameter(9,0);
    f2->FixParameter(10,0);
    f2->SetParLimits(1,-1.,1.);
    f2->SetParLimits(2,-1.,1.);


    const Int_t np20 = 9;
    Double_t f20params[np20] = {a0,SigFit,amp,y1,y2,y3,x1,x2,x3};
    TF2 *f20 = new TF2("f20",gaus20,-MaxGrad,MaxGrad,-MaxGrad,MaxGrad, np20);
    f20->SetParameters(f20params);
    //f20->SetParameters(100,0.1,1.,0.001,0.001,0.001,0.001,0.001,0.001);
    f20->SetParNames("a0","sigma","amp","y1","y2","y3","x1","x2","x3");    
    f20->SetParLimits(2,0.,10.);
    //f20->SetParameter(1,SigFit);
    f20->SetParLimits(1,0.1,1.00);
    f20->FixParameter(5,0);
    f20->FixParameter(6,0);
    f20->FixParameter(7,0);
    f20->FixParameter(8,0);


    const Int_t npbkg = 7;
    Double_t b2params[npbkg] = {a0,y1,y2,y3,x1,x2,x3};
    TF2 *bkg2 = new TF2("bkg2",bkg,-MaxGrad,MaxGrad,-MaxGrad,MaxGrad, npbkg);
    bkg2->SetParameters(b2params);
    bkg2->SetParNames("a0","y1","y2","y3","x1","x2","x3");
    bkg2->FixParameter(3,0);
    bkg2->FixParameter(4,0);
    bkg2->FixParameter(5,0);
    bkg2->FixParameter(6,0);


   // ---------- 2Dim Fit on h2 histogram of event distribution
    cout << "Before 2D fits... " << endl;
    
    // Full fit of best MS position
    h2_shad->Fit("f2","BNL"); // "B" option if fixed Sigma !
    double chi2Fit = f2->GetChisquare();
    double Norm   = f2->GetParameter(0);
    double XS = f2->GetParameter(1);
    double YS = f2->GetParameter(2);
    double Sig = f2->GetParameter(3);           
           amp    = f2->GetParameter(4);
    double A1     = f2->GetParameter(5);
    double A2     = f2->GetParameter(6);
    double A3     = f2->GetParameter(7);
    double B1     = f2->GetParameter(8);
    double B2     = f2->GetParameter(9);
    double B3     = f2->GetParameter(10);
    

    //chi2Fit = f2->GetChisquare();
    //cout << "MS fit at best pos gaus = " << chi2Fit <<endl;
    
    TVirtualFitter *fitter2 = TVirtualFitter::GetFitter();
    Double_t amin2,edm,errdef;
    Int_t nvpar,nparx;
    
    fitter2->GetStats(amin2,edm,errdef,nvpar,nparx);
    chi2Fit = amin2*2;
    cout << "MS fit best pos. LogL_chi2  = "<< amin2*2 << endl;


    // fit of MS at nominal position
    h2_shad->Fit("f20","NL");
    double amp0 = f20->GetParameter(2);
    //chi2Fit0 = f20->GetChisquare();
    //cout << "MS fit gaus at nominal pos. = " << chi2Fit0 <<endl;
    double Sig0 = f20->GetParameter(1);


    TVirtualFitter *fitter20 = TVirtualFitter::GetFitter();
    Double_t amin20;
    

    fitter20->PrintResults(2,0.);
    fitter20->GetStats(amin20,edm,errdef,nvpar,nparx);
    chi2Fit0 = amin20*2;
    cout << "MS fit at nominal position (XS=0, YS=0) LogL_chi2  = "<< amin20*2 <<endl;



    // background only fit
    
    h2_shad->Fit("bkg2","NL");
    double chi2Fitbkg = bkg2->GetChisquare();
    TVirtualFitter *fitterbkg = TVirtualFitter::GetFitter();
    Double_t aminbkg;
    fitterbkg->GetStats(aminbkg,edm,errdef,nvpar,nparx);
    //cout << "MS fit bkg classical chi  = "<< chi2Fitbkg <<endl;
    chi2Fitbkg = aminbkg*2;
    cout << "MS fit bkg logL_chi2  = " << aminbkg*2 <<endl;

    // ---------- output of various delta chi2
    double lambda_0 = chi2Fit0 - chi2Fitbkg;
    double lambda_min = chi2Fit - chi2Fitbkg;
    double theta = lambda_0 - lambda_min;

    cout << "Events: "<<n_h1<<endl;
    cout << "Norm [sq deg] = " << Norm/Bin/Bin << endl;
    cout << "--------- BEST SHADOW POSITION FIT (XS, YS) --------  "<<endl;
    cout << "Sigma resolution: "<< f2->GetParameter(3)<<"  with error: "<<f2->GetParError(3)<<endl;
    cout << "amp = " << amp <<" +- "<<f2->GetParError(4)<< endl;
    cout << "Chi2 (best fit 0.1 deg bin, +-3 grad) = " << chi2Fit <<endl;
    cout << "XS = " << XS << endl;
    cout << "YS =  " << YS << endl;
    cout << "P-value with classical chi2 method = " <<f2->GetProb()<< endl;
    cout << "f2 chi2 / NDOF = " << chi2Fit << " / "<< f2->GetNDF()<<endl;
    cout << "p(chi2f2) = " << TMath::Prob(chi2Fit,f2->GetNDF())<<endl;
    cout << "--------- NOMINAL POSITION FIT (XS=0, YS=0) --------  "<<endl;
    cout << "Sigma resolution: "<< Sig0<<" +- "<<f20->GetParError(1)<<endl; // comment  if sigma fixed for nominal fit
    cout << "amp(0,0) = " << amp0 <<" +- "<<f20->GetParError(2)<< endl;
    cout << "P-value with classical chi2 method = " <<f20->GetProb()<< endl;
    cout << "f20 chi2 / NDOF = " << chi2Fit0 << " / "<< f20->GetNDF()<<endl;
    cout << "p(chi2f20) = " << TMath::Prob(chi2Fit0,f20->GetNDF())<<endl;
    cout << "--------- BACKGROUND ONLY FIT ----------  "<<endl;
    cout << "P-value with classical chi2 method = " << bkg2->GetProb()<<endl;
    cout << "Bkg chi2 / NDOF = " << chi2Fitbkg << " / "<< bkg2->GetNDF()<<endl;
    cout << "p(chi2bkg) = " << TMath::Prob(chi2Fitbkg,bkg2->GetNDF())<<endl;
    cout << "------------- LAMBDA -------------  "<<endl;
    cout << "lambda_0  = " << lambda_0 <<endl;
    cout << "lambda_min  = " << lambda_min <<endl;
    cout << "Theta = " << theta <<endl;
    cout << "p-value(0,0) = " << TMath::Prob(-lambda_0,2) <<endl;// 2 becomes 1 if for nominal fit sigma is fixed
    cout << "sigma = "<<sqrt(2)*TMath::ErfInverse(1.-TMath::Prob(-lambda_0,2)) << endl; // 2 becomes 1 if for nominal fit sigma is fixed
    cout << "p-value(theta) = " << TMath::Prob(theta,2) <<endl;// 2 becomes 1 if for nominal fit sigma is fixed
    cout << "sigma = "<<sqrt(2)*TMath::ErfInverse(1.-TMath::Prob(theta,2)) << endl; // 2 becomes 1 if for nominal fit sigma is Fixed
     
    

    //double MaxGrad2 = 3;//6;
    //double BinChi2 = 0.5;//0.1;
    //int NbinChi2 = (2.*MaxGrad2 + BinChi2/2. ) / BinChi2;

    double max_deg = std::stof(argv[4]);
    double BinChi2 = std::stof(argv[5]);
    cout<<"binchi2 = "<<BinChi2<<endl;
    double MinGrad2 = -max_deg+BinChi2;
    double MaxGrad2 = +max_deg+BinChi2;
    int NbinChi2 = (MaxGrad2 - MinGrad2)/BinChi2;
    cout<<"NbinChi2 = "<<NbinChi2<<endl;

    TH2D *hChi2 = new TH2D("hChi2","hChi2",NbinChi2,MinGrad2,MaxGrad2,NbinChi2,MinGrad2,MaxGrad2);


    //TH2D *hChi2 = new TH2D("hChi2","",NbinChi2,-MaxGrad2,MaxGrad2,NbinChi2,-MaxGrad2,MaxGrad2);

    hChi2->SetStats(kFALSE);
    hChi2->SetXTitle("x [deg]");
    hChi2->SetYTitle("y [deg]");
    hChi2->GetXaxis()->SetTitleSize(size);
    hChi2->GetXaxis()->SetTitleOffset(0.9);
    hChi2->GetXaxis()->SetTitleFont(62);
    hChi2->GetXaxis()->SetLabelFont(62);
    hChi2->GetXaxis()->SetLabelSize(size);
    hChi2->GetYaxis()->SetTitleSize(size);
    hChi2->GetYaxis()->SetTitleOffset(0.9);
    hChi2->GetYaxis()->SetTitleFont(62);
    hChi2->GetYaxis()->SetLabelFont(62);
    hChi2->GetYaxis()->SetLabelSize(size);
    hChi2->GetXaxis()->SetRangeUser(-max_deg,max_deg);
    hChi2->GetYaxis()->SetRangeUser(-max_deg,max_deg);
    hChi2->SetMinimum(0.);

    gPad->SetGridx();
    gPad->SetGridy();
    //gStyle->SetPalette(kLightTemperature);
    //TColor::InvertPalette();
    // --------------- fill the delta chi2 histogram from refitting amplitude
    
    for (int i=0;i<NbinChi2;i++) { // loop over x-axis in chi2 histogram
      for (int j=0;j<NbinChi2;j++) { // loop over y-axis in chi2 histogram
	int ibin = hChi2->GetBin(i,j);
	double XS = hChi2->GetXaxis()->GetBinCenter(i);
	double YS = hChi2->GetYaxis()->GetBinCenter(j);
	amp = 0.; // amplitude start value zero
	f2->SetParameter(0,Norm);
	f2->FixParameter(1,XS);
	f2->FixParameter(2,YS);
	//f2->FixParameter(3,SigFit);
	f2->FixParameter(3,Sig0);
	f2->SetParameter(4,amp);
	f2->SetParameter(5,A1);
	f2->SetParameter(6,A2);// <----------------------------------------- here change fix/set for real lambda
	f2->FixParameter(7,0);
	f2->FixParameter(8,0);
	f2->FixParameter(9,0);
	f2->FixParameter(10,0);
	f2->SetParLimits(4,0.,2.);
	//f2->SetParLimits(3,0.,1.00);
	h2_shad->Fit("f2","QBNL");
	double chi2; // = f2->GetChisquare();
	TVirtualFitter *fitter_chi = TVirtualFitter::GetFitter();
	Double_t amin_chi;
	fitter_chi->GetStats(amin_chi,edm,errdef,nvpar,nparx);
	chi2 = amin_chi*2;
	Double_t lam;
	lam = chi2-chi2Fitbkg;
	//lam = chi2-chi2Fit; -->profile likelihood
	//hChi2->SetBinContent(ibin,chi2-chi2Fit);
	hChi2->SetBinContent(ibin,-lam);
	
	//cout<<"XS = "<<XS<<" YS = "<<YS<<" l = "<<lam << " amp = " << f2->GetParameter(4) << " +- " << f2->GetParError(4) << " sigma = " << f2->GetParameter(3) << " +- " << f2->GetParError(3) << endl;
      }
    }

    /*
    Int_t min_bin = hChi2->GetMinimumBin();
    Int_t x,y,z;
    Double_t lamb_min=hChi2->GetBinContent(hChi2->GetMinimumBin());
    hChi2->GetBinXYZ(min_bin, x, y, z);
    cout << "The bin having the minimum value is "<< hChi2->GetXaxis()->GetBinCenter(x)<<","<<hChi2->GetXaxis()->GetBinCenter(y)<<" with "<< lamb_min <<endl;
    double th = lambda_0 - lamb_min;
    */

    //Int_t min_bin = hChi2->GetMinimumBin();
    Int_t max_bin = hChi2->GetMaximumBin();
    Int_t x,y,z;
    //Double_t lamb_min=hChi2->GetBinContent(hChi2->GetMinimumBin());
    Double_t lamb_max=hChi2->GetBinContent(hChi2->GetMaximumBin());
    Double_t lamb_min = -lamb_max;
    //hChi2->GetBinXYZ(min_bin, x, y, z);
    hChi2->GetBinXYZ(max_bin, x, y, z);
    //cout << "The bin having the minimum value is "<< hChi2->GetXaxis()->GetBinCenter(x)<<","<<hChi2->GetXaxis()->GetBinCenter(y)<<" with "<< lamb_min <<endl;
    Double_t x_max=hChi2->GetXaxis()->GetBinCenter(x);
    Double_t y_max=hChi2->GetYaxis()->GetBinCenter(y);
    //cout << "The bin having the maximum value is "<< hChi2->GetXaxis()->GetBinCenter(x)<<","<<hChi2->GetYaxis()->GetBinCenter(y)<<" with "<< lamb_max <<endl;
    double th = lambda_0 - lamb_min;
    

    cout << "lambda_minimum_bin  = " << lamb_min <<endl;
    cout << "Theta = " << th <<endl;
    cout << "p-value(theta) = " << TMath::Prob(th,1) <<endl;// 2 becomes 1 if for nominal fit sigma is fixed
    cout << "sigma = "<<sqrt(2)*TMath::ErfInverse(1.-TMath::Prob(th,1)) << endl; // 2 becomes 1 if for nominal fit sigma is fixed
    
    
    TCanvas *c4 = new TCanvas("c4", "lambda_map", 200, 10, 600, 400);
    c4->cd();
    hChi2->DrawCopy("colz0");
    
    TCanvas *c5 = new TCanvas("c5", "contour", 200, 10, 600, 400);
    c5->cd();
    double contours[3];
    contours[0] = lamb_max-2.30;
    contours[1] = lamb_max-4.61;
    contours[2] = lamb_max-5.99;
    
    hChi2->SetContour(3,contours);
    hChi2->SetLineWidth(2); // to draw the contours
    hChi2->GetXaxis()->SetRangeUser(-1.2,1.2);
    hChi2->GetYaxis()->SetRangeUser(-1.2,1.2);
    hChi2->Draw("cont1");
    
    // Draw a marker
    Double_t BM = gPad->GetBottomMargin();
    Double_t LM = gPad->GetLeftMargin();
    Double_t RM = gPad->GetRightMargin();
    Double_t TM = gPad->GetTopMargin();
    Double_t X1 = hChi2->GetXaxis()->GetXmin()-BinChi2; //---> needed to resize the pad and compute assign proper marker position
    Double_t Y1 = hChi2->GetYaxis()->GetXmin()-BinChi2;
    Double_t X2 = hChi2->GetXaxis()->GetXmax()-BinChi2;
    Double_t Y2 = hChi2->GetYaxis()->GetXmax()-BinChi2;
    
    TPad *null=new TPad("null","null",0,0,1,1);
    
    null->SetFillStyle(0);
    null->SetFrameFillStyle(0);
    null->Draw();
    null->cd();
    null->Range(X1-(X2-X1)*(LM/(1-RM-LM)),
		Y1-(Y2-Y1)*(BM/(1-TM-LM)),
		X2+(X2-X1)*(RM/(1-RM-LM)),
		Y2+(Y2-Y1)*(TM/(1-TM-LM)));
    
    gPad->Update();
    TMarker *mark = new TMarker(0,0,8);
    TMarker *mark2 = new TMarker(x_max,y_max,34);
    
    /*
      TCanvas *c4 = new TCanvas("c4", "lambda_map", 200, 10, 600, 400);
 c4-cd();
 
 TColor::InvertPalette();
 hChi2->Draw("colz");
 c4->SaveAs(argv[2]);

 TCanvas *c5 = new TCanvas("c5", "lambda_map", 200, 10, 600, 400);
 c5->cd();
 
double contours[3];
 contours[0] = lamb_min+2.30;
 contours[1] = lamb_min+4.61;
 contours[2] = lamb_min+5.99;

 hcont->DrawCopy("colz");
 hcont->SetContour(3,contours);
 hcont->SetLineWidth(2); // to draw the contours 
 hcont->Draw("cont1 same");
 
 //hChi2->DrawCopy("colz");
 c4->SaveAs(argv[2]);
 //c4->SaveAs("/sps/km3net/users/fbenfe/Moon_shadow/combined_shadow/data/true_azi/arca19_moon_6deg_-8h.pdf");
 */

    std::string hist_name = std::string(argv[2])+std::string("_bin_")+std::string(argv[4])+std::string("_2d_histogram.pdf");
    const char* hist_name_char = hist_name.c_str();
    std::string contour_name = std::string(argv[2])+std::string("_bin_")+std::string(argv[4])+std::string("_contour.pdf");
    const char* contour_name_char = contour_name.c_str();

    c4->SaveAs(hist_name_char);
    c5->SaveAs(contour_name_char);


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<minutes>(stop - start);
    cout << "-------------------------------------"<<endl;
    cout << "Time elapsed: "
	 << duration.count() << " minutes" << endl;
    cout << "-------------------------------------"<<endl;
    
    //app.Run(kTRUE); 
    cout<<"finished"<<endl;
    return 0;
}

