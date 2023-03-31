#include "TVirtualFitter.h"
#include "TFile.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TColor.h"
 
using namespace std;

Double_t bkg(Double_t *x, Double_t *par) {
  return par[0];
}

Double_t gaus1(Double_t *x, Double_t *par) {
  Double_t R = 0.26;
  Double_t shad = TMath::Pi()*R*R;
  return par[0]*(1-(par[2]*shad/(2*TMath::Pi()*par[1]*par[1])*TMath::Exp(-(x[0]*x[0])/(2*par[1]*par[1]))));
  //return par[0]*(1-(R*R/(2*par[1]*par[1])*TMath::Exp(-0.5*(x[0]*x[0])/(2*par[1]*par[1]))));
}

Double_t gaus_fit(Double_t x, Double_t k, Double_t shad, Double_t sig) {
  Double_t R = 0.26;
    return k*(1-(shad*R*R/(2*sig*sig)*TMath::Exp(-0.5*(x*x)/(2*sig*sig))));
}


int main(int argc, char** argv) {

  if (argc < 3) {
    cout << "3 arguments are required: <input_file> <output_plots_name> <beta_cut> "<<endl;
    exit(1);
  } 
  
 //TApplication app("app",NULL,NULL); 
 TH1F* h[4];
 TF1* f[3];
 TCanvas* c[3];
  
 //gROOT->SetStyle("Plain");
 //gStyle->SetOptFit(1111);
 //gStyle->SetOptStat("nemr");
 gStyle->SetOptStat(0);

 //double MinGrad = 0.1, MaxGrad = 4.1;
 //double bin_width = 0.1;
 //int Nbin = (MaxGrad - MinGrad ) / bin_width;
 double MinGrad = 0.1, MaxGrad = 4.0;
 double Nbin = 39, bin_width = 0.1;
 
 string daz, dalt, fitinf0, sundist, MC_sun_dist, lik, len, time, run;
 double r=0;
 double densities[40], errors[40];

 double beta_cut = std::stof(argv[3]);//0.32;//5.;//0.4 for arca8 moon, 0.34 for arca8 sun;
 
 float size = 0.05; // label and title size for various axes

 cout<<"Analyzing file "<<argv[1]<<endl;
 cout<<"Saving plot in "<<argv[2]<<endl;
 cout<<"Chosen Beta0 cut = "<<argv[3]<<endl;
 
 
 // ---------- define hist of 1D event density distribution

 h[1] = new TH1F("h1","",Nbin,MinGrad,MaxGrad);
 //h[1] = new TH1F("h1","",40,0,4);
 h[1]->SetXTitle("Distance [deg]");
 h[1]->SetYTitle("Event density [1/deg^2] ");
 h[1]->SetMarkerStyle(kOpenCircle);
 h[1]->SetMarkerSize(0.5);
 h[1]->SetMarkerColor(kBlue);
 h[1]->SetLineColor(kBlue);
 h[1]->SetFillColor(0);
 h[1]->GetXaxis()->SetTitleSize(size);
 h[1]->GetXaxis()->SetTitleOffset(1);
 h[1]->GetXaxis()->SetTitleFont(62);
 h[1]->GetXaxis()->SetLabelFont(62);
 h[1]->GetXaxis()->SetLabelSize(size);
 h[1]->GetYaxis()->SetTitleSize(size);
 h[1]->GetYaxis()->SetTitleOffset(1);
 h[1]->GetYaxis()->SetTitleFont(62);
 h[1]->GetYaxis()->SetLabelFont(62);
 h[1]->GetYaxis()->SetLabelSize(size);
 //h[1]->GetYaxis()->SetRangeUser(0.,200);
 h[1]->GetXaxis()->SetRangeUser(0.,4);
 
 
  // ---------- loop over ntuple muon events --------------------------
   ifstream fstream;
   int ibin=0,i=0;
   double area=0;
   double ntot=0;
   for (r=0;r<MaxGrad;r=r+bin_width){
     double R=r+bin_width; 
     double n1=0,n1dens=0;
     //int bin1 = h[1]->GetBin(ibin);
     //fstream.open("/sps/km3net/users/fbenfe/Moon_shadow/csv/arca19/data/arca19_data_pre+anoise+sun_fake_12.txt");
     fstream.open(argv[1]);

 {
      while ( ! fstream.eof() )
	{//aware of n. of columns: remove time and run_id if not present in data !!!!
	  fstream >> daz >> dalt >> fitinf0 >> sundist >> lik >> len >> time >> run;
     double DeltaAzi = std::stof(daz);
     double DeltaAlt = std::stof(dalt);
     double beta0 = std::stof(fitinf0)*TMath::RadToDeg();
     double sun_dist = std::stof(sundist);
     double likelihood = std::stof(lik);
     double track_length = std::stof(len);
     int time_since_epoch = std::stoi(time);
     int run_id = std::stoi(run);
    
     if ( beta0 < beta_cut ) // <------------------------------------------------------ CUT ----------------
       { 
	 if (sun_dist < R && sun_dist > r) {
	   //if (track_length > 200 && likelihood > 60) {
	   //if (track_length > 200) {
	   n1++;
	   }
	 //}
     }
    }
      area = (TMath::Pi()*(R*R-r*r));
      n1dens=n1/area;
      h[1]->SetBinContent(ibin,n1dens);
      h[1]->SetBinError(ibin,TMath::Sqrt(n1)/area);
      //cout<<"bin center = "<< h[1]->GetBinCenter(ibin)<<endl;
      //cout<<"bin low edge = "<< h[1]->GetBinLowEdge(ibin)<<endl;
      densities[i] = n1dens;
      errors[i] = TMath::Sqrt(n1)/area;
      ibin++;
      cout<<"r = "<<r<<endl;
      cout<<"n = "<<n1<<endl;
      cout<<"area = "<<area<<endl;
      cout<<"density = "<<n1dens<<endl;
      i++;
      ntot+=n1;
 }
   // end event loop
 fstream.close();
   } 

   //set y axis limits
   double min = densities[0]-errors[0],max = densities[0]+errors[0];
   for(i = 1;i < 40; i++) {
     // Change < to > if you want to find the smallest element
     if(densities[i]-errors[i] < min) {
       min = densities[i]-errors[i];
     }
     if(densities[i]+errors[i] > max) {
       max = densities[i]+errors[i];
     }
   }

   h[1]->GetYaxis()->SetRangeUser(min-50,max+50);
   //h[1]->GetYaxis()->SetRangeUser(0,max+50);

 
 // fit histograms
   
   c[0] = new TCanvas("c[0]", "c[0]", 10, 20, 600, 400);
   
   
   //fit gaus on S1
   f[0] = new TF1("f[0]", gaus1,MinGrad-0.1,MaxGrad,3);
   f[0]->SetParameter(0,100);
   f[0]->SetParameter(1,0.29);
   //f[0]->FixParameter(2,1.);
   f[0]->SetParameter(2,0.5);
   f[0]->SetParLimits(1,0.,1.);
   f[0]->SetParLimits(2,0.,10.);
   //   f[0]->SetParLimits(0,0.,800.);
   f[0]->SetParNames ("k","sigma","amp");
   f[0]->SetLineColor(kRed);
   f[0]->SetLineStyle(1);
   h[1]->Fit("f[0]","");
   
   /*
   Double_t amin_gaus,amin_bkg,amin_s0,chi2Fitgaus,chi2Fitbkg,chi2_s0;
   Double_t edm,errdef;
   Int_t nvpar,nparx;

   
   TVirtualFitter *fitter_gaus = TVirtualFitter::GetFitter();

   cout<< "TVirtualFitter Gaus S1 stats: " <<endl;
   fitter_gaus->GetStats(amin_gaus,edm,errdef,nvpar,nparx);
   chi2Fitgaus = amin_gaus*2;
   cout << "chi2Fitgaus = "<< chi2Fitgaus<<endl;
   */

   double chi2Fitgaus = f[0]->GetChisquare();
   //cout<<"classical chi = "<<chi2Fitgaus<<endl;
   
   //fit bkg on S1
   //f[1] = new TF1("f[1]",bkg,MinGrad-0.5,MaxGrad,1);
   f[1] = new TF1("f[1]",bkg,MinGrad-0.1,MaxGrad,1);
   f[1]->SetParameter(0,f[0]->GetParameter(0));  //<----- k is the constant on gaus1 fit!
   f[1]->SetParNames("k");
   f[1]->SetLineColor(kRed);
   f[1]->SetLineStyle(2);
   h[1]->Fit("f[1]","");

   /*
   TVirtualFitter *fitter_bkg = TVirtualFitter::GetFitter();
   cout<< "TVirtualFitter bkg S1 stats: " <<endl;
   fitter_bkg->GetStats(amin_bkg,edm,errdef,nvpar,nparx);
   chi2Fitbkg = amin_bkg*2;
   cout << "chi2bkg = "<< chi2Fitbkg<<endl;
   */
   double chi2Fitbkg = f[1]->GetChisquare();
   double lambda = chi2Fitgaus - chi2Fitbkg;

   
   cout<<"-----------------------------S1 fit------------------------------"<<endl;
   cout<<"Total events: "<<ntot<<endl;
   cout<<"Gaussian Fit : " <<'\n'
       <<"k: "<<f[0]->GetParameter(0)<<"+-"<<f[0]->GetParError(0)<<'\n';
   cout<<"sigma: "<<f[0]->GetParameter(1)<<"+-"<<f[0]->GetParError(1)<<'\n';
   cout<<"amp: "<<f[0]->GetParameter(2)<<"+-"<<f[0]->GetParError(2)<<'\n';
   cout<<"Chi square: "<<chi2Fitgaus<<'\n';
   cout<<"NDF: "<<f[0]->GetNDF()<<'\n';
   cout<<"Prob.H1: "<<f[0]->GetProb()<<"\n";
    
   cout<<"BKG fit : " <<'\n'
       <<"k: "<<f[1]->GetParameter(0)<<"+-"<<f[0]->GetParError(0)<<'\n';
   cout<<"Chi square: "<<chi2Fitbkg<<'\n';
   cout<<"NDF: "<<f[1]->GetNDF()<<'\n';
   cout<<"Prob.H0: "<<f[1]->GetProb()<<"\n";
   cout<<"Delta Chi2 (H1-H0) = "<<lambda<<endl;
   cout <<"p-value(0,0) = " << TMath::Prob(-lambda,2) <<endl;
   cout <<"sigma = "<<sqrt(2)*TMath::ErfInverse(1.-TMath::Prob(-lambda,2)) << endl; // 2 becomes 1 if for nominal fit sigma is fixed
    
   c[0]->cd();
   h[1]->Draw("E,P,SAME");
   f[0]->Draw("SAME");
   f[1]->Draw("SAME");
   c[0]->Draw();

   
   TLegend *legend = new TLegend(0.65,0.2,0.85,0.4);
   legend->AddEntry("f[0]","Gaussian fit","l");
   legend->AddEntry("f[1]","Background fit","l");
   legend->Draw();

   c[0]->SaveAs(argv[2]);
   //gStyle->SetOptFit(1111);
   

   
   //app.Run(kTRUE);
 cout<<"finished"<<endl;
 return 0;
   }
