#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TColor.h"
 
using namespace std;

int scan_beta_plot() {
 
  TCanvas* c[3];
  TGraph* gr[3];

  gr[0] = new TGraph();
  gr[0]->SetTitle("Beta scan; beta_cut (deg) ; lambda");
  gr[0]->SetLineColor(46);
  gr[0]->SetLineWidth(4);
  gr[0]->SetMarkerColor(4);
  gr[0]->SetMarkerSize(1);
  gr[0]->SetMarkerStyle(23);

  gr[1] = new TGraph();
  gr[1]->SetTitle("Beta scan; beta_cut (deg) ; sigma");
  gr[1]->SetLineColor(46);
  gr[1]->SetLineWidth(4);
  gr[1]->SetMarkerColor(4);
  gr[1]->SetMarkerSize(1);
  gr[1]->SetMarkerStyle(23);
  

  gr[2] = new TGraph();
  gr[2]->SetTitle("Beta scan; beta_cut (deg) ; amplitude");
  gr[2]->SetLineColor(46);
  gr[2]->SetLineWidth(4);
  gr[2]->SetMarkerColor(4);
  gr[2]->SetMarkerSize(1);
  gr[2]->SetMarkerStyle(23);
 
 string beta0, lamb, sig, amp;
 double i=0;

 float size = 0.05; // label and title size for various axes


 
 ifstream fstream;
 fstream.open("csv/arca8/arca8rbr_mu_pre+anoise+sun_betas_5.txt");
 {
   while (!fstream.eof()) {
	
     fstream >> beta0 >> lamb >> sig >> amp;
   
     double beta_cut = std::stof(beta0);
     double lambda = std::stof(lamb);
     double sigma = std::stof(sig);
     double a = std::stof(amp);		      
     gr[0]->SetPoint(i, beta_cut, lambda);
     gr[1]->SetPoint(i, beta_cut, sigma);
     gr[2]->SetPoint(i, beta_cut, a);
     i++;
   }
 }
 fstream.close();
      
 c[0] = new TCanvas("c[0]", "c[0]", 10, 20, 600, 400);
 c[1] = new TCanvas("c[1]", "c[1]", 10, 20, 600, 400);
 c[2] = new TCanvas("c[2]", "c[2]", 10, 20, 600, 400);
 
 c[0]->cd();
 c[0]->SetGrid();
 gr[0]->Draw("ACP");
 c[0]->Draw();
 
 c[1]->cd();
 c[1]->SetGrid();
 gr[1]->Draw("ACP");
 c[1]->Draw();
 
 c[2]->cd();
 c[2]->SetGrid();
 gr[2]->Draw("ACP");
 c[2]->Draw();


 cout<<"finished"<<endl;
 return 0;
 }
