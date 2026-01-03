#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include <vector>
#include <iostream>

void paper_data(){

    TFile *f=new TFile("/data03/tianye/pbar/root/paper_data/paper_data.root","recreate");


    // Define the normal distribution function
    TF1 *normal = new TF1("normal", "TMath::Gaus(x, 1, 0.1)", 0, 2);

    // Define the exponential distribution function
    TF1 *exponential = new TF1("exponential", "[0]*TMath::Exp(-x/[1])", 0, 2);
    exponential->SetParameters(1, 2); // Set initial parameters for the exponential function

    std::vector<int>nmc= {100};
    for(int k=0;k<5;k++){
    gRandom->SetSeed(k); // Set the random seed to 0 for reproducibility
    for(int i=0;i<nmc.size();i++){
        for(int j=0;j<2000;j++){
            TH1D* sig_hist= new TH1D(Form("sig_hist_seed%d_nmc%d_%d",k,nmc[i],j), Form("sig_hist_seed%d_nmc%d_%d",k,nmc[i],j), 15, 0, 2);
            TH1D* bkg_hist= new TH1D(Form("bkg_hist_seed%d_nmc%d_%d",k,nmc[i],j), Form("bkg_hist_seed%d_nmc%d_%d",k,nmc[i],j), 15, 0, 2);
            TH1D* data_hist= new TH1D(Form("data_hist_seed%d_nmc%d_%d",k,nmc[i],j), Form("data_hist_seed%d_nmc%d_%d",k,nmc[i],j), 15, 0, 2);
            gRandom->SetSeed(0); // Set the random seed to 0 for reproducibility
            sig_hist->FillRandom("normal", nmc[i]);
            bkg_hist->FillRandom("exponential", nmc[i]);
            data_hist->FillRandom("normal", 250);
            data_hist->FillRandom("exponential", 750);
            f->cd();
            sig_hist->Write();
            bkg_hist->Write();
            data_hist->Write();
            delete sig_hist;
            delete bkg_hist;
            delete data_hist;
        }
        
  }
  cout<<"k="<<k<<endl;
    }
    f->Close();
}