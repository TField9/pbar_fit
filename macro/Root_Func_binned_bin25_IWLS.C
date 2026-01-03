#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"
#include <cstdlib>
#include <vector>
#include <string>
#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TH3F.h"
#include "TMath.h"
#include "RooCrystalBall.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include <iostream>
#include <map>
#include <stdio.h>
#include "advfit_unbinned.h"
#include "Math/MinimizerOptions.h"
#include "RooMsgService.h"
#include "RooFit.h"
#include <fstream>
#include <chrono>
#include <thread>

// void setParameterWithDebug(TF1* func, int parIndex, double value, double low, double up) {
//     if (value < low || value > up) {
//         std::cerr << func->GetName()<<" "<<"Warning: Parameter " << parIndex << " value " << value << " is outside the range [" << low << ", " << up << "]. Setting to (" << low << " + " << up << ") / 2." << std::endl;
//         value = (low + up) / 2;
//     }
//     func->SetParameter(parIndex, value);
//     func->SetParLimits(parIndex, low, up);
// }


std::vector<double> Root_chi2(TH1D* hpr,TH1D* hpi, TH1D* hneg){

    TH1D* hpr_clone = (TH1D*)hpr->Clone();
    hpr_clone->Scale(1.0 / hpr->Integral() / hpr->GetBinWidth(1));

    TH1D* hpi_clone = (TH1D*)hpi->Clone();
    hpi_clone->Scale(1.0 / hpi->Integral() / hpi->GetBinWidth(1));

    TH1D* hneg_clone = (TH1D*)hneg->Clone();
    hneg_clone->Scale(1.0 / hneg->Integral() / hneg->GetBinWidth(1));

    TF1* fit1=new TF1("fpr",funcExpGausExp,-1,3,7);
    TF1* fit2=new TF1("fpi",funcExpGausExp,-1,3,7);
    TF1* fit3=new TF1("fneg",combine_chi2,-1,3,15);

    setParameterWithDebug(fit1, 0, hpr_clone->GetMean(), 0.5, 1.1);
    setParameterWithDebug(fit1, 1, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit1, 3, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 4, 1, 0.01, 6.0);
    fit1->FixParameter(5, -1);
    fit1->FixParameter(6, 3);

    setParameterWithDebug(fit2, 0, hpi_clone->GetMean(), 0.001, 0.35);
    setParameterWithDebug(fit2, 1, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit2, 3, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 4, 1, 0.01, 5.0);
    fit2->FixParameter(5, -1);
    fit2->FixParameter(6, 3);

    hpr_clone->Fit(fit1,"MNQ", "R");
    hpi_clone->Fit(fit2, "MNQ", "R");


    double comb_params[15];
    for (int i = 0; i < 7; ++i) {
        comb_params[i] = fit1->GetParameter(i);
    }
    for (int i = 0; i < 7; ++i) {
        comb_params[i + 7] = fit2->GetParameter(i);
    }
    comb_params[14] = 0.1;
    fit3->SetParameters(comb_params);
    for (int i = 0; i < 14; ++i) {
      fit3->FixParameter(i, comb_params[i]);
    }

    fit3->SetParLimits(14,0.,1.);
    fit3->SetNpx(150);

    hneg_clone->Fit(fit3, "MNQ", "R");

    std::vector<double> params;
    params.push_back(fit1->GetParameter(0));
    params.push_back(fit1->GetParameter(1));
    params.push_back(fit1->GetParameter(2));
    params.push_back(fit1->GetParameter(3));
    params.push_back(fit1->GetParameter(4));
    params.push_back(fit2->GetParameter(0));
    params.push_back(fit2->GetParameter(1));
    params.push_back(fit2->GetParameter(2));
    params.push_back(fit2->GetParameter(3));
    params.push_back(fit2->GetParameter(4));
    params.push_back(fit3->GetParameter(14) * hneg->GetEntries());
    params.push_back((1 - fit3->GetParameter(14)) * hneg->GetEntries());

    fit1->Delete();
    fit2->Delete();
    fit3->Delete();
    // hpr->Delete();
    // hpi->Delete();
    // hneg->Delete();
    hpr_clone->Delete();
    hpi_clone->Delete();
    hneg_clone->Delete();


    return params;
}

std::vector<double> Root_mll(TH1D* hpr,TH1D* hpi, TH1D* hneg){

    TH1D* hpr_clone = (TH1D*)hpr->Clone();
    hpr_clone->Scale(1.0 / hpr->Integral() / hpr->GetBinWidth(1));

    TH1D* hpi_clone = (TH1D*)hpi->Clone();
    hpi_clone->Scale(1.0 / hpi->Integral() / hpi->GetBinWidth(1));

    TH1D* hneg_clone = (TH1D*)hneg->Clone();
    hneg_clone->Scale(1.0 / hneg->Integral() / hneg->GetBinWidth(1));

    TF1* fit1=new TF1("fpr",funcExpGausExp,-1,3,7);
    TF1* fit2=new TF1("fpi",funcExpGausExp,-1,3,7);
    TF1* fit3=new TF1("fneg",combine_chi2,-1,3,15);

    setParameterWithDebug(fit1, 0, hpr_clone->GetMean(), 0.5, 1.1);
    setParameterWithDebug(fit1, 1, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit1, 3, hpr_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit1, 4, 1, 0.01, 6.0);
    fit1->FixParameter(5, -1);
    fit1->FixParameter(6, 3);

    setParameterWithDebug(fit2, 0, hpi_clone->GetMean(), 0.001, 0.35);
    setParameterWithDebug(fit2, 1, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 2, 1, 0.01, 5.0);
    setParameterWithDebug(fit2, 3, hpi_clone->GetRMS(), 0.05, 1.1);
    setParameterWithDebug(fit2, 4, 1, 0.01, 5.0);
    fit2->FixParameter(5, -1);
    fit2->FixParameter(6, 3);

    hpr_clone->Fit(fit1,"LMNQ", "R");
    hpi_clone->Fit(fit2, "LMNQ", "R");


    double comb_params[15];
    for (int i = 0; i < 7; ++i) {
        comb_params[i] = fit1->GetParameter(i);
    }
    for (int i = 0; i < 7; ++i) {
        comb_params[i + 7] = fit2->GetParameter(i);
    }
    comb_params[14] = 0.1;
    fit3->SetParameters(comb_params);
    for (int i = 0; i < 14; ++i) {
      fit3->FixParameter(i, comb_params[i]);
    }

    fit3->SetParLimits(14,0.,1.);
    fit3->SetNpx(150);

    hneg_clone->Fit(fit3, "LMNQ", "R");

    std::vector<double> params;
    params.push_back(fit1->GetParameter(0));
    params.push_back(fit1->GetParameter(1));
    params.push_back(fit1->GetParameter(2));
    params.push_back(fit1->GetParameter(3));
    params.push_back(fit1->GetParameter(4));
    params.push_back(fit2->GetParameter(0));
    params.push_back(fit2->GetParameter(1));
    params.push_back(fit2->GetParameter(2));
    params.push_back(fit2->GetParameter(3));
    params.push_back(fit2->GetParameter(4));
    params.push_back(fit3->GetParameter(14) * hneg->GetEntries());
    params.push_back((1 - fit3->GetParameter(14)) * hneg->GetEntries());

    fit1->Delete();
    fit2->Delete();
    fit3->Delete();
    // hpr->Delete();
    // hpi->Delete();
    // hneg->Delete();
    hpr_clone->Delete();
    hpi_clone->Delete();
    hneg_clone->Delete();


    return params;
}



void RooFit_binned(int npr,int npi){
        // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    // TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_multi_bin/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npr),"read");
    TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_IWLS_seq0001_npr%04d_npi%04d.root",npr,npi),"recreate");
    std::cout<<"----open file----"<<std::endl;
    std::cout<<"IWLS npr="<<npr<<"npi="<<npi<<std::endl;
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    // std::vector<double> params1(12), params2(12), params3(12), params4(12);
    // fitParamsTree->Branch("bin20_chi2_pararr", &params1);
    // fitParamsTree->Branch("bin100_chi2_pararr", &params2);
    // fitParamsTree->Branch("bin20_mll_pararr", &params3);
    // fitParamsTree->Branch("bin100_mll_pararr", &params4);

    std::vector<double> nsig(2);
    fitParamsTree->Branch("nsig", &nsig);
    

    for(int j=0;j<10000;j++){

        
        TH1D *data_pr30=(TH1D*)f->Get(Form("data_pr30_%04d",j));
        TH1D *data_pi30=(TH1D*)f1->Get(Form("data_pi30_%04d",j));
        TH1D *data_neg30=(TH1D*)f1->Get(Form("data_neg30_%04d",j));
  


        IWLS_HistFit(data_pr30,data_pi30,data_neg30,nsig[0],nsig[1],5);



    

        fitParamsTree->Fill();



    

        data_pi30->Delete();

        data_pr30->Delete();

        data_neg30->Delete();

        if(j%100==0){
        std::cout<<"j = "<<j<<std::endl;
        std::cout<<"nsig[0] = "<<nsig[0]<<std::endl;

        }

    }

    fOutput->cd();
    fitParamsTree->Write();
    fitParamsTree->Delete();
    fOutput->Close();

   
}

int main(int argc,char* argv[]){
    RooFit_binned(atoi(argv[1]),atoi(argv[2]));
    std::cout<<"----done----"<<std::endl;
    return 0;
}