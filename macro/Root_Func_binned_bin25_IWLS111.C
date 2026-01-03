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



void RooFit_binned(int seq, int npi, int low, int high){
        // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    // TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_multi_bin/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_morebin/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    //TFile *f1 = new TFile("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq0001_npi5000.root","read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_IWLS1_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
    std::cout<<"----open file----"<<std::endl;
    std::cout<<seq<<" "<<npi<<" "<<low<<" "<<high<<std::endl;
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    // std::vector<double> params1(12), params2(12), params3(12), params4(12);
    // fitParamsTree->Branch("bin20_chi2_pararr", &params1);
    // fitParamsTree->Branch("bin100_chi2_pararr", &params2);
    // fitParamsTree->Branch("bin20_mll_pararr", &params3);
    // fitParamsTree->Branch("bin100_mll_pararr", &params4);

    std::vector<double> nsig(4);
    fitParamsTree->Branch("nsig", &nsig);
    
    // TH1D *data_pr25=(TH1D*)f1->Get("data_pr25_9999");
    // TH1D *data_pr100=(TH1D*)f1->Get("data_pr100_9999");
    // TH1D *data_neg25=(TH1D*)f1->Get("data_neg25_9999");
    // TH1D *data_neg100=(TH1D*)f1->Get("data_neg100_9999");



    // TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    // c1->Divide(2, 2);
    // c1->Print(Form("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/fit_plots_seq%04d_npi%04d_%04d.pdf[", seq, npi, low));


    

    for(int j=low;j<high;j++){

        TH1D *data_pi25=(TH1D*)f->Get(Form("data_pi30_%04d",j));
        TH1D *data_pi100=(TH1D*)f->Get(Form("data_pi100_%04d",j));
        TH1D *data_pr25=(TH1D*)f->Get(Form("data_pr30_%04d",j));
        TH1D *data_pr100=(TH1D*)f->Get(Form("data_pr100_%04d",j));
        TH1D *data_neg25=(TH1D*)f->Get(Form("data_neg30_%04d",j));
        TH1D *data_neg100=(TH1D*)f->Get(Form("data_neg100_%04d",j));


        // for(int iter=0;iter<10;iter++){
        //     IWLS_HistFit(data_pr25,data_pi25,data_neg25,neg_error_bin25,nsig[0]);
        //     IWLS_HistFit(data_pr100,data_pi100,data_neg100,neg_error_bin25,nsig[1]);
        // }

        IWLS_FunctionFit_Minuit2(data_pr25,data_pi25,data_neg25,nsig[0],nsig[2],5);
        IWLS_FunctionFit_Minuit2(data_pr100,data_pi100,data_neg100,nsig[1],nsig[3],5);

    



        if(!nsig.empty()){
            //std::cout<<"----fill----"<<std::endl;
            fitParamsTree->Fill();
        }
        else{
            std::cout<<"----empty----"<<std::endl;
            continue;
        }

    
        // if(j%100==0){
        //     TF1* fit_pi1=new TF1("fit_pi1",funcExpGausExp,-1,3,7);
        //     TF1* fit_neg1=new TF1("fit_neg1",combine_chi2,-1,3,15);
        //     TF1* fit_pi2=new TF1("fit_pi2",funcExpGausExp,-1,3,7);
        //     TF1* fit_neg2=new TF1("fit_neg2",combine_chi2,-1,3,15);

        //     Double_t params1_pi1[] = {params1[5], params1[6], params1[7], params1[8], params1[9], -1, 3};
        //     fit_pi1->SetParameters(params1_pi1);
        //     Double_t params1_neg1[] = {params1[0], params1[1], params1[2], params1[3], params1[4], -1, 3, params1[5], params1[6], params1[7], params1[8], params1[9], -1, 3, params1[10]/data_neg25->GetEntries()};
        //     fit_neg1->SetParameters(params1_neg1);
        //     Double_t params2_pi2[] = {params2[5], params2[6], params2[7], params2[8], params2[9], -1, 3};
        //     fit_pi2->SetParameters(params2_pi2);
        //     Double_t params2_neg2[] = {params2[0], params2[1], params2[2], params2[3], params2[4], -1, 3, params2[5], params2[6], params2[7], params2[8], params2[9], -1, 3, params2[10]/data_neg100->GetEntries()};
        //     fit_neg2->SetParameters(params2_neg2);

        //     TCanvas *c1 = new TCanvas(Form("j%04d",j), Form("j%04d",j), 800, 600);
        //     c1->Divide(2, 2);

        //     c1->cd(1);
        //     data_pi25->Scale(1.0 / data_pi25->GetSumOfWeights() / data_pi25->GetBinWidth(1));
        //     data_pi25->Draw();
        //     fit_pi1->Draw("same");

        //     c1->cd(2);
        //     data_pi100->Scale(1.0 / data_pi100->GetSumOfWeights() / data_pi100->GetBinWidth(1));
        //     data_pi100->Draw();
        //     fit_pi2->Draw("same");

        //     c1->cd(3);
        //     data_neg25->Scale(1.0 / data_neg25->GetSumOfWeights() / data_neg25->GetBinWidth(1));
        //     data_neg25->Draw();
        //     fit_neg1->Draw("same");

        //     c1->cd(4);
        //     data_neg100->Scale(1.0 / data_neg100->GetSumOfWeights() / data_neg100->GetBinWidth(1));
        //     data_neg100->Draw();
        //     fit_neg2->Draw("same");


        //     c1->Print(Form("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/fit_plots_seq%04d_npi%04d_%04d.pdf", seq, npi, low));

        //     delete fit_pi1;
        //     delete fit_neg1;
        //     delete fit_pi2;
        //     delete fit_neg2;
        // }



        // params1.clear();         
        // params1.shrink_to_fit();
        // params2.clear();         
        // params2.shrink_to_fit();
        // params3.clear();         
        // params3.shrink_to_fit(); 
        // params4.clear();         
        // params4.shrink_to_fit(); 
        

        data_pi25->Delete();
        data_pi100->Delete();
        data_pr25->Delete();
        data_pr100->Delete();
        data_neg25->Delete();
        data_neg100->Delete();


        std::cout<<"j = "<<j<<std::endl;
        std::cout<<"nsig[0] = "<<nsig[0]<<std::endl;
        std::cout<<"nsig[1] = "<<nsig[1]<<std::endl;
        std::cout<<"nsig[2] = "<<nsig[2]<<std::endl;
        std::cout<<"nsig[3] = "<<nsig[3]<<std::endl;
    }
    //c1->Print(Form("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/fit_plots_seq%04d_npi%04d_%04d.pdf]", seq, npi, low));
    // data_pr25->Delete();
    // data_pr100->Delete();
    // data_neg25->Delete();
    // data_neg100->Delete();
    fOutput->cd();
    fitParamsTree->Write();
    fitParamsTree->Delete();
    fOutput->Close();
    std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/Root_Func_binned_bin25_IWLS_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "Root_Func_binned_bin25_IWLS - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "Root_Func_binned_bin25_IWLS - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
        timeFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
   
}

int main(int argc,char* argv[]){
    RooFit_binned(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    std::cout<<"----done----"<<std::endl;
    return 0;
}