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

void RooFit_binned(int seq, int npi, int low, int high){
       TStopwatch timer;
    // 启动计时器
    timer.Start();
    //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    //TFile *f1 = new TFile("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq0001_npi5000.root","read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc_crys/binned_IWLS1_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
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

        TH1D *data_pi25=(TH1D*)f->Get(Form("data_pi25_%04d",j));
        TH1D *data_pi100=(TH1D*)f->Get(Form("data_pi100_%04d",j));
        TH1D *data_pr25=(TH1D*)f->Get(Form("data_pr25_%04d",j));
        TH1D *data_pr100=(TH1D*)f->Get(Form("data_pr100_%04d",j));
        TH1D *data_neg25=(TH1D*)f->Get(Form("data_neg25_%04d",j));
        TH1D *data_neg100=(TH1D*)f->Get(Form("data_neg100_%04d",j));


        // for(int iter=0;iter<10;iter++){
        //     IWLS_HistFit(data_pr25,data_pi25,data_neg25,neg_error_bin25,nsig[0]);
        //     IWLS_HistFit(data_pr100,data_pi100,data_neg100,neg_error_bin25,nsig[1]);
        // }

        IWLS_CrysFunctionFit_Minuit2(data_pr25,data_pi25,data_neg25,nsig[0],nsig[2],5);
        IWLS_CrysFunctionFit_Minuit2(data_pr100,data_pi100,data_neg100,nsig[1],nsig[3],5);

    



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
    //delete c1;
        std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/Root_Func_binned_bin25_IWLS_crys_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "Root_Func_binned_bin25_IWLS_crys - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "Root_Func_binned_bin25_IWLS_crys - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
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