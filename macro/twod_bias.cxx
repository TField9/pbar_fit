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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>

void twod_bias(){
    // 配置
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,};
    const char* names[] = {
        "RooComb","TFF","IWLS","DA","JSC","asy","URComb","URFit", // 

    };

    int n =8;
    TH1D *h_fit[n];
    TF1 *f_fit[n];
    TH2D *h_bias[n], *h_sigma[n];




    for(int i=0;i<n;++i){
        h_bias[i]  = new TH2D(Form("h_%s_bias", names[i]), Form("h_%s_bias;npr;npi", names[i]), 9, 0, 9,9,0,9);
        h_sigma[i] = new TH2D(Form("h_%s_sigma", names[i]), Form("h_%s_sigma;npr;npi", names[i]), 9, 0, 9,9,0,9);
    }

    TFile *f =new TFile("/data03/tianye/pbar/root/bias_2d/2d_bias.root","recreate");


    for(int npr=0;npr<9;npr++){
        for(int npi=0;npi<9;npi++){

            for(int i=0;i<n;++i){
                h_fit[i]=new TH1D(Form("h_%s_fit_%d_%d", names[i],size[npr],size[npi]), Form("h_%s_fit_%d_%d", names[i],npr,npi), 300, 0, 300);
                f_fit[i] = new TF1(Form("f_%s_fit_%d_%d", names[i],size[npr],size[npi]), "gaus", 0, 300);
            }


            TFile *f1=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_roocomb_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f2=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_TFF_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f3=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_IWLS_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f4=new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_OM_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f5=new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_roocomb_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");
            TFile *f6=new TFile(Form("/data03/tianye/pbar/root/2d_result/unbinned_roofit_seq0001_npr%04d_npi%04d.root", size[npr], size[npi]),"read");

            

                    // 取树
            TTree *t_roocomb = (TTree*)f1->Get("fitParamsTree");
            TTree *t_TFF     = (TTree*)f2    ->Get("fitParamsTree");
            TTree *t_IWLS    = (TTree*)f3   ->Get("fitParamsTree");
            TTree *t_OM      = (TTree*)f4    ->Get("fitParamsTree");
            TTree *t_unRC    = (TTree*)f5   ->Get("fitParamsTree");
            TTree *t_unRF    = (TTree*)f6   ->Get("fitParamsTree");

                        // 参数容器
            std::vector<double>* p_roocomb = new std::vector<double>(2);
            std::vector<double>* p_TFF     = new std::vector<double>(2);
            std::vector<double>* p_IWLS    = new std::vector<double>(2);
            Double_t pDA[2]={0}, pJSC[2]={0}, pasy[2]={0};
            std::vector<double>* p_unRC    = new std::vector<double>(2);
            std::vector<double>* p_unRF    = new std::vector<double>(2);


            // 绑定
            t_roocomb->SetBranchAddress("nsig", &p_roocomb);
            t_TFF    ->SetBranchAddress("nsig", &p_TFF);
            t_IWLS   ->SetBranchAddress("nsig", &p_IWLS);
            t_OM->SetBranchAddress("Ysig_DA",   pDA);
            t_OM->SetBranchAddress("Ysig_JSC",  pJSC);
            t_OM->SetBranchAddress("Ysig_asy",  pasy);
            t_unRC   ->SetBranchAddress("nsig", &p_unRC);
            t_unRF   ->SetBranchAddress("nsig", &p_unRF);


            //////////must reset//////////

            for(int i = 0; i < n; ++i){
                h_fit[i]->Reset();
            }


            for(int j=0;j<10000;++j){
                t_roocomb->GetEntry(j);
                t_TFF    ->GetEntry(j);
                t_IWLS   ->GetEntry(j);
                t_OM     ->GetEntry(j);
                t_unRC   ->GetEntry(j);
                t_unRF   ->GetEntry(j);

                // 填充
                h_fit[0]->Fill((*p_roocomb)[0]);
                h_fit[1]->Fill((*p_TFF)[0]);
                h_fit[2]->Fill((*p_IWLS)[0]);
                h_fit[3]->Fill(pDA[0]);
                h_fit[4]->Fill(pJSC[0]);
                h_fit[5]->Fill(pasy[0]);
                h_fit[6]->Fill((*p_unRC)[0]);
                h_fit[7]->Fill((*p_unRF)[0]);

            }

            for(int i = 0; i < n; ++i) {
                f_fit[i]->SetParameters(h_fit[i]->GetMean(), h_fit[i]->GetRMS());
                h_fit[i]->Fit(f_fit[i], "Q0");
                h_bias[i]->SetBinContent(npr+1, npi+1, f_fit[i]->GetParameter(1)/100);
                h_bias[i]->SetBinError(npr+1, npi+1, f_fit[i]->GetParError(1)/100);
                h_sigma[i]->SetBinContent(npr+1, npi+1, f_fit[i]->GetParameter(2)/100);
                h_sigma[i]->SetBinError(npr+1,npi+1,f_fit[i]->GetParError(2)/100);
            }
            f->cd();
            
            
            for(int i=0;i<n;++i){
                h_fit[i]->Write();
                f_fit[i]->Write();
            }




        }
        cout<<"npr="<<npr<<endl;
        
    }
    f->cd();
    for(int i=0;i<n;i++){
        h_bias[i]->Write();
        h_sigma[i]->Write();
    }
    f->Close();


}