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

void oned_bias_check(){
    // 配置 - 只保留roocomb方法
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,200,400};
    int nsize = size.size();
    const char* names = "RooComb";  // 只保留roocomb方法
    
    int met = 1;  // 只需要一个方法
    TH1D* h_fit;
    TF1 *f_fit;
    TGraphErrors *g_bias[nsize], *g_sigma[nsize]; // 维度调整为nsize

    // 初始化图
    for(int j=0; j<nsize; j++){
        TString bias = Form("g_%s_bias_Signal_%d", names, size[j]);
        TString sigma = Form("g_%s_sigma_Signal_%d", names, size[j]);
        
        g_bias[j] = new TGraphErrors();
        g_bias[j]->SetName(bias);
        g_bias[j]->SetTitle(Form("%s;N_{B};Bias(%%)", bias.Data()));
        
        g_sigma[j] = new TGraphErrors();
        g_sigma[j]->SetName(sigma);
        g_sigma[j]->SetTitle(Form("%s;N_{B};Sigma", sigma.Data()));
    }

    TFile *f = new TFile("/data03/tianye/pbar/root/bias_2d/1d_bias_check.root","recreate");

    for(int npr=0; npr<nsize; npr++){
        for(int npi=0; npi<nsize; npi++){
            // 只创建roocomb需要的对象
            h_fit = new TH1D(Form("h_%s_fit_%d_%d", names, size[npr], size[npi]), 
                            Form("h_%s_fit_%d_%d", names, npr, npi), 300, 0, 300);
            f_fit = new TF1(Form("f_%s_fit_%d_%d", names, size[npr], size[npi]), "gaus", 0, 300);
            
            // 只打开roocomb文件
            TFile *f_roocomb = new TFile(Form("/data03/tianye/pbar/root/2d_result/binned_roocomb_npr%04d_npi%04d.root", 
                                             size[npr], size[npi]), "read");
            
            // 只获取roocomb树
            TTree *t_roocomb = (TTree*)f_roocomb->Get("fitParamsTree");
            std::vector<double>* p_roocomb = new std::vector<double>(2);
            t_roocomb->SetBranchAddress("nsig", &p_roocomb);
            
            // 重置直方图
            h_fit->Reset();
            
            // 填充直方图
            for(int j=0; j<10000; ++j){
                t_roocomb->GetEntry(j);
                h_fit->Fill((*p_roocomb)[0]);
            }
            
            // 拟合并设置图表
            f_fit->SetParameters(h_fit->GetMean(), h_fit->GetRMS());
            h_fit->Fit(f_fit, "Q0");
            g_bias[npr]->SetPoint(npi, size[npi], (f_fit->GetParameter(1)/100 - 1)*100);
            g_bias[npr]->SetPointError(npi, size[npi], f_fit->GetParError(1)/100);
            g_sigma[npr]->SetPoint(npi, size[npi], f_fit->GetParameter(2)/100);
            g_sigma[npr]->SetPointError(npi, size[npi], f_fit->GetParError(2)/100);
            
            // 写入并清理
            f->cd();
            h_fit->Write();
            f_fit->Write();
            f_roocomb->Close();
            delete f_roocomb;
        }
        cout << "npr=" << npr << endl;
    }
    
    // 写入最终结果
    f->cd();
    for(int j=0; j<nsize; j++){
        g_bias[j]->Write();
        g_sigma[j]->Write();
    }
    f->Close();
}