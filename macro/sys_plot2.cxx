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

void sys_plot(){
    std::vector<int> size = {10,15,20,25,30,40,50, 75, 100, 150, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 4000, 5000, 8000};
    int n = size.size();
    int met = 8; // 7种方法 (4+3)
    int nSys = 4; // 两种背景系统atics情况

    // 三类数据 × 九种方法 × n × 两种背景情况
    double sig_yield[nSys][met][n];
    double sig_yield_error[nSys][met][n];
    double fit_sigma[nSys][met][n];
    double fit_sigma_error[nSys][met][n];

    // ROOT 样式
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TFile* bias_graph= new TFile("/data03/tianye/pbar/root/bias_graph/bias_graph.root", "RECREATE");

    // 逐个事例数读取并拟合
    for(int i = 0; i < n; ++i){
        // [1] 确保文件存在
        TString fname;
        
        fname = Form("/data03/tianye/pbar/root/sys_result/binned_roocomb_seq0001_npi%04d.root", size[i]);
        if(gSystem->AccessPathName(fname)) {
            cout << "Missing file: " << fname << endl;
            continue;
        }
        TFile *f_roocomb = new TFile(fname, "read");
        
        fname = Form("/data03/tianye/pbar/root/sys_result/binned_TFF_seq0001_npi%04d.root", size[i]);
        if(gSystem->AccessPathName(fname)) {
            cout << "Missing file: " << fname << endl;
            continue;
        }
        TFile *f_TFF = new TFile(fname, "read");
        
        // 类似检查其他文件...
        //TFile *f_SPD     = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_SPD_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_IWLS    = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_IWLS_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_OM= new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_OM_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_unbinned_roocomb= new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roocomb_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_roofit = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roofit_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_unbinned_rookeyspdf = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_rookeyspdf_diffwidth_seq0001_npi%04d.root", size[i]), "read");
        // [2] 检查Tree是否存在
        TTree *t_roocomb = (TTree*)f_roocomb->Get("fitParamsTree");
        if(!t_roocomb) {
            cout << "Missing tree in f_roocomb" << endl;
            continue;
        }
        
        TTree *t_TFF     = (TTree*)f_TFF->Get("fitParamsTree");
        //TTree *t_SPD     = (TTree*)f_SPD->Get("fitParamsTree");
        TTree *t_IWLS    = (TTree*)f_IWLS->Get("fitParamsTree");
        TTree *t_OM      = (TTree*)f_OM->Get("fitParamsTree");
        TTree *t_un_roocomb = (TTree*)f_unbinned_roocomb->Get("fitParamsTree");
        TTree *t_roofit = (TTree*)f_roofit->Get("fitParamsTree");
        TTree *t_un_rookeyspdf = (TTree*)f_unbinned_rookeyspdf->Get("fitParamsTree");

        // [3] 修复OM方法的读取类型
        // 前四种方法的参数
        std::vector<double>* params_roocomb = new std::vector<double>(8);
        std::vector<double>* params_TFF     = new std::vector<double>(8);
        //std::vector<double>* params_SPD     = new std::vector<double>(8);
        std::vector<double>* params_IWLS    = new std::vector<double>(8);
        std::vector<double>* params_un_roocomb = new std::vector<double>(8);
        std::vector<double>* params_roofit = new std::vector<double>(8);
        std::vector<std::vector<double>>* params_un_rookeyspdf = new std::vector<std::vector<double>> nsig(5, std::vector<double>(8));
        
        // OM方法参数 (使用Double_t数组而不是vector)
        Double_t params_DA[2] = {0};
        Double_t params_JSC[2] = {0};
        Double_t params_asy[2] = {0};
        Double_t params_DA1[2] = {0};
        Double_t params_JSC1[2] = {0};
        Double_t params_asy1[2] = {0};

        t_roocomb->SetBranchAddress("nsig", &params_roocomb);
        t_TFF    ->SetBranchAddress("nsig", &params_TFF);
        //t_SPD    ->SetBranchAddress("nsig", &params_SPD);
        t_IWLS   ->SetBranchAddress("nsig", &params_IWLS);
        t_un_roocomb->SetBranchAddress("nsig", &params_un_roocomb);
        t_roofit->SetBranchAddress("nsig", &params_roofit);
        t_un_rookeyspdf->SetBranchAddress("nsig", &params_un_rookeyspdf);
        t_OM->SetBranchAddress("Ysig_DA", params_DA);     // 直接用数组
        t_OM->SetBranchAddress("Ysig_JSC", params_JSC);   // 直接用数组
        t_OM->SetBranchAddress("Ysig_asy", params_asy);    // 直接用数组
        t_OM->SetBranchAddress("Ysig1_DA", params_DA1);     // 直接用数组
        t_OM->SetBranchAddress("Ysig1_JSC", params_JSC1);   // 直接用数组
        t_OM->SetBranchAddress("Ysig1_asy", params_asy1);    // 直接用数组

        int nEntries = t_roocomb->GetEntries();

        // 对于每种背景系统atics情况分别处理
        for(int j = 0; j < nSys; ++j){
            // 根据系统atics情况调整直方图范围
            double hist_min, hist_max, fit_min, fit_max;
            if(j == 0||j==2) {
                // 第一种背景情况：信号数在140-280范围内
                hist_min = 0;
                hist_max = 400;
                fit_min = 0;
                fit_max = 400;
            } else {
                // 第二种背景情况：信号数在0-100范围内
                hist_min = 0;
                hist_max = 200;
                fit_min = 0;
                fit_max = 200;
            }
            int nBins = 300;
            
            TH1D* h_roocomb = new TH1D(Form("h_roocomb_%d_%d", i, j), Form("h_roocomb_%d_%d;Fitted Signal Event;Count", i, j),nBins, hist_min, hist_max);
            TH1D* h_TFF     = new TH1D(Form("h_TFF_%d_%d", i, j),     Form("h_TFF_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            //TH1D* h_SPD     = new TH1D(Form("h_SPD_%d_%d", i, j),     Form("h_SPD_%d_%d", i, j), nBins, hist_min, hist_max);
            TH1D* h_IWLS    = new TH1D(Form("h_IWLS_%d_%d", i, j),    Form("h_IWLS_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_DA= new TH1D(Form("h_DA_%d_%d", i, j), Form("h_DA_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_JSC= new TH1D(Form("h_JSC_%d_%d", i, j), Form("h_JSC_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_asy= new TH1D(Form("h_asy_%d_%d", i, j), Form("h_asy_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_un_roocomb=new TH1D(Form("h_un_roocomb_%d_%d", i, j), Form("h_un_roocomb_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_roofit = new TH1D(Form("h_roofit_%d_%d", i, j), Form("h_roofit_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);
            TH1D* h_un_rookeyspdf = new TH1D(Form("h_un_rookeyspdf_%d_%d", i, j), Form("h_un_rookeyspdf_%d_%d;Fitted Signal Event;Count", i, j), nBins, hist_min, hist_max);

            for(int k = 0; k < nEntries; ++k){
                t_roocomb->GetEntry(k);
                t_TFF    ->GetEntry(k);
                //t_SPD    ->GetEntry(k);
                t_IWLS   ->GetEntry(k);
                t_OM->GetEntry(k);
                t_un_roocomb->GetEntry(k);
                t_roofit->GetEntry(k);
               
                // 填充直方图，使用对应背景情况的参数
                h_roocomb->Fill((*params_roocomb)[j] );
                h_TFF    ->Fill((*params_TFF)[j] * ((j == 0||j==2 )? 400 : 250));
                //h_SPD    ->Fill((*params_SPD)[j] );
                h_IWLS   ->Fill((*params_IWLS)[j] );
                h_un_roocomb->Fill((*params_un_roocomb)[j] );
                h_roofit->Fill((*params_roofit)[j] );
                
                if(j==0||j==1){
                    // 填充OM方法
                    h_DA->Fill(params_DA[j]);
                    h_JSC->Fill(params_JSC[j]);
                    h_asy->Fill(params_asy[j]);
                }
                else{
                    h_DA->Fill(params_DA1[j-2]);
                    h_JSC->Fill(params_JSC1[j-2]);
                    h_asy->Fill(params_asy1[j-2]);
                }
            }
            bias_graph->cd();
            h_roocomb->Write();
            h_TFF->Write();
            //h_SPD->Write();
            h_IWLS->Write();
            h_DA->Write();
            h_JSC->Write();
            h_asy->Write();
            h_un_roocomb->Write();
            h_roofit->Write();

            TF1* f1 = new TF1("f1", "gaus", fit_min, fit_max);
            TF1* f2 = new TF1("f2", "gaus", fit_min, fit_max);
            TF1* f3 = new TF1("f3", "gaus", fit_min, fit_max);
            TF1* f4 = new TF1("f4", "gaus", fit_min, fit_max);
            TF1* f5 = new TF1("f5", "gaus", fit_min, fit_max);
            TF1* f6 = new TF1("f6", "gaus", fit_min, fit_max);
            TF1* f7 = new TF1("f7", "gaus", fit_min, fit_max);
            TF1* f8 = new TF1("f8", "gaus", fit_min, fit_max);

            // 只对有数据的直方图进行拟合
            if(h_roocomb->GetEntries() > 10) h_roocomb->Fit(f1, "Q");
            if(h_TFF->GetEntries() > 10) h_TFF->Fit(f2, "Q");
            //if(h_SPD->GetEntries() > 10) h_SPD->Fit(f3, "Q");
            if(h_IWLS->GetEntries() > 10) h_IWLS->Fit(f3, "Q");
            if(h_DA->GetEntries() > 10) h_DA->Fit(f4, "Q");
            if(h_JSC->GetEntries() > 10) h_JSC->Fit(f5, "Q");
            if(h_asy->GetEntries() > 10) h_asy->Fit(f6, "Q");
            if(h_un_roocomb->GetEntries() > 10) h_un_roocomb->Fit(f7, "Q");
            if(h_roofit->GetEntries() > 10) h_roofit->Fit(f8, "Q");

            // 存储7种方法的拟合结果
            double scale_factor = (j == 0||j==2) ? 200.0 : 50.0;
            
            // 前4种方法
            sig_yield[j][0][i]       = h_roocomb->GetEntries() > 10 ? f1->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][0][i] = h_roocomb->GetEntries() > 10 ? f1->GetParError(1) / scale_factor : 0;
            sig_yield[j][1][i]       = h_TFF->GetEntries() > 10 ? f2->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][1][i] = h_TFF->GetEntries() > 10 ? f2->GetParError(1) / scale_factor : 0;
            //sig_yield[j][2][i]       = h_SPD->GetEntries() > 10 ? f3->GetParameter(1) / scale_factor : 0;
            //sig_yield_error[j][2][i] = h_SPD->GetEntries() > 10 ? f3->GetParError(1) / scale_factor : 0;
            sig_yield[j][2][i]       = h_IWLS->GetEntries() > 10 ? f3->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][2][i] = h_IWLS->GetEntries() > 10 ? f3->GetParError(1) / scale_factor : 0;
            
            // OM方法
            sig_yield[j][3][i]       = h_DA->GetEntries() > 10 ? f4->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][3][i] = h_DA->GetEntries() > 10 ? f4->GetParError(1) / scale_factor : 0;
            sig_yield[j][4][i]       = h_JSC->GetEntries() > 10 ? f5->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][4][i] = h_JSC->GetEntries() > 10 ? f5->GetParError(1) / scale_factor : 0;
            sig_yield[j][5][i]       = h_asy->GetEntries() > 10 ? f6->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][5][i] = h_asy->GetEntries() > 10 ? f6->GetParError(1) / scale_factor : 0;

            sig_yield[j][6][i]       = h_un_roocomb->GetEntries() > 10 ? f7->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][6][i] = h_un_roocomb->GetEntries() > 10 ? f7->GetParError(1) / scale_factor : 0;

            sig_yield[j][7][i]       = h_roofit->GetEntries() > 10 ? f8->GetParameter(1) / scale_factor : 0;
            sig_yield_error[j][7][i] = h_roofit->GetEntries() > 10 ? f8->GetParError(1) / scale_factor : 0;
            
            
            
            
            fit_sigma[j][0][i]       = h_roocomb->GetEntries() > 10 ? f1->GetParameter(2)/ scale_factor  : 0;
            fit_sigma_error[j][0][i] = h_roocomb->GetEntries() > 10 ? f1->GetParError(2) / scale_factor : 0;
            fit_sigma[j][1][i]       = h_TFF->GetEntries() > 10 ? f2->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][1][i] = h_TFF->GetEntries() > 10 ? f2->GetParError(2) / scale_factor : 0;
            //fit_sigma[j][2][i]       = h_SPD->GetEntries() > 10 ? f3->GetParameter(2) / scale_factor : 0;
            //fit_sigma_error[j][2][i] = h_SPD->GetEntries() > 10 ? f3->GetParError(2) / scale_factor : 0;
            fit_sigma[j][2][i]       = h_IWLS->GetEntries() > 10 ? f3->GetParameter(2)/ scale_factor  : 0;
            fit_sigma_error[j][2][i] = h_IWLS->GetEntries() > 10 ? f3->GetParError(2) / scale_factor : 0;
            fit_sigma[j][3][i]       = h_DA->GetEntries() > 10 ? f4->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][3][i] = h_DA->GetEntries() > 10 ? f4->GetParError(2) / scale_factor : 0;
            fit_sigma[j][4][i]       = h_JSC->GetEntries() > 10 ? f5->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][4][i] = h_JSC->GetEntries() > 10 ? f5->GetParError(2) / scale_factor : 0;
            fit_sigma[j][5][i]       = h_asy->GetEntries() > 10 ? f6->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][5][i] = h_asy->GetEntries() > 10 ? f6->GetParError(2) / scale_factor : 0;

            fit_sigma[j][6][i]       = h_un_roocomb->GetEntries() > 10 ? f7->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][6][i] = h_un_roocomb->GetEntries() > 10 ? f7->GetParError(2) / scale_factor : 0;
            fit_sigma[j][7][i]       = h_roofit->GetEntries() > 10 ? f8->GetParameter(2) / scale_factor : 0;
            fit_sigma_error[j][7][i] = h_roofit->GetEntries() > 10 ? f8->GetParError(2) / scale_factor : 0;


            // 清理内存
            delete h_roocomb;
            delete h_TFF;
            //delete h_SPD;
            delete h_IWLS;
            delete h_DA;
            delete h_JSC;
            delete h_asy;
            delete h_un_roocomb;
            delete h_roofit;
            delete f1;
            delete f2;
            delete f3;
            delete f4;
            delete f5;
            delete f6;
            delete f7;
            delete f8;

        }

        delete params_roocomb;
        delete params_TFF;
        //delete params_SPD;
        delete params_IWLS;
        delete params_un_roocomb;
        delete params_roofit;
        f_roocomb->Close();
        f_TFF    ->Close();
        //f_SPD    ->Close();
        f_IWLS   ->Close();
        f_OM->Close();
        f_unbinned_roocomb->Close();
        f_roofit->Close();
        delete f_roocomb;
        delete f_TFF;
        //delete f_SPD;
        delete f_IWLS;
        delete f_OM;
        delete f_unbinned_roocomb;
        delete f_roofit;
    }

    // 7种方法的名称和颜色
    const char* method_name[8] = {
        "RooCombineFit", 
        "TFF", 
        "IWLS", 
        "DA", 
        "JSC", 
        "asy",
        "unbinned_RooComb",
        "unbinned_RooFit"
    };
    
    int method_color[8] = {kRed, kBlue,  kMagenta, kCyan, kViolet, kGreen,kOrange+1, kPink};
    
    const char* sys[4] = {
        "Sig:Bkg=200:200 (bkg variation)",  // 背景情况1 - 信号数~200
        "Sig:Bkg=50:200 (bkg variation)" ,  // 背景情况2 - 信号数~50
        "Sig:Bkg=200:200 (sig variation)",  // 背景情况1 - 信号数~200
        "Sig:Bkg=50:200 (sig variation)" ,  // 背景情况2 - 信号数~50

    };
    bias_graph->cd();
    // 打开多页 PDF
    TCanvas* c_dummy = new TCanvas();
    c_dummy->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf[");
    delete c_dummy;

    // 对于每种背景系统atics情况
    for(int j = 0; j < nSys; ++j){
        // —— 信号产额散点图 —— 
        TCanvas* c_sy = new TCanvas(Form("c_sy_sys%d", j), Form("Bias Systematics %d", j), 800, 600);
        c_sy->SetLogx();
        c_sy->SetGridy();
        
        TLegend* leg_sy = new TLegend(0.6, 0.6, 0.9, 0.9);
        
        for(int m = 0; m < met; ++m){

            TGraphErrors* gr = new TGraphErrors(n);
            int validPoints = 0;
            
            for(int i = 0; i < n; ++i){
                // 只添加有效数据点
                if(sig_yield_error[j][m][i] > 0){
                    gr->SetPoint(validPoints, size[i], sig_yield[j][m][i]);
                    gr->SetPointError(validPoints, 0, sig_yield_error[j][m][i]);
                    validPoints++;
                }
            }
            
            if(validPoints > 0){
                gr->SetMarkerStyle(20 + m);
                gr->SetMarkerColor(method_color[m]);
                gr->SetLineColor(method_color[m]);
                
                if(m == 0){
                    // 根据背景情况设置不同的Y轴范围
                    if(j == 0) {
                        gr->GetYaxis()->SetRangeUser(0.9, 1.35);
                    } else {
                        gr->GetYaxis()->SetRangeUser(0.5, 2.3);
                    }
                    gr->Draw("AP");
                    gr->SetTitle(Form("Signal Yield Bias (%s);N_{MC};Bias", sys[j]));
                    gr->GetXaxis()->SetTitleSize(0.05);
                    gr->GetXaxis()->SetTitleOffset(0.85);
                    gr->GetYaxis()->SetTitleSize(0.05);
                    gr->GetYaxis()->SetTitleOffset(0.9);
                } else {
                    gr->SetTitle(Form("Signal Yield Bias (%s);N_{MC};Bias", sys[j]));
                    gr->Draw("P SAME");
                }
                leg_sy->AddEntry(gr, method_name[m], "p");
            }
            gr->SetName(Form("g_bias_%s_sys%d", method_name[m], j));
            gr->Write();
        }
        
        leg_sy->Draw();
        c_sy->Update();
        c_sy->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf");
        delete leg_sy;
        delete c_sy;

        // —— 高斯宽度散点图 —— 
        TCanvas* c_sig = new TCanvas(Form("c_sigma_sys%d", j), Form("Sigma Systematics %d", j), 800, 600);
        c_sig->SetLogx();
        c_sig->SetGridy();
        
        TLegend* leg_sig = new TLegend(0.6, 0.6, 0.9, 0.9);
        
        for(int m = 0; m < met; ++m){

            TGraphErrors* grs = new TGraphErrors(n);
            int validPoints = 0;
            
            for(int i = 0; i < n; ++i){
                // 只添加有效数据点
                if(fit_sigma_error[j][m][i] > 0){
                    grs->SetPoint(validPoints, size[i], fit_sigma[j][m][i]);
                    grs->SetPointError(validPoints, 0, fit_sigma_error[j][m][i]);
                    validPoints++;
                }
            }
            
            if(validPoints > 0){
                grs->SetMarkerStyle(24 + m);
                grs->SetMarkerColor(method_color[m]);
                grs->SetLineColor(method_color[m]);
                
                if(m == 0){
                    // 根据背景情况设置不同的Y轴范围
                    if(j == 0) {
                        grs->GetYaxis()->SetRangeUser(0, 0.5);
                    } else {
                        grs->GetYaxis()->SetRangeUser(0, 0.5);
                    }
                    grs->Draw("AP");
                    grs->SetTitle(Form("Signal Yield Sigma(Normalized) (%s);N_{MC};Sigma", sys[j]));
                    grs->GetXaxis()->SetTitleSize(0.05);
                    grs->GetXaxis()->SetTitleOffset(0.85);
                    grs->GetYaxis()->SetTitleSize(0.05);
                    grs->GetYaxis()->SetTitleOffset(0.9);
                } else {
                    grs->SetTitle(Form("Signal Yield Sigma(Normalized) (%s);N_{MC};Sigma", sys[j]));
                    grs->Draw("P SAME");
                }
                leg_sig->AddEntry(grs, method_name[m], "p");
            }
            grs->SetName(Form("g_sigma_%s_sys%d", method_name[m], j));
            grs->Write();
        }
        
        leg_sig->Draw();
        c_sig->Update();
        c_sig->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf");
        delete leg_sig;
        delete c_sig;
    }

    // 关闭 PDF
    TCanvas* c_end = new TCanvas();
    c_end->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf]");
    delete c_end;

    // // 保存结果到ROOT文件
    // bias_graph->cd();
    // for(int j = 0; j < nSys; ++j){
        
    //     for(int m = 0; m < met; ++m){

    //         TGraphErrors* gr_bias = new TGraphErrors(n);
    //         TGraphErrors* gr_sigma = new TGraphErrors(n);
            
    //         int validPoints = 0;
    //         for(int i = 0; i < n; ++i){
    //             if(sig_yield_error[j][m][i] > 0){
    //                 gr_bias->SetPoint(validPoints, size[i], sig_yield[j][m][i]);
    //                 gr_bias->SetPointError(validPoints, 0, sig_yield_error[j][m][i]);
    //                 gr_sigma->SetPoint(validPoints, size[i], fit_sigma[j][m][i]);
    //                 gr_sigma->SetPointError(validPoints, 0, fit_sigma_error[j][m][i]);
    //                 validPoints++;
    //             }
    //         }
            
    //         if(validPoints > 0){
    //             gr_bias->SetName(Form("g_bias_%s_sys%d", method_name[m], j));
    //             gr_sigma->SetName(Form("g_sigma_%s_sys%d", method_name[m], j));
    //             gr_bias->Write();
    //             gr_sigma->Write();
    //         }
            
    //         delete gr_bias;
    //         delete gr_sigma;
    //     }
    // }
    bias_graph->Close();
}