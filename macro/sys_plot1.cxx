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
    int met = 9; // 9种方法

    // 三类数据 × 九种方法 × n
    double sig_yield[1][met][n];
    double sig_yield_error[1][met][n];
    double fit_sigma[1][met][n];
    double fit_sigma_error[1][met][n];

    // ROOT 样式
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TFile* bias_graph= new TFile("/data03/tianye/pbar/root/bias_graph/bias_graph.root", "RECREATE");

    // 逐个事例数读取并拟合
    for(int i = 0; i < n; ++i){
        TFile *f_roocomb = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_roocomb_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_TFF     = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_TFF_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_SPD     = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_SPD_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_IWLS    = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_IWLS_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_roocomb_unbinned = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roocomb_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_roofit_unbinned = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roofit_seq0001_npi%04d.root", size[i]), "read");
        TFile *f_OM= new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_OM_seq0001_npi%04d.root", size[i]), "read");

        TTree *t_roocomb = (TTree*)f_roocomb->Get("fitParamsTree");
        TTree *t_TFF     = (TTree*)f_TFF    ->Get("fitParamsTree");
        TTree *t_SPD     = (TTree*)f_SPD    ->Get("fitParamsTree");
        TTree *t_IWLS    = (TTree*)f_IWLS   ->Get("fitParamsTree");
        TTree *t_roocomb_unbinned = (TTree*)f_roocomb_unbinned->Get("fitParamsTree");
        TTree *t_roofit_unbinned = (TTree*)f_roofit_unbinned->Get("fitParamsTree");
        TTree *t_OM = (TTree*)f_OM->Get("fitParamsTree");

        std::vector<double>* params_roocomb = new std::vector<double>(4);
        std::vector<double>* params_TFF     = new std::vector<double>(4);
        std::vector<double>* params_SPD     = new std::vector<double>(4);
        std::vector<double>* params_IWLS    = new std::vector<double>(4);
        std::vector<double>* params_roocomb_unbinned = new std::vector<double>(4);
        std::vector<double>* params_roofit_unbinned = new std::vector<double>(4);
        std::vector<double>* params_DA = new std::vector<double>(2);
        std::vector<double>* params_JSC = new std::vector<double>(2);
        std::vector<double>* params_asy = new std::vector<double>(2);


        t_roocomb->SetBranchAddress("nsig", &params_roocomb);
        t_TFF    ->SetBranchAddress("nsig", &params_TFF);
        t_SPD    ->SetBranchAddress("nsig", &params_SPD);
        t_IWLS   ->SetBranchAddress("nsig", &params_IWLS);
        t_roocomb_unbinned->SetBranchAddress("nsig", &params_roocomb_unbinned);
        t_roofit_unbinned->SetBranchAddress("nsig", &params_roofit_unbinned);
        t_OM->SetBranchAddress("Ysig_DA", &params_DA);
        t_OM->SetBranchAddress("Ysig_JSC", &params_JSC);
        t_OM->SetBranchAddress("Ysig_asy", &params_asy);

        int nEntries = t_roocomb->GetEntries();

        for(int j = 0; j < 1; ++j){
            TH1D* h_roocomb = new TH1D("h_roocomb", "", 200, 140, 280);
            TH1D* h_TFF     = new TH1D("h_TFF",     "", 200, 140, 280);
            TH1D* h_SPD     = new TH1D("h_SPD",     "", 200, 140, 280);
            TH1D* h_IWLS    = new TH1D("h_IWLS",    "", 200, 140, 280);
            TH1D* h_roocomb_unbinned = new TH1D("h_roocomb_unbinned", "", 200, 140, 280);
            TH1D* h_roofit_unbinned = new TH1D("h_roofit_unbinned", "", 200, 140, 280);
            TH1D* h_DA= new TH1D("h_DA", "h_DA", 200, 140, 280);
            TH1D* h_JSC= new TH1D("h_JSC", "h_JSC", 200, 140, 280);
            TH1D* h_asy= new TH1D("h_asy", "h_asy", 200, 140, 280);

            for(int k = 0; k < nEntries; ++k){
                t_roocomb->GetEntry(k);
                t_TFF    ->GetEntry(k);
                t_SPD    ->GetEntry(k);
                t_IWLS   ->GetEntry(k);
                t_roocomb_unbinned->GetEntry(k);
                t_roofit_unbinned->GetEntry(k);
                t_OM->GetEntry(k);
                h_roocomb->Fill((*params_roocomb)[j] );
                h_TFF    ->Fill((*params_TFF)[j]  *400 );
                h_SPD    ->Fill((*params_SPD)[j]     );
                h_IWLS   ->Fill((*params_IWLS)[j]   );
                h_roocomb_unbinned->Fill((*params_roocomb_unbinned)[j] * 400);
                h_roofit_unbinned->Fill((*params_roofit_unbinned)[j]);
                h_DA->Fill(params_DA );
                h_JSC->Fill(params_JSC );
                h_asy->Fill(params_asy );
            }

            TF1* f1 = new TF1("f1", "gaus", 100, 300);
            TF1* f2 = new TF1("f2", "gaus", 100, 300);
            TF1* f3 = new TF1("f3", "gaus", 100, 300);
            TF1* f4 = new TF1("f4", "gaus", 100, 300);
            TF1* f5 = new TF1("f5", "gaus", 100, 300);
            TF1* f6 = new TF1("f6", "gaus", 100, 300);
            TF1* f7 = new TF1("f7", "gaus", 100, 300);
            TF1* f8 = new TF1("f8", "gaus", 100, 300);
            TF1* f9 = new TF1("f9", "gaus", 100, 300);

            h_roocomb->Fit(f1, "Q");
            h_TFF    ->Fit(f2, "Q");
            h_SPD    ->Fit(f3, "Q");
            h_IWLS   ->Fit(f4, "Q");
            h_roocomb_unbinned->Fit(f5, "Q");
            h_roofit_unbinned->Fit(f6, "Q");
            h_DA->Fit(f7, "Q");
            h_JSC->Fit(f8, "Q");
            h_asy->Fit(f9, "Q");

            // 存储9种方法的拟合结果
            sig_yield[j][0][i]       = f1->GetParameter(1) / 200.0;
            sig_yield_error[j][0][i] = f1->GetParError(1) / 200.0;
            sig_yield[j][1][i]       = f2->GetParameter(1) / 200.0;
            sig_yield_error[j][1][i] = f2->GetParError(1) / 200.0;
            sig_yield[j][2][i]       = f3->GetParameter(1) / 200.0;
            sig_yield_error[j][2][i] = f3->GetParError(1) / 200.0;
            sig_yield[j][3][i]       = f4->GetParameter(1) / 200.0;
            sig_yield_error[j][3][i] = f4->GetParError(1) / 200.0;
            sig_yield[j][4][i]       = f5->GetParameter(1) / 200.0;
            sig_yield_error[j][4][i] = f5->GetParError(1) / 200.0;
            sig_yield[j][5][i]       = f6->GetParameter(1) / 200.0;
            sig_yield_error[j][5][i] = f6->GetParError(1) / 200.0;
            sig_yield[j][6][i]       = f7->GetParameter(1) / 200.0;
            sig_yield_error[j][6][i] = f7->GetParError(1) / 200.0;
            sig_yield[j][7][i]       = f8->GetParameter(1) / 200.0;
            sig_yield_error[j][7][i] = f8->GetParError(1) / 200.0;
            sig_yield[j][8][i]       = f9->GetParameter(1) / 200.0;
            sig_yield_error[j][8][i] = f9->GetParError(1) / 200.0;

            fit_sigma[j][0][i]       = f1->GetParameter(2);
            fit_sigma_error[j][0][i] = f1->GetParError(2);
            fit_sigma[j][1][i]       = f2->GetParameter(2);
            fit_sigma_error[j][1][i] = f2->GetParError(2);
            fit_sigma[j][2][i]       = f3->GetParameter(2);
            fit_sigma_error[j][2][i] = f3->GetParError(2);
            fit_sigma[j][3][i]       = f4->GetParameter(2);
            fit_sigma_error[j][3][i] = f4->GetParError(2);
            fit_sigma[j][4][i]       = f5->GetParameter(2);
            fit_sigma_error[j][4][i] = f5->GetParError(2);
            fit_sigma[j][5][i]       = f6->GetParameter(2);
            fit_sigma_error[j][5][i] = f6->GetParError(2);
            fit_sigma[j][6][i]       = f7->GetParameter(2);
            fit_sigma_error[j][6][i] = f7->GetParError(2);
            fit_sigma[j][7][i]       = f8->GetParameter(2);
            fit_sigma_error[j][7][i] = f8->GetParError(2);
            fit_sigma[j][8][i]       = f9->GetParameter(2);
            fit_sigma_error[j][8][i] = f9->GetParError(2);

            // 清理内存
            delete h_roocomb;
            delete h_TFF;
            delete h_SPD;
            delete h_IWLS;
            delete h_roocomb_unbinned;
            delete h_roofit_unbinned;
            delete h_DA;
            delete h_JSC;
            delete h_asy;
            delete f1;
            delete f2;
            delete f3;
            delete f4;
            delete f5;
            delete f6;
            delete f7;
            delete f8;
            delete f9;
        }

        delete params_roocomb;
        delete params_TFF;
        delete params_SPD;
        delete params_IWLS;
        delete params_roocomb_unbinned;
        delete params_roofit_unbinned;
        f_roocomb->Close();
        f_TFF    ->Close();
        f_SPD    ->Close();
        f_IWLS   ->Close();
        f_roocomb_unbinned->Close();
        f_roofit_unbinned->Close();
        f_OM->Close();
        delete f_roocomb;
        delete f_TFF;
        delete f_SPD;
        delete f_IWLS;
        delete f_roocomb_unbinned;
        delete f_roofit_unbinned;
        delete f_OM;
    }

    // 9种方法的名称和颜色
    const char* method_name[9] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS", "U_RooCombineFit", "U_RooFit", "DA", "JSC", "asy"};
    int method_color[9] = {kRed, kBlue, kGreen, kMagenta, kPink, kOrange, kCyan, kViolet, kGray};
    const char* sys[5] = {
        "Original Background ",
        "alpha_R: 0.530 -> 0.400",
        "alpha_R: 0.530 -> 0.660",
        "sigma_R: 0.115 -> 0.095",
        "sigma_R: 0.115 -> 0.140"
    };

    // 打开多页 PDF
    TCanvas* c_dummy = new TCanvas();
    c_dummy->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf[");
    delete c_dummy;

    std::vector<TGraphErrors*> g_sy;
    std::vector<TGraphErrors*> g_sig;

    for(int j = 0; j < 1; ++j){
        // —— 信号产额散点图 —— 
        TCanvas* c_sy = new TCanvas(Form("c_sy_cat%d", j), "", 800, 600);
        c_sy->SetLogx();
        c_sy->SetGridy();
        
        for(int m = 0; m < met; ++m){
            TGraphErrors* gr = new TGraphErrors(n);
            for(int i = 0; i < n; ++i){
                gr->SetPoint(i, size[i], sig_yield[j][m][i]);
                gr->SetPointError(i, 0, sig_yield_error[j][m][i]);
            }
            gr->SetMarkerStyle(20 + m);
            gr->SetMarkerColor(method_color[m]);
            gr->SetLineColor(method_color[m]);
            gr->SetTitle(method_name[m]);
            if(m == 0){
                gr->GetYaxis()->SetRangeUser(0.9, 1.3);
                gr->Draw("APE");
                gr->SetTitle(Form("SignalYield Bias vs Nmc;Nmc;Bias", j, sys[j]));
                gr->GetXaxis()->SetTitleSize(0.05);
                gr->GetXaxis()->SetTitleOffset(0.85);
                gr->GetYaxis()->SetTitleSize(0.05);
                gr->GetYaxis()->SetTitleOffset(0.9);
            } else {
                gr->Draw("P E SAME");
            }
            g_sy.push_back(gr);
        }
        TLegend* leg_sy = new TLegend(0.5, 0.5, 0.9, 0.9);
        for(int m = 0; m < met; ++m){
            leg_sy->AddEntry(g_sy[m], method_name[m], "p");
        }
        leg_sy->Draw();
        c_sy->Update();
        c_sy->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf");
        delete leg_sy;
        delete c_sy;

        // —— 高斯宽度散点图 —— 
        TCanvas* c_sig = new TCanvas(Form("c_sig_cat%d", j), "", 800, 600);
        c_sig->SetLogx();
        c_sig->SetGridy();
        
        for(int m = 0; m < met; ++m){
            TGraphErrors* grs = new TGraphErrors(n);
            for(int i = 0; i < n; ++i){
                grs->SetPoint(i, size[i], fit_sigma[j][m][i]);
                grs->SetPointError(i, 0, fit_sigma_error[j][m][i]);
            }
            grs->SetMarkerStyle(24 + m);
            grs->SetMarkerColor(method_color[m]);
            grs->SetLineColor(method_color[m]);
            grs->SetTitle(method_name[m]);
            if(m == 0){
                grs->GetYaxis()->SetRangeUser(4.0, 24.0);
                grs->Draw("APE");
                grs->SetTitle(Form("SignalYield Sigma vs Nmc;Nmc;Sigma", j));
                grs->GetXaxis()->SetTitleSize(0.05);
                grs->GetXaxis()->SetTitleOffset(0.85);
                grs->GetYaxis()->SetTitleSize(0.05);
                grs->GetYaxis()->SetTitleOffset(0.9);
            } else {
                grs->Draw("P E SAME");
            }
            g_sig.push_back(grs);
        }
        TLegend* leg_sig = new TLegend(0.5, 0.5, 0.9, 0.9);
        for(int m = 0; m < met; ++m){
            leg_sig->AddEntry(g_sig[m], method_name[m], "p");
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

    // 保存结果到ROOT文件
    bias_graph->cd();
    for(int m=0; m<met; ++m){
        g_sy[m]->Write(Form("g_bias_%s", method_name[m]));
        g_sig[m]->Write(Form("g_sigma_%s", method_name[m]));
    }
    bias_graph->Close();

    // 清理内存
    for(auto gr : g_sy) delete gr;
    for(auto gr : g_sig) delete gr;
}