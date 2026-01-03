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
    // 配置
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,150,200,300,400,500,700,1000,1500,2000,3000,4000,5000,8000};
    const char* names[] = {
        "RooComb","TFF","IWLS","DA","JSC","asy","URComb","URFit",
        "KP0.7","KP0.9","KP1.0","KP1.1","KP1.3"
    };
    int n = size.size();
    const int base_methods = 8;
    const std::vector<double> widths = {0.7, 0.9, 1.0, 1.1, 1.3};
    int key_pdf_widths = widths.size();
    int met = base_methods + key_pdf_widths;
    int nSys = 4;

    // 数据容器
    auto Make3 = [&]( ){ return std::vector<std::vector<std::vector<double>>>(
        nSys, std::vector<std::vector<double>>(met, std::vector<double>(n, 0.))
    );};
    auto sig_yield     = Make3();
    auto sig_yield_err = Make3();
    auto fit_sigma     = Make3();
    auto fit_sigma_err = Make3();

    // ROOT 样式
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);

    // 输出文件
    TFile* bias_graph = new TFile("/data03/tianye/pbar/root/bias_graph/bias_graph.root", "RECREATE");
    bias_graph->cd();

    // 循环样本规模
    for(int i = 0; i < n; ++i){
        // 输入文件路径
        TString fn_binned_roocomb = Form("/data03/tianye/pbar/root/sys_result/binned_roocomb_seq0001_npi%04d.root", size[i]);
        TString fn_binned_TFF     = Form("/data03/tianye/pbar/root/sys_result/binned_TFF_seq0001_npi%04d.root",     size[i]);
        TString fn_binned_IWLS    = Form("/data03/tianye/pbar/root/sys_result/binned_IWLS_seq0001_npi%04d.root",   size[i]);
        TString fn_binned_OM      = Form("/data03/tianye/pbar/root/sys_result/binned_OM_seq0001_npi%04d.root",     size[i]);
        TString fn_un_roocomb     = Form("/data03/tianye/pbar/root/sys_result/unbinned_roocomb_seq0001_npi%04d.root", size[i]);
        TString fn_un_roofit      = Form("/data03/tianye/pbar/root/sys_result/unbinned_roofit_seq0001_npi%04d.root", size[i]);
        TString fn_un_keypdf      = Form("/data03/tianye/pbar/root/sys_result/unbinned_rookeyspdf_diffwidth_seq0001_npi%04d.root", size[i]);

        // 检查文件存在
        if(gSystem->AccessPathName(fn_binned_roocomb) ||
           gSystem->AccessPathName(fn_binned_TFF)     ||
           gSystem->AccessPathName(fn_binned_IWLS)    ||
           gSystem->AccessPathName(fn_binned_OM)      ||
           gSystem->AccessPathName(fn_un_roocomb)     ||
           gSystem->AccessPathName(fn_un_roofit)      ||
           gSystem->AccessPathName(fn_un_keypdf) ){
            std::cout<<"缺失文件，N="<<size[i]<<std::endl;
            continue;
        }

        // 打开输入文件
        TFile *f_roocomb = new TFile(fn_binned_roocomb, "READ");
        TFile *f_TFF     = new TFile(fn_binned_TFF,     "READ");
        TFile *f_IWLS    = new TFile(fn_binned_IWLS,    "READ");
        TFile *f_OM      = new TFile(fn_binned_OM,      "READ");
        TFile *f_unRC    = new TFile(fn_un_roocomb,     "READ");
        TFile *f_unRF    = new TFile(fn_un_roofit,      "READ");
        TFile *f_keypdf  = new TFile(fn_un_keypdf,      "READ");

        // 取树
        TTree *t_roocomb = (TTree*)f_roocomb->Get("fitParamsTree");
        TTree *t_TFF     = (TTree*)f_TFF    ->Get("fitParamsTree");
        TTree *t_IWLS    = (TTree*)f_IWLS   ->Get("fitParamsTree");
        TTree *t_OM      = (TTree*)f_OM     ->Get("fitParamsTree");
        TTree *t_unRC    = (TTree*)f_unRC   ->Get("fitParamsTree");
        TTree *t_unRF    = (TTree*)f_unRF   ->Get("fitParamsTree");
        TTree *t_keypdf  = (TTree*)f_keypdf ->Get("fitParamsTree");
        if(!t_roocomb||!t_TFF||!t_IWLS||!t_OM||!t_unRC||!t_unRF||!t_keypdf){
            std::cout<<"缺失树，N="<<size[i]<<std::endl;
            for(auto f:{f_roocomb,f_TFF,f_IWLS,f_OM,f_unRC,f_unRF,f_keypdf}) f->Close();
            continue;
        }

        // 参数容器
        std::vector<double>* p_roocomb = new std::vector<double>(4);
        std::vector<double>* p_TFF     = new std::vector<double>(4);
        std::vector<double>* p_IWLS    = new std::vector<double>(4);
        Double_t pDA[2]={0}, pJSC[2]={0}, pasy[2]={0}, pDA1[2]={0}, pJSC1[2]={0}, pasy1[2]={0};
        std::vector<double>* p_unRC    = new std::vector<double>(4);
        std::vector<double>* p_unRF    = new std::vector<double>(4);
        std::vector<std::vector<double>>* p_keypdf = new std::vector<std::vector<double>>(widths.size());

        // 绑定
        t_roocomb->SetBranchAddress("nsig", &p_roocomb);
        t_TFF    ->SetBranchAddress("nsig", &p_TFF);
        t_IWLS   ->SetBranchAddress("nsig", &p_IWLS);
        t_unRC   ->SetBranchAddress("nsig", &p_unRC);
        t_unRF   ->SetBranchAddress("nsig", &p_unRF);
        t_keypdf ->SetBranchAddress("nsig", &p_keypdf);
        t_OM->SetBranchAddress("Ysig_DA",   pDA);
        t_OM->SetBranchAddress("Ysig_JSC",  pJSC);
        t_OM->SetBranchAddress("Ysig_asy",  pasy);
        t_OM->SetBranchAddress("Ysig1_DA",  pDA1);
        t_OM->SetBranchAddress("Ysig1_JSC", pJSC1);
        t_OM->SetBranchAddress("Ysig1_asy", pasy1);

        int entries = t_roocomb->GetEntries();

        // 系统atics 循环
        for(int j = 0; j < nSys; ++j){
            double hmin = 0, hmax = ((j==0||j==2)?400:200), fmin = hmin, fmax = hmax;
            int bins = 300;

            // 建立所有 TH1D
            std::vector<TH1D*> hs(met);
            for(int m=0; m<met; ++m){
                TString tag = (m<base_methods
                    ? Form("%s", names[m])
                    : Form("KP%02d", int(widths[m - base_methods] * 10)));
                hs[m] = new TH1D(Form("h_%s_%d_%d", tag.Data(), i, j),
                                 Form("%s;Fitted Signal;Count", tag.Data()),
                                 bins, hmin, hmax);
            }

            // 填充
            for(int k=0; k<entries; ++k){
                t_roocomb->GetEntry(k);
                t_TFF    ->GetEntry(k);
                t_IWLS   ->GetEntry(k);
                t_unRC   ->GetEntry(k);
                t_unRF   ->GetEntry(k);
                t_OM     ->GetEntry(k);
                t_keypdf ->GetEntry(k);

                // 前八种
                hs[0]->Fill((*p_roocomb)[j]);
                hs[1]->Fill((*p_TFF)[j]*((j==0||j==2)?400:250));
                hs[2]->Fill((*p_IWLS)[j]);
                double *dA = (j<2? pDA:pDA1), *jsc=(j<2? pJSC:pJSC1), *asyarr=(j<2? pasy:pasy1);
                int idx = (j<2? j: j-2);
                hs[3]->Fill(dA[idx]);
                hs[4]->Fill(jsc[idx]);
                hs[5]->Fill(asyarr[idx]);
                hs[6]->Fill((*p_unRC)[j]);
                hs[7]->Fill((*p_unRF)[j]);

                // 五个 width
                for(int w=0; w<key_pdf_widths; ++w){
                    hs[base_methods + w]->Fill((*p_keypdf)[w][j]);
                }
            }

            // 写入并拟合
            for(int m=0; m<met; ++m){
                bias_graph->cd();
                hs[m]->Write();
                TF1 *f = new TF1(Form("f_%d_%d", m, j), "gaus", fmin, fmax);
                if(hs[m]->GetEntries()>10) hs[m]->Fit(f,"Q");
                double scale = ((j==0||j==2)?200:50);
                double mean = (hs[m]->GetEntries()>10? f->GetParameter(1)/scale: 0.);
                double errm = (hs[m]->GetEntries()>10? f->GetParError(1)/scale: 0.);
                double sigm = (hs[m]->GetEntries()>10? f->GetParameter(2)/scale: 0.);
                double errs = (hs[m]->GetEntries()>10? f->GetParError(2)/scale: 0.);
                sig_yield[j][m][i]     = mean;
                sig_yield_err[j][m][i] = errm;
                fit_sigma[j][m][i]     = sigm;
                fit_sigma_err[j][m][i] = errs;
                delete f;
                delete hs[m];
            }
        }

        // 清理
        delete p_roocomb; delete p_TFF; delete p_IWLS;
        delete p_unRC; delete p_unRF; delete p_keypdf;
        for(auto f:{f_roocomb,f_TFF,f_IWLS,f_OM,f_unRC,f_unRF,f_keypdf}) f->Close();
    }

    // 绘制多页 PDF

    Int_t cols[] = {
        kRed,kBlue,kMagenta,kCyan,kViolet,kGreen,kOrange+1,kPink,
        kRed-7,kBlue-7,kMagenta-7,kCyan-7,kViolet-7
    };
    const char* slabel[] = {
        "200:200(b)","50:200(b)","200:200(s)","50:200(s)"
    };

    TCanvas* c0 = new TCanvas(); 
    c0->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf[");
    delete c0;

    for(int j=0; j<nSys; ++j){
        // Bias
        TCanvas* cb = new TCanvas(); cb->SetLogx(); cb->SetGridy();
        TLegend* lb = new TLegend(0.6,0.6,0.9,0.9);
        for(int m=0; m<met; ++m){
            TGraphErrors* g = new TGraphErrors();
            int pt=0;
            for(int i=0;i<n;++i){
                if(sig_yield_err[j][m][i]>0){
                    g->SetPoint(pt, size[i], sig_yield[j][m][i]);
                    g->SetPointError(pt, 0, sig_yield_err[j][m][i]);
                    ++pt;
                }
            }
            if(pt>0){
                g->SetMarkerStyle(20+m);
                g->SetMarkerColor(cols[m]);
                g->SetLineColor(cols[m]);
                if(m==0){
                    if(j==0) g->GetYaxis()->SetRangeUser(0.96,1.35);
                    else if(j==1) g->GetYaxis()->SetRangeUser(0.9,2.3);
                    else if(j==2) g->GetYaxis()->SetRangeUser(0.5,1.1);
                    else     g->GetYaxis()->SetRangeUser(0.5,1.2);
                    g->Draw("AP");
                    g->SetTitle(Form("Bias (%s);N_{MC};Bias", slabel[j]));
                } else {
                    g->SetTitle(Form("Bias (%s);N_{MC};Bias", slabel[j]));
                    g->Draw("P SAME");
                }
                lb->AddEntry(g, names[m], "p");
            }
            g->SetName(Form("g_bias_%s_sys%d", names[m], j));
            g->Write();
        }
        lb->Draw();
        cb->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf");
        delete lb; delete cb;

        // Sigma
        TCanvas* cs = new TCanvas(); cs->SetLogx(); cs->SetGridy();
        TLegend* ls = new TLegend(0.6,0.6,0.9,0.9);
        for(int m=0; m<met; ++m){
            TGraphErrors* g = new TGraphErrors();
            int pt=0;
            for(int i=0;i<n;++i){
                if(fit_sigma_err[j][m][i]>0){
                    g->SetPoint(pt, size[i], fit_sigma[j][m][i]);
                    g->SetPointError(pt, 0, fit_sigma_err[j][m][i]);
                    ++pt;
                }
            }
            if(pt>0){
                g->SetMarkerStyle(24+m);
                g->SetMarkerColor(cols[m]);
                g->SetLineColor(cols[m]);
                if(m==0){
                    g->GetYaxis()->SetRangeUser(0,0.5);
                    g->Draw("AP");
                    g->SetTitle(Form("Sigma (%s);N_{MC};Sigma", slabel[j]));
                } else {
                    g->SetTitle(Form("Sigma (%s);N_{MC};Sigma", slabel[j]));
                    g->Draw("P SAME");
                }
                ls->AddEntry(g, names[m], "p");
            }
            g->SetName(Form("g_sigma_%s_sys%d", names[m], j));
        }
        ls->Draw();
        cs->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf");
        delete ls; delete cs;
    }

    TCanvas* c1 = new TCanvas();
    c1->Print("/data03/tianye/pbar/syscode/pdf/binned_bias_plots_with_syserror.pdf]");
    delete c1;

    bias_graph->Close();
}
