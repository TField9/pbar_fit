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

void sig_yield1(){
    std::vector<int> size = {100, 150, 200, 500, 700, 1000, 1500, 2000, 5000, 8000};
    std::vector<std::vector<TH1D*>> all_sig_yield;
    std::vector<std::vector<TH1D*>> all_errors;
    std::vector<std::vector<double>> all_bias;
    std::vector<std::vector<double>> all_bias_error;
    std::vector<std::vector<TF1*>> all_fit;
    std::vector<std::vector<TH1D*>> all_pulls;
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);
    c1->Divide(3,4);
    c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf[");

    for (int i : size){ 
            // gStyle->SetLabelSize(0.1, "X");
            std::vector<TH1D*> sig_yield;
            std::vector<TH1D*> sig_yield_error;
            std::vector<TH1D*> sig_yield_pull;
            std::vector<double> bias;
            std::vector<double> bias_error;
            //TFile *f=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_roocomb_seq0001_npi%04d.root",i),"read");
            //TFile *f=new TFile(Form("/data03/tianye/pbar/root/bin25_seed1_nllre/binned_roocomb_seq0001_npi%04d.root",i),"read");
            TFile *f=new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_roocomb_seq0001_npi%04d.root",i),"read");
            //TFile *f=new TFile(Form("/data03/tianye/pbar/root/bin25_seed1_nllre1/binned_roocomb_seq0001_npi%04d.root",i),"read");
            TFile *f1=new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_TFF_seq0001_npi%04d.root",i),"read");
            // TFile *f2=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_SPD_seq0001_npi%04d.root",i),"read");
            TFile *f2=new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_SPD_seq0001_npi%04d.root",i),"read");
            TFile *f3=new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_IWLS1_seq0001_npi%04d.root",i),"read");
            TFile *f4=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/unbinned_roocomb_seq0001_npi%04d.root",i),"read");
            TFile *f5=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/unbinned_roofit_seq0001_npi%04d.root",i),"read");

            std::vector<double>* params1 = new std::vector<double>(22);
            std::vector<double>* params2 = new std::vector<double>(22);
            std::vector<double>* params3 = new std::vector<double>(4);
            // std::vector<double>* params4 = new std::vector<double>(4);
            // std::vector<double>* params5 = new std::vector<double>(4);
            float params4[4];
            float params5[4];
            std::vector<double>* params6 = new std::vector<double>(4);
            std::vector<double>* params7 = new std::vector<double>(22);
            std::vector<double>* params8 = new std::vector<double>(2);

            TTree *t = (TTree*)f->Get("fitParamsTree");
            TTree *t1 = (TTree*)f1->Get("fitParamsTree");
            TTree *t2 = (TTree*)f2->Get("fitParamsTree");
            TTree* t3 = (TTree*)f3->Get("fitParamsTree");
            TTree* t4 = (TTree*)f4->Get("fitParamsTree");
            TTree* t5 = (TTree*)f5->Get("fitParamsTree");

            t->SetBranchAddress("fitparams_bin30", &params1);
            t->SetBranchAddress("fitparams_bin100", &params2);
            t1->SetBranchAddress("nsig", &params3);
            t2->SetBranchAddress("bin30_SPD_pararr", &params4);
            t2->SetBranchAddress("bin100_SPD_pararr", &params5);
            t3->SetBranchAddress("nsig", &params6);
            t4->SetBranchAddress("fitparams_unbinned", &params7);
            t5->SetBranchAddress("fitparams_unbinned", &params8);

            TH1D* bin30_roocomb = new TH1D("roocomb", "RooCombineFit", 200, 140, 280);
            TH1D* bin100_roocomb = new TH1D("bin100_roocomb", "bin100_roocomb", 200, 140, 280);
            TH1D* bin30_TFF = new TH1D("TFF", "TFF", 200, 140, 280);
            TH1D* bin100_TFF = new TH1D("bin100_TFF", "bin100_TFF", 200, 140, 280);
            TH1D* bin30_SPD = new TH1D("Dembinski", "ALM-TF", 200, 140, 280);
            TH1D* bin100_SPD = new TH1D("bin100_Dembinski", "bin100_Dembinski", 200, 160, 300);
            TH1D* bin30_IWLS = new TH1D("IWLS", "IWLS", 200, 140, 280);
            TH1D* bin100_IWLS = new TH1D("bin100_IWLS", "bin100_IWLS", 200, 140, 280);
            TH1D* unbinned_comb = new TH1D("unbinned_comb", "Unbinned RooComb", 200, 140, 280);
            TH1D* unbinned_roofit = new TH1D("unbinned_roofit", "Unbinned RooFit", 200, 140, 280);

            TH1D* bin30_roocomb_error= new TH1D("roocomb_error", "RooCombineFit Error", 200, 10, 25);
            TH1D* bin100_roocomb_error= new TH1D("bin100_roocomb_error", "bin100_roocomb_error", 200, 10, 30);
            TH1D* bin30_TFF_error= new TH1D("TFF_error", "TFF Error", 200, 110, 150);
            TH1D* bin100_TFF_error= new TH1D("bin100_TFF_error", "bin100_TFF_error", 200, 50, 100);
            TH1D* bin30_SPD_error= new TH1D("Dembinski_error", "ALM-TF Error", 200, 10, 20);
            TH1D* bin100_SPD_error= new TH1D("bin100_Dembinski_error", "bin100_Dembinski_error", 200, 20, 30);
            TH1D* bin30_IWLS_error= new TH1D("IWLS_error", "IWLS Error", 200, 65, 100);
            TH1D* bin100_IWLS_error= new TH1D("bin100_IWLS_error", "bin100_IWLS_error", 200, 35, 55);
            TH1D* unbinned_comb_error= new TH1D("unbinned_comb_error", "Unbinned RooComb Error", 200, 10, 30);
            TH1D* unbinned_roofit_error= new TH1D("unbinned_roofit_error", "Unbinned RooFit Error", 200, 10, 18);

            TH1D* bin30_roocomb_pull= new TH1D("RooComb Pull", "RooComb Pull", 200, -3, 3);
            TH1D* bin100_roocomb_pull= new TH1D("bin100_roocomb_pull", "bin100_roocomb_pull", 200, -3, 3);
            TH1D* bin30_TFF_pull= new TH1D("TFF Pull", "TFF Pull", 200, -3, 3);
            TH1D* bin100_TFF_pull= new TH1D("bin100_TFF_pull", "bin100_TFF_pull", 200, -3, 3);
            TH1D* bin30_SPD_pull= new TH1D("ALM-TF Pull", "ALM-TF Pull", 200, -3, 3);
            TH1D* bin100_SPD_pull= new TH1D("bin100_Dembinski_pull", "bin100_Dembinski_pull", 200, -3, 3);
            TH1D* bin30_IWLS_pull= new TH1D("IWLS Pull", "IWLS Pull", 200, -3, 3);
            TH1D* bin100_IWLS_pull= new TH1D("bin100_IWLS_pull", "bin100_IWLS_pull", 200, -3, 3);
            TH1D* unbinned_comb_pull= new TH1D("Unbinned RooComb Pull", "Unbinned RooComb Pull", 200, -3, 3);
            TH1D* unbinned_roofit_pull= new TH1D("Unbinned RooFit Pull", "Unbinned RooFit Pull", 200, -3, 3);


                int nEntries = t->GetEntries();
                for (int i = 0; i < nEntries; ++i) {
                    t->GetEntry(i);
                    t1->GetEntry(i);
                    t2->GetEntry(i);
                    t3->GetEntry(i);
                    t4->GetEntry(i);
                    t5->GetEntry(i);

                    bin30_roocomb->Fill((*params1)[10] * 400);
                    bin100_roocomb->Fill((*params2)[10] * 400);
                    bin30_TFF->Fill((*params3)[0] * 400);
                    bin100_TFF->Fill((*params3)[1] * 400);
                    bin30_SPD->Fill(params4[0] );
                    bin100_SPD->Fill(params5[0] );
                    bin30_IWLS->Fill((*params6)[0] * 400);
                    bin100_IWLS->Fill((*params6)[1] * 400);
                    unbinned_comb->Fill((*params7)[10] * 400);
                    unbinned_roofit->Fill((*params8)[0]);


                    bin30_roocomb_error->Fill((*params1)[21] * 400);
                    bin30_TFF_error->Fill((*params3)[2] * 400);
                    bin30_SPD_error->Fill(params4[2] );
                    bin30_IWLS_error->Fill((*params6)[2] * 400);
                    unbinned_comb_error->Fill((*params7)[21] * 400);
                    unbinned_roofit_error->Fill((*params8)[1]);
                    bin100_roocomb_error->Fill((*params2)[21] * 400);
                    bin100_TFF_error->Fill((*params3)[3] * 400);
                    bin100_SPD_error->Fill(params5[2] * 400);
                    bin100_IWLS_error->Fill((*params6)[3] * 400);

                    bin30_roocomb_pull->Fill(((*params1)[10] - 0.5)/(*params1)[21] );
                    bin30_TFF_pull->Fill(((*params3)[0] - 0.5)/(*params3)[2] );
                    bin30_SPD_pull->Fill((params4[0] - 200)/params4[2] );
                    bin30_IWLS_pull->Fill(((*params6)[0] - 0.5)/(*params6)[2] );
                    unbinned_comb_pull->Fill(((*params7)[10] - 0.5)/(*params7)[21] );
                    unbinned_roofit_pull->Fill(((*params8)[0] - 200)/(*params8)[1] );
                    bin100_roocomb_pull->Fill(((*params2)[10] - 0.5)/(*params2)[21] );
                    bin100_TFF_pull->Fill(((*params3)[1] - 0.5)/(*params3)[3] );
                    bin100_SPD_pull->Fill((params5[0] - 200)/params5[2] );
                    bin100_IWLS_pull->Fill(((*params6)[1] - 0.5)/(*params6)[3] );
                }

                // Fit histograms with Gaussian
                TF1* gaussFit1 = new TF1("gaussFit1", "gaus", 100, 300);
                TF1* gaussFit2 = new TF1("gaussFit2", "gaus", 100, 300);
                TF1* gaussFit3 = new TF1("gaussFit3", "gaus", 100, 300);
                TF1* gaussFit4 = new TF1("gaussFit4", "gaus", 100, 300);
                TF1* gaussFit5 = new TF1("gaussFit5", "gaus", 100, 300);
                TF1* gaussFit6 = new TF1("gaussFit6", "gaus", 100, 300);
                TF1* gaussFit7 = new TF1("gaussFit7", "gaus", 100, 300);
                TF1* gaussFit8 = new TF1("gaussFit8", "gaus", 100, 300);
                TF1* gaussFit9 = new TF1("gaussFit9", "gaus", 100, 300);
                TF1* gaussFit10 = new TF1("gaussFit10", "gaus", 100, 300);

                bin30_roocomb->Fit(gaussFit1, "Q");
                bin30_TFF->Fit(gaussFit2, "Q");
                bin30_SPD->Fit(gaussFit3, "Q");
                bin30_IWLS->Fit(gaussFit4, "Q");
                unbinned_comb->Fit(gaussFit5, "Q");
                unbinned_roofit->Fit(gaussFit6, "Q");
                bin100_roocomb->Fit(gaussFit7, "Q");
                bin100_TFF->Fit(gaussFit8, "Q");
                bin100_SPD->Fit(gaussFit9, "Q");
                bin100_IWLS->Fit(gaussFit10, "Q");

                sig_yield.push_back(bin30_roocomb);
                bias.push_back(gaussFit1->GetParameter(1) / 200);
                bias_error.push_back(gaussFit1->GetParError(1) / 200);

                sig_yield.push_back(bin30_TFF);
                bias.push_back(gaussFit2->GetParameter(1) / 200);
                bias_error.push_back(gaussFit2->GetParError(1) / 200);

                sig_yield.push_back(bin30_SPD);
                bias.push_back(gaussFit3->GetParameter(1) / 200);
                bias_error.push_back(gaussFit3->GetParError(1) / 200);

                sig_yield.push_back(bin30_IWLS);
                bias.push_back(gaussFit4->GetParameter(1) / 200);
                bias_error.push_back(gaussFit4->GetParError(1) / 200);

                sig_yield.push_back(unbinned_comb);
                bias.push_back(gaussFit5->GetParameter(1) / 200);
                bias_error.push_back(gaussFit5->GetParError(1) / 200);

                sig_yield.push_back(unbinned_roofit);
                bias.push_back(gaussFit6->GetParameter(1) / 200);
                bias_error.push_back(gaussFit6->GetParError(1) / 200);

                sig_yield.push_back(bin100_roocomb);
                bias.push_back(gaussFit7->GetParameter(1) / 200);
                bias_error.push_back(gaussFit7->GetParError(1) / 200);

                sig_yield.push_back(bin100_TFF);
                bias.push_back(gaussFit8->GetParameter(1) / 200);
                bias_error.push_back(gaussFit8->GetParError(1) / 200);

                sig_yield.push_back(bin100_SPD);
                bias.push_back(gaussFit9->GetParameter(1) / 200);
                bias_error.push_back(gaussFit9->GetParError(1) / 200);

                sig_yield.push_back(bin100_IWLS);
                bias.push_back(gaussFit10->GetParameter(1) / 200);
                bias_error.push_back(gaussFit10->GetParError(1) / 200);

                all_sig_yield.push_back(sig_yield);
                all_bias.push_back(bias);
                all_bias_error.push_back(bias_error);
                all_fit.push_back({gaussFit1, gaussFit2, gaussFit3, gaussFit4, gaussFit5, gaussFit6, gaussFit7, gaussFit8, gaussFit9, gaussFit10});
                all_pulls.push_back({bin30_roocomb_pull, bin30_TFF_pull, bin30_SPD_pull, bin30_IWLS_pull, unbinned_comb_pull, unbinned_roofit_pull,bin100_roocomb_pull, bin100_TFF_pull, bin100_SPD_pull, bin100_IWLS_pull});
                all_errors.push_back({bin30_roocomb_error, bin30_TFF_error, bin30_SPD_error, bin30_IWLS_error, unbinned_comb_error, unbinned_roofit_error, bin100_roocomb_error, bin100_TFF_error, bin100_SPD_error, bin100_IWLS_error});
    }
    gStyle->SetOptStat(0);
    TCanvas *c2 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c2->SetLogx();
    c2->SetGridy();
    const char *titles[10] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS", "U_RooComb", "U_RooFit", "bin100_roocomb", "bin100_TFF", "bin100_Dembinski", "bin100_IWLS"};
    int colors[10] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kRed+2, kBlue+2, kGreen+2, kMagenta+2};
    std::vector<TGraphErrors*> graphs;

    for (size_t j = 0; j < 6; ++j) {
        TGraphErrors *graph = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            if(i==0){cout<<all_bias[i][j]<<endl;}
            graph->SetPoint(i, size[i], all_bias[i][j]);
            graph->SetPointError(i, 0, all_bias_error[i][j]);
        }

        graph->SetMarkerStyle(20 + j);
        graph->SetMarkerColor(colors[j]);
        graph->SetTitle(titles[j]);

        if (j == 0) {
            graph->GetYaxis()->SetRangeUser(0.98, 1.04);
            graph->Draw("APE");
            graph->SetTitle("Bias vs Template Background N;Template Background N;Bias");
        } else {
            graph->Draw("P E SAME");
        }
        graphs.push_back(graph);
    }

    auto legend = new TLegend(0.5, 0.6, 0.7, 0.9);
    for (size_t j = 0; j < graphs.size(); ++j) {
        legend->AddEntry(graphs[j], titles[j], "p");
    }
    legend->Draw();
    c2->Update();
    c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
    //c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf]");


    for (size_t j = 0; j < size.size(); ++j) {
        TCanvas *c3 = new TCanvas(Form("Nmc%04d",int(j)), "signal_yield", 1200, 500);
        c3->Divide(3, 2);
        for (int i = 1; i <= 6; i++) {
            TPad *pad = (TPad*)c3->GetPad(i);

            // 获取原始位置
            double xlow, ylow, xup, yup;
            pad->GetPadPar(xlow, ylow, xup, yup);

            // 计算缩小后的中心和新尺寸
            double xcenter = (xlow + xup) / 2.0;
            double ycenter = (ylow + yup) / 2.0;
            double xwidth = (xup - xlow) * 0.9; // 宽度缩小到90%
            double yheight = (yup - ylow) * 0.9; // 高度缩小到90%

            // 设置新的Pad位置
            pad->SetPad(xcenter - xwidth / 2.0, ycenter - yheight / 2.0,
                        xcenter + xwidth / 2.0, ycenter + yheight / 2.0);

        }

         TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.03);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.98, Form("Signal Yield for Template Background N = %d", size[j]));
        for(size_t i = 0; i < 6; ++i){
            c3->cd(i+1);
            all_sig_yield[j][i]->SetTitleSize(0.06,"XYT");
            all_sig_yield[j][i]->SetLineColor(colors[i]);
            all_sig_yield[j][i]->GetXaxis()->SetNdivisions(504,"+");
            all_sig_yield[j][i]->GetYaxis()->SetNdivisions(504,"+");
            all_sig_yield[j][i]->GetXaxis()->SetLabelSize(0.06);
            all_sig_yield[j][i]->GetYaxis()->SetLabelSize(0.06);
            all_sig_yield[j][i]->GetXaxis()->SetTitle("");
            all_sig_yield[j][i]->GetYaxis()->SetTitle("Count");
            all_sig_yield[j][i]->GetXaxis()->SetTitleOffset(0.9);
            all_sig_yield[j][i]->GetYaxis()->SetTitleOffset(0.8);


            // all_sig_yield[j][i]->GetXaxis()->SetTitleSize(0.045);
            // all_sig_yield[j][i]->GetYaxis()->SetTitleSize(0.045);
            all_sig_yield[j][i]->Draw();
            all_fit[j][i]->SetLineColor(colors[i]+2);
            all_fit[j][i]->Draw("same");
            TLegend *leg = new TLegend(0.55, 0.5, 0.88, 0.8);
            leg->AddEntry(all_sig_yield[j][i], titles[i], "l");
            leg->AddEntry(
                all_fit[j][i],
                Form("#splitline{mean=%.2f #pm %.2f}{sigma=%.2f}",
                        all_fit[j][i]->GetParameter(1),
                        all_fit[j][i]->GetParError(1),
                        all_fit[j][i]->GetParameter(2)),
                    "l"
                );  
            leg->SetTextSize(0.05);
            leg->SetBorderSize(0);
            leg->Draw();
            TLatex tex;
            tex.SetNDC();               // 使用归一化设备坐标
            tex.SetTextSize(0.06);      // 根据需要调整文字大小
            tex.SetTextAlign(33);       // 右下对齐：3表示右对齐，3表示下对齐
            tex.DrawLatex(0.94, 0.05, "Signal Yield");
             

        }
        c3->Update();
        c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
        // if(j==size.size()-1){
        //     c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf]");
        // }
    }

    for (size_t j = 0; j < size.size(); ++j) {
        TCanvas *c4 = new TCanvas(Form("Signal Yields Error_Nmc%04d", int(j)), "error", 1200, 600);
        c4->Divide(3, 2);
        for (int i = 1; i <= 6; i++) {
            TPad *pad = (TPad*)c4->GetPad(i);

            double xlow, ylow, xup, yup;
            pad->GetPadPar(xlow, ylow, xup, yup);

            double xcenter = (xlow + xup) / 2.0;
            double ycenter = (ylow + yup) / 2.0;
            double xwidth = (xup - xlow) * 0.9;
            double yheight = (yup - ylow) * 0.9;

            pad->SetPad(xcenter - xwidth / 2.0, ycenter - yheight / 2.0,
                        xcenter + xwidth / 2.0, ycenter + yheight / 2.0);
        }

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.SetTextAlign(22);
        latex.DrawLatex(0.5, 0.98, Form("Signal Yield Error for Template Background N = %d", size[j]));

        for (size_t i = 0; i < 6; ++i) {
            c4->cd(i + 1);
            all_errors[j][i]->SetTitleSize(0.05, "XYT");
            all_errors[j][i]->SetLineColor(colors[i]);
            all_errors[j][i]->GetXaxis()->SetNdivisions(504,"+");
            all_errors[j][i]->GetYaxis()->SetNdivisions(504,"+");
            all_errors[j][i]->GetXaxis()->SetLabelSize(0.06);
            all_errors[j][i]->GetYaxis()->SetLabelSize(0.05);
            all_errors[j][i]->GetXaxis()->SetTitle("Error");
            all_errors[j][i]->GetYaxis()->SetTitle("Count");
            all_errors[j][i]->GetYaxis()->SetTitleOffset(1.1);
            // all_errors[j][i]->GetXaxis()->SetTitleSize(0.045);
            // all_errors[j][i]->GetYaxis()->SetTitleSize(0.045);
            all_errors[j][i]->Draw();
        }
        c4->Update();
        c4->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");

        gStyle->SetOptStat(0);
        TCanvas *c5 = new TCanvas(Form("Pull_Nmc%04d", int(j)), "pull", 1200, 600);
        c5->Divide(3, 2);
        for (int i = 1; i <= 6; i++) {
            TPad *pad = (TPad*)c5->GetPad(i);

            double xlow, ylow, xup, yup;
            pad->GetPadPar(xlow, ylow, xup, yup);

            double xcenter = (xlow + xup) / 2.0;
            double ycenter = (ylow + yup) / 2.0;
            double xwidth = (xup - xlow) * 0.9;
            double yheight = (yup - ylow) * 0.9;

            pad->SetPad(xcenter - xwidth / 2.0, ycenter - yheight / 2.0,
                        xcenter + xwidth / 2.0, ycenter + yheight / 2.0);
        }

        latex.DrawLatex(0.5, 0.98, Form("Pull for Template Background N = %d", size[j]));

        for (size_t i = 0; i < 6; ++i) {
            c5->cd(i + 1);
            all_pulls[j][i]->SetTitleSize(0.05, "XYT");
            all_pulls[j][i]->SetLineColor(colors[i]);
            all_pulls[j][i]->GetXaxis()->SetNdivisions(504, kTRUE);
            all_pulls[j][i]->GetYaxis()->SetNdivisions(505, kTRUE);
            all_pulls[j][i]->GetXaxis()->SetLabelSize(0.06);
            all_pulls[j][i]->GetYaxis()->SetLabelSize(0.06);
            all_pulls[j][i]->GetXaxis()->SetTitle("Pull");
            all_pulls[j][i]->GetYaxis()->SetTitle("");
            // all_pulls[j][i]->GetXaxis()->SetTitleSize(0.045);
            // all_pulls[j][i]->GetYaxis()->SetTitleSize(0.045);
            all_pulls[j][i]->Draw();
            TLegend *leg = new TLegend(0.6, 0.6, 0.88, 0.8);
            leg->AddEntry(all_pulls[j][i], titles[i], "l");
            leg->AddEntry((TObject*)0, Form("mean=%.2f", all_pulls[j][i]->GetMean()), "");
            leg->AddEntry((TObject*)0, Form("std dev=%.2f", all_pulls[j][i]->GetStdDev()), "");
            leg->SetTextSize(0.045);
            leg->SetBorderSize(0);
            leg->Draw("same");
            TLatex tex;
            tex.SetTextAngle(90);
            tex.SetTextSize(0.05);  // 根据需要调整文字大小
            tex.SetNDC();               // 使用归一化设备坐标
            tex.SetTextAlign(13);       // 左上角对齐
            // x, y 坐标需根据画布实际情况进行调整
            tex.DrawLatex(0.03, 0.85, "Count");
        }
        c5->Update();
        c5->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
    }


    TCanvas *c6 = new TCanvas("c6", "Signal Yield Comparison", 1200, 600);
    c6->Divide(3, 1);
     gStyle->SetOptStat(0);
            for (int i = 1; i <= 3; i++) {
            TPad *pad = (TPad*)c6->GetPad(i);

            double xlow, ylow, xup, yup;
            pad->GetPadPar(xlow, ylow, xup, yup);

            double xcenter = (xlow + xup) / 2.0;
            double ycenter = (ylow + yup) / 2.0;
            double xwidth = (xup - xlow) * 1;
            double yheight = (yup - ylow) * 1;

            pad->SetPad(xcenter - xwidth /2, ycenter - yheight / 2.0,
                        xcenter + xwidth / 2, ycenter + yheight / 2.0);
            pad->SetLeftMargin(0.15); // 增加左边距
            pad->SetRightMargin(0.05); // 适度缩小右边距
        }
    std::vector<int> selected_sizes = {100, 1500, 5000};
    for (size_t k = 0; k < selected_sizes.size(); ++k) {
    c6->cd(k + 1);
    TH1F *h1 = new TH1F(Form("h1_%zu", k), 
                        Form("Signal Yield for template background N = %d", selected_sizes[k]), 
                        6, 0.5, 6.5); // 6 bins, X 轴从 0.5 到 6.5
    h1->SetTitleSize(0.05, "XYT"); // 设置标题字体大小
    //h1->SetTitleSize(0.1, "T");
    int index = 0;
    if (k == 0) {
        index = 0;
    } else if (k == 1) {
        index = 1;
    } else if (k == 2) {
        index = 8;
    }

    // 设置每个 bin 的内容和误差
    for (size_t i = 0; i < 6; ++i) {
        int k = 0;
        if(i<3){
            k=i+1;}
        else if(i==3){
            k=0;}
        else{k=i;}
        h1->SetBinContent(i + 1, all_fit[index][k]->GetParameter(1));
        h1->SetBinError(i + 1, all_fit[index][k]->GetParError(1));
    }

    // 设置 X 轴的标签
    TAxis *xAxis = h1->GetXaxis();
    xAxis->SetBinLabel(1, "TFF");
    xAxis->SetBinLabel(2, "ALM-TF");
    xAxis->SetBinLabel(3, "IWLS");
    xAxis->SetBinLabel(4, "RooCombineFit");
    xAxis->SetBinLabel(5, "U_RooComb");
    xAxis->SetBinLabel(6, "U_Roofit");
        // 设置标签字体大小
    xAxis->SetLabelSize(0.04); // 字体大小

    // 设置标签字体粗细
    xAxis->SetLabelFont(62); // 42 表示正常字体，62 表示粗体
    //xAxis->SetLabelOffset(-0.02); // 负值表示左移，正值表示右移
    // 旋转标签
    //xAxis->LabelsOption("v");

    // 设置 Y 轴范围和标题
    h1->GetYaxis()->SetTitle("Signal Yield");
    h1->GetYaxis()->SetRangeUser(199, 206.5);

    // 设置 X 轴标题
    xAxis->SetTitle("Method");
    xAxis->SetTitleOffset(1);

     TAxis *yAxis = h1->GetYaxis();
    yAxis->SetLabelSize(0.04); // 字体大小
    yAxis->SetLabelFont(62); // 42 表示正常字体，62 表示粗体
    yAxis->SetNdivisions(505, kTRUE); // 设置 Y 轴刻度数量

    // 设置样式
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.8);
    h1->SetMarkerColor(kBlack);

    // 绘制图表
    gStyle->SetErrorX(0); // 不显示 X 轴误差
    h1->SetLineColor(kBlack); // 设置误差棒颜色为黑色
    h1->Draw("E"); // 带误差条的点图

    // 添加水平虚线 y=200
    TLine *line = new TLine(0.5, 200, 6.5, 200);
    line->SetLineStyle(2); // Dashed line
    line->SetLineColor(kBlack);
    line->Draw("same");
}
    c6->Update();
    c6->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
    

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TCanvas *c7 = new TCanvas("c7", "Signal Yields Sigma vs Nmc", 800, 600);
    c7->SetLogx();
    c7->SetGridy();
    const char *titles1[10] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS", "U_RooComb", "U_RooFit", "bin100_roocomb", "bin100_TFF", "bin100_Dembinski", "bin100_IWLS"};
    int colors1[10] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kRed+2, kBlue+2, kGreen+2, kMagenta+2};
    std::vector<TGraphErrors*> graphs1;

    for (size_t j = 0; j < 6; ++j) {
        TGraphErrors *graph1 = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            graph1->SetPoint(i, size[i], all_fit[i][j]->GetParameter(2));
            graph1->SetPointError(i, 0, all_fit[i][j]->GetParError(2));
            //cout<<all_bias_error[j][i]<<endl;
        }

        graph1->SetMarkerStyle(20 + j);
        graph1->SetMarkerColor(colors1[j]);
        graph1->SetTitle(titles1[j]);

        if (j == 0) {
            graph1->GetYaxis()->SetRangeUser(6, 9.5);
            graph1->Draw("APE");
            graph1->SetTitle("Signal Yields Sigma vs Template Background N;Template Background N;Signal Yields Sigma");
        } else {
            graph1->Draw("P E SAME");
        }
        graphs1.push_back(graph1);
    }

    auto legend1 = new TLegend(0.5, 0.6, 0.7, 0.9);
    for (size_t j = 0; j < graphs1.size(); ++j) {
        legend1->AddEntry(graphs1[j], titles1[j], "p");
    }
    legend1->Draw();
    c7->Update();
    c7->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
    
    

    
    
    c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf]");
        
    

}