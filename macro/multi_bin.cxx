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

void multi_bin(){
    std::vector<int> size = {100, 150, 200, 500, 700, 1000, 1500, 2000, 5000, 8000};
    std::vector<std::vector<TH1D*>> all_sig_yield;
    std::vector<std::vector<TH1D*>> all_errors;
    std::vector<std::vector<double>> all_bias;
    std::vector<std::vector<double>> all_bias_error;
    std::vector<std::vector<TF1*>> all_fit;
    std::vector<std::vector<TH1D*>> all_pulls;
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);

    c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_multi_bin.pdf[");

for (int i : size){ 
            // gStyle->SetLabelSize(0.1, "X");
            std::vector<TH1D*> sig_yield;
            std::vector<TH1D*> sig_yield_error;
            std::vector<TH1D*> sig_yield_pull;
            std::vector<double> bias;
            std::vector<double> bias_error;

            TFile *f=new TFile(Form("/data03/tianye/pbar/root/multi_bin/binned_roocomb_seq0001_npi%04d.root",i),"read");

            std::vector<double>* params1 = new std::vector<double>(22);
            std::vector<double>* params2 = new std::vector<double>(22);
            std::vector<double>* params3 = new std::vector<double>(22);
            std::vector<double>* params4 = new std::vector<double>(22);
            std::vector<double>* params5 = new std::vector<double>(22);

            TTree *t = (TTree*)f->Get("fitParamsTree");


            t->SetBranchAddress("fitparams_bin20", &params1);
            t->SetBranchAddress("fitparams_bin30", &params2);  
            t->SetBranchAddress("fitparams_bin40", &params3);
            t->SetBranchAddress("fitparams_bin50", &params4);          
            t->SetBranchAddress("fitparams_bin100", &params5);


            TH1D* bin20_roocomb = new TH1D("bin20_roocomb", "bin20_roocomb", 200, 140, 280);
            TH1D* bin30_roocomb = new TH1D("bin30_roocomb", "bin30_roocomb", 200, 140, 280);
            TH1D* bin40_roocomb = new TH1D("bin40_roocomb", "bin40_roocomb", 200, 140, 280);
            TH1D* bin50_roocomb = new TH1D("bin50_roocomb", "bin50_roocomb", 200, 140, 280);
            TH1D* bin100_roocomb = new TH1D("bin100_roocomb", "bin100_roocomb", 200, 140, 280);

                int nEntries = t->GetEntries();
                for (int i = 0; i < nEntries; ++i) {
                    t->GetEntry(i);
                   

                    bin20_roocomb->Fill((*params1)[10] * 400);
                    bin30_roocomb->Fill((*params2)[10] * 400);
                    bin40_roocomb->Fill((*params3)[10] * 400);
                    bin50_roocomb->Fill((*params4)[10] * 400);
                    bin100_roocomb->Fill((*params5)[10] * 400);
                  
                }

                // Fit histograms with Gaussian
                TF1* gaussFit1 = new TF1("gaussFit1", "gaus", 100, 300);
                TF1* gaussFit2 = new TF1("gaussFit2", "gaus", 100, 300);
                TF1* gaussFit3 = new TF1("gaussFit3", "gaus", 100, 300);
                TF1* gaussFit4 = new TF1("gaussFit4", "gaus", 100, 300);
                TF1* gaussFit5 = new TF1("gaussFit5", "gaus", 100, 300);

                bin20_roocomb->Fit(gaussFit1, "Q");
                bin30_roocomb->Fit(gaussFit2, "Q");
                bin40_roocomb->Fit(gaussFit3, "Q");
                bin50_roocomb->Fit(gaussFit4, "Q");
                bin100_roocomb->Fit(gaussFit5, "Q");


                sig_yield.push_back(bin20_roocomb);
                bias.push_back(gaussFit1->GetParameter(1) / 200);
                bias_error.push_back(gaussFit1->GetParError(1) / 200);

                sig_yield.push_back(bin30_roocomb);
                bias.push_back(gaussFit2->GetParameter(1) / 200);
                bias_error.push_back(gaussFit2->GetParError(1) / 200);

                sig_yield.push_back(bin40_roocomb);
                bias.push_back(gaussFit3->GetParameter(1) / 200);
                bias_error.push_back(gaussFit3->GetParError(1) / 200);

                sig_yield.push_back(bin50_roocomb);
                bias.push_back(gaussFit4->GetParameter(1) / 200);
                bias_error.push_back(gaussFit4->GetParError(1) / 200);

                sig_yield.push_back(bin100_roocomb);
                bias.push_back(gaussFit5->GetParameter(1) / 200);
                bias_error.push_back(gaussFit5->GetParError(1) / 200);

            

                all_sig_yield.push_back(sig_yield);
                all_bias.push_back(bias);
                all_bias_error.push_back(bias_error);
                all_fit.push_back({gaussFit1, gaussFit2, gaussFit3, gaussFit4, gaussFit5});
               
    }

    TCanvas *c2 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c2->SetLogx();
    c2->SetGridy();
    const char *titles[5] = {"bin20_roocomb", "bin30_roocomb", "bin40_roocomb", "bin50_roocomb", "bin100_roocomb"};
    int colors[5] = {kRed, kBlue, kGreen, kMagenta, kCyan};
    std::vector<TGraphErrors*> graphs;

    for (size_t j = 0; j < 5; ++j) {
        TGraphErrors *graph = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            graph->SetPoint(i, size[i], all_bias[i][j]);
            graph->SetPointError(i, 0, all_bias_error[i][j]);
            cout << all_bias_error[i][j] << endl;
        }

        graph->SetMarkerStyle(20 + j);
        graph->SetMarkerColor(colors[j]);
        graph->SetTitle(titles[j]);

        if (j == 0) {
            graph->GetYaxis()->SetRangeUser(0.98, 1.12);
            graph->Draw("APE");
            graph->SetTitle("Bias vs N_template_bkg;N_template_bkg;Bias");
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
    c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_multi_bin.pdf");


    for (size_t j = 0; j < size.size(); ++j) {
        TCanvas *c3 = new TCanvas(Form("Nmc%04d", int(j)), "signal_yield", 1200, 300);
        c3->Divide(3, 2);
        for (int i = 1; i <= 5; i++) {
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
        latex.DrawLatex(0.5, 0.98, Form("Signal Yield for N_template_background = %d", size[j]));

        for (size_t i = 0; i < 5; ++i) {
            c3->cd(i + 1);
            all_sig_yield[j][i]->SetLineColor(colors[i]);
            all_sig_yield[j][i]->GetXaxis()->SetNdivisions(504, kTRUE);
            all_sig_yield[j][i]->GetYaxis()->SetNdivisions(505, kTRUE);
            all_sig_yield[j][i]->GetXaxis()->SetLabelSize(0.06);
            all_sig_yield[j][i]->GetYaxis()->SetLabelSize(0.06);
            all_sig_yield[j][i]->GetXaxis()->SetTitle("Signal Yield");
            all_sig_yield[j][i]->GetYaxis()->SetTitle("Count");
            all_sig_yield[j][i]->GetYaxis()->SetTitleOffset(1.2);
            all_sig_yield[j][i]->GetXaxis()->SetTitleSize(0.045);
            all_sig_yield[j][i]->GetYaxis()->SetTitleSize(0.045);
            all_sig_yield[j][i]->Draw();
            all_fit[j][i]->SetLineColor(colors[i] + 2);
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
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->Draw();
        }
        c3->Update();
        c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_multi_bin.pdf");
    }
    c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_multi_bin.pdf]");
}