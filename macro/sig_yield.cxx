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

void sig_yield(){
    std::vector<int> size = {100, 150, 200, 500, 700, 1000, 1500, 2000, 5000, 8000};
    std::vector<std::vector<TH1D*>> all_sig_yield;
    std::vector<std::vector<double>> all_bias;
    std::vector<std::vector<double>> all_bias_error;
    std::vector<std::vector<TF1*>> all_fit;
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    c1->Divide(3,4);
    c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf[");

    for (int i : size){ 
            // gStyle->SetLabelSize(0.1, "X");
            std::vector<TH1D*> sig_yield;
            std::vector<double> bias;
            std::vector<double> bias_error;
            TFile *f=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_roocomb_seq0001_npi%04d.root",i),"read");
            TFile *f1=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_TFF_seq0001_npi%04d.root",i),"read");
            //TFile *f2=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_SPD_seq0001_npi%04d.root",i),"read");
            TFile *f2=new TFile(Form("/data03/tianye/pbar/root/paper_data/binned_SPD_seq0001_npi%04d.root",i),"read");
            TFile *f3=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/binned_IWLS1_seq0001_npi%04d.root",i),"read");
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

            t->SetBranchAddress("fitparams_bin25", &params1);
            t->SetBranchAddress("fitparams_bin100", &params2);
            t1->SetBranchAddress("nsig", &params3);
            t2->SetBranchAddress("bin25_SPD_pararr", &params4);
            t2->SetBranchAddress("bin100_SPD_pararr", &params5);
            t3->SetBranchAddress("nsig", &params6);
            t4->SetBranchAddress("fitparams_unbinned", &params7);
            t5->SetBranchAddress("fitparams_unbinned", &params8);

            TH1D* bin20_roocomb = new TH1D("bin20_roocomb", "bin20_roocomb", 200, 140, 280);
            TH1D* bin100_roocomb = new TH1D("bin100_roocomb", "bin100_roocomb", 200, 140, 280);
            TH1D* bin20_TFF = new TH1D("bin20_TFF", "bin20_TFF", 200, 140, 280);
            TH1D* bin100_TFF = new TH1D("bin100_TFF", "bin100_TFF", 200, 140, 280);
            TH1D* bin20_SPD = new TH1D("bin20_SPD", "bin20_paper", 200, 140, 280);
            TH1D* bin100_SPD = new TH1D("bin100_SPD", "bin100_paper", 200, 140, 280);
            TH1D* bin20_IWLS = new TH1D("bin20_IWLS", "bin20_IWLS", 200, 140, 280);
            TH1D* bin100_IWLS = new TH1D("bin100_IWLS", "bin100_IWLS", 200, 140, 280);
            TH1D* unbinned_comb = new TH1D("unbinned_comb", "unbinned_comb", 200, 140, 280);
            TH1D* unbinned_roofit = new TH1D("unbinned_roofit", "unbinned_roofit", 200, 140, 280);

            TH1D* bin20_roocomb_error= new TH1D("bin20_roocomb_error", "bin20_roocomb_error", 200, 10, 30);
            TH1D* bin100_roocomb_error= new TH1D("bin100_roocomb_error", "bin100_roocomb_error", 200, 10, 30);
            TH1D* bin20_TFF_error= new TH1D("bin20_TFF_error", "bin20_TFF_error", 200, 120, 160);
            TH1D* bin100_TFF_error= new TH1D("bin100_TFF_error", "bin100_TFF_error", 200, 50, 100);
            TH1D* bin20_SPD_error= new TH1D("bin20_paper_error", "bin20_paper_error", 200, 20, 30);
            TH1D* bin100_SPD_error= new TH1D("bin100_paper_error", "bin100_paper_error", 200, 20, 30);
            TH1D* bin20_IWLS_error= new TH1D("bin20_IWLS_error", "bin20_IWLS_error", 200, 50, 110);
            TH1D* bin100_IWLS_error= new TH1D("bin100_IWLS_error", "bin100_IWLS_error", 200, 35, 55);
            TH1D* unbinned_comb_error= new TH1D("unbinned_comb_error", "unbinned_comb_error", 200, 10, 30);
            TH1D* unbinned_roofit_error= new TH1D("unbinned_roofit_error", "unbinned_roofit_error", 200, 10, 30);


            TH1D* bin20_roocomb_pull= new TH1D("bin20_roocomb_pull", "bin20_roocomb_pull", 200, -3, 3);
            TH1D* bin100_roocomb_pull= new TH1D("bin100_roocomb_pull", "bin100_roocomb_pull", 200, -3, 3);
            TH1D* bin20_TFF_pull= new TH1D("bin20_TFF_pull", "bin20_TFF_pull", 200, -3, 3);
            TH1D* bin100_TFF_pull= new TH1D("bin100_TFF_pull", "bin100_TFF_pull", 200, -3, 3);
            TH1D* bin20_SPD_pull= new TH1D("bin20_SPD_pull", "bin20_paper_pull", 200, -3, 3);
            TH1D* bin100_SPD_pull= new TH1D("bin100_SPD_pull", "bin100_paper_pull", 200, -3, 3);
            TH1D* bin20_IWLS_pull= new TH1D("bin20_IWLS_pull", "bin20_IWLS_pull", 200, -3, 3);
            TH1D* bin100_IWLS_pull= new TH1D("bin100_IWLS_pull", "bin100_IWLS_pull", 200, -3, 3);
            TH1D* unbinned_comb_pull= new TH1D("unbinned_comb_pull", "unbinned_comb_pull", 200, -3, 3);
            TH1D* unbinned_roofit_pull= new TH1D("unbinned_roofit_pull", "unbinned_roofit_pull", 200, -3, 3);


                int nEntries = t->GetEntries();
                for (int i = 0; i < nEntries; ++i) {
                    t->GetEntry(i);
                    t1->GetEntry(i);
                    t2->GetEntry(i);
                    t3->GetEntry(i);
                    t4->GetEntry(i);
                    t5->GetEntry(i);

                    bin20_roocomb->Fill((*params1)[10] * 400);
                    bin100_roocomb->Fill((*params2)[10] * 400);
                    bin20_TFF->Fill((*params3)[0] * 400);
                    bin100_TFF->Fill((*params3)[1] * 400);
                    // bin20_SPD->Fill((*params4)[0] * 400);
                    // bin100_SPD->Fill((*params5)[0] * 400);
                    bin20_SPD->Fill(params4[0] );
                    bin100_SPD->Fill(params5[0] );
                    bin20_IWLS->Fill((*params6)[0] * 400);
                    bin100_IWLS->Fill((*params6)[1] * 400);
                    unbinned_comb->Fill((*params7)[10] * 400);
                    unbinned_roofit->Fill((*params8)[0]);


                    // bin20_roocomb_error->Fill((*params1)[21] * 400);
                    // bin100_roocomb_error->Fill((*params2)[21] * 400);
                    // bin20_TFF_error->Fill((*params3)[2] * 400);
                    // bin100_TFF_error->Fill((*params3)[3] * 400);
                    // bin20_SPD_error->Fill((*params4)[2] * 400);
                    // bin100_SPD_error->Fill((*params5)[2] * 400);
                    // bin20_IWLS_error->Fill((*params6)[2] * 400);
                    // bin100_IWLS_error->Fill((*params6)[3] * 400);
                    // unbinned_comb_error->Fill((*params7)[21] * 400);
                    // unbinned_roofit_error->Fill((*params8)[1]);

                    // bin20_roocomb_pull->Fill(((*params1)[10] - 0.5)/(*params1)[21] );
                    // bin100_roocomb_pull->Fill(((*params2)[10] - 0.5)/(*params2)[21] );
                    // bin20_TFF_pull->Fill(((*params3)[0] - 0.5)/(*params3)[2] );
                    // bin100_TFF_pull->Fill(((*params3)[1] - 0.5)/(*params3)[3] );
                    // bin20_SPD_pull->Fill(((*params4)[0] - 0.5)/(*params4)[2] );
                    // bin100_SPD_pull->Fill(((*params5)[0] - 0.5)/(*params5)[2] );
                    // bin20_IWLS_pull->Fill(((*params6)[0] - 0.5)/(*params6)[2] );
                    // bin100_IWLS_pull->Fill(((*params6)[1] - 0.5)/(*params6)[3] );
                    // unbinned_comb_pull->Fill(((*params7)[10] - 0.5)/(*params7)[21] );
                    // unbinned_roofit_pull->Fill(((*params8)[0] - 200)/(*params8)[1] );
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

                bin20_roocomb->Fit(gaussFit1, "Q");
                bin100_roocomb->Fit(gaussFit2, "Q");
                bin20_TFF->Fit(gaussFit3, "Q");
                bin100_TFF->Fit(gaussFit4, "Q");
                bin20_SPD->Fit(gaussFit5, "Q");
                bin100_SPD->Fit(gaussFit6, "Q");
                bin20_IWLS->Fit(gaussFit7, "Q");
                bin100_IWLS->Fit(gaussFit8, "Q");
                unbinned_comb->Fit(gaussFit9, "Q");
                unbinned_roofit->Fit(gaussFit10, "Q");

                sig_yield.push_back(bin20_roocomb);
                bias.push_back(gaussFit1->GetParameter(1) / 200);
                bias_error.push_back(gaussFit1->GetParError(1) / 200);
                
                sig_yield.push_back(bin100_roocomb);
                bias.push_back(gaussFit2->GetParameter(1) / 200);
                bias_error.push_back(gaussFit2->GetParError(1) / 200);
                
                sig_yield.push_back(bin20_TFF);
                bias.push_back(gaussFit3->GetParameter(1) / 200);
                bias_error.push_back(gaussFit3->GetParError(1) / 200);
                
                sig_yield.push_back(bin100_TFF);
                bias.push_back(gaussFit4->GetParameter(1) / 200);
                bias_error.push_back(gaussFit4->GetParError(1) / 200);
                
                sig_yield.push_back(bin20_SPD);
                bias.push_back(gaussFit5->GetParameter(1) / 200);
                bias_error.push_back(gaussFit5->GetParError(1) / 200);
                
                sig_yield.push_back(bin100_SPD);
                bias.push_back(gaussFit6->GetParameter(1) / 200);
                bias_error.push_back(gaussFit6->GetParError(1) / 200);
                
                sig_yield.push_back(bin20_IWLS);
                bias.push_back(gaussFit7->GetParameter(1) / 200);
                bias_error.push_back(gaussFit7->GetParError(1) / 200);
                
                sig_yield.push_back(bin100_IWLS);
                bias.push_back(gaussFit8->GetParameter(1) / 200);
                bias_error.push_back(gaussFit8->GetParError(1) / 200);
                
                sig_yield.push_back(unbinned_comb);
                bias.push_back(gaussFit9->GetParameter(1) / 200);
                bias_error.push_back(gaussFit9->GetParError(1) / 200);
                
                sig_yield.push_back(unbinned_roofit);
                bias.push_back(gaussFit10->GetParameter(1) / 200);
                bias_error.push_back(gaussFit10->GetParError(1) / 200);

                all_sig_yield.push_back(sig_yield);
                all_bias.push_back(bias);
                all_bias_error.push_back(bias_error);
                all_fit.push_back({gaussFit1, gaussFit2, gaussFit3, gaussFit4, gaussFit5, gaussFit6, gaussFit7, gaussFit8, gaussFit9, gaussFit10});
    
    }

    TCanvas *c2 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c2->SetLogx();
    c2->SetGridy();
    const char *titles[10] = {"bin20_roocomb", "bin100_roocomb", "bin20_TFF", "bin100_TFF", "bin20_SPD", "bin100_SPD", "bin20_IWLS", "bin100_IWLS", "unbinned_comb", "unbinned_roofit"};
    int colors[10] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kRed+2, kBlue+2, kGreen+2, kMagenta+2};
    std::vector<TGraphErrors*> graphs;

    for (size_t j = 0; j < 10; ++j) {
        TGraphErrors *graph = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            graph->SetPoint(i, size[i], all_bias[i][j]);
            graph->SetPointError(i, 0, all_bias_error[j][i]);
            cout<<all_bias_error[j][i]<<endl;
        }

        graph->SetMarkerStyle(20 + j);
        graph->SetMarkerColor(colors[j]);
        graph->SetTitle(titles[j]);

        if (j == 0) {
            graph->GetYaxis()->SetRangeUser(0.98, 1.12);
            graph->Draw("APE");
            graph->SetTitle("Bias vs Nmc;Nmc;Bias");
        } else {
            graph->Draw("P E SAME");
        }
        graphs.push_back(graph);
    }

    auto legend = new TLegend(0.5, 0.5, 0.7, 0.9);
    for (size_t j = 0; j < graphs.size(); ++j) {
        legend->AddEntry(graphs[j], titles[j], "p");
    }
    legend->Draw();
    c2->Update();
    c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
    //c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf]");


    for (size_t j = 0; j < size.size(); ++j) {
        TCanvas *c3 = new TCanvas(Form("Nmc%04d",int(j)), "signal_yield", 800, 600);
        c3->Divide(2, 5);
        c3->SetTitle(Form("Signal Yield for pi_Nmc = %d", size[j]));
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.92, c3->GetTitle());
        for(size_t i = 0; i < 10; ++i){
            c3->cd(i+1);
            all_sig_yield[j][i]->Draw();
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.05);
            latex.DrawLatex(0.65, 0.85, Form("Mean: %.2f", all_fit[j][i]->GetParameter(1)));
            latex.DrawLatex(0.65, 0.80, Form("Sigma: %.2f", all_fit[j][i]->GetParameter(2)));
            latex.DrawLatex(0.65, 0.75, Form("Mean Error: %.2f", all_fit[j][i]->GetParError(1)));
            all_fit[j][i]->Draw("SAME");
        }
        c3->Update();
        c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf");
        if(j==size.size()-1){
            c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/bias_plot8000.pdf]");
        }
    }
        
    

}