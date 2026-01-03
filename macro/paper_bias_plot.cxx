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

void paper_bias_plot(){
    std::vector<int> size = {50,75,100, 150, 200,300,400,500, 700, 1000, 1500, 2000,3000,4000,5000, 8000};
    int n=size.size();
    double sig_yield[3][7][n];
    double sig_yield_error[3][7][n]; 
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(0);
    c1->Divide(3,4);
    c1->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf[");

    for (int i=0; i < size.size(); ++i){ 



            TFile *f=new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_roocomb_seq0001_npi%04d.root",size[i]),"read");
            TFile *f1=new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_TFF_seq0001_npi%04d.root",size[i]),"read");
            TFile *f2=new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_SPD_seq0001_npi%04d.root",size[i]),"read");
            TFile *f3=new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_IWLS1_seq0001_npi%04d.root",size[i]),"read");
            // TFile *f4=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/unbinned_roocomb_seq0001_npi%04d.root",i),"read");
            // TFile *f5=new TFile(Form("/data03/tianye/pbar/root/bin25_Nmc8000/unbinned_rookeyspdf1_seq0001_npi%04d.root",i),"read");

            std::vector<double>* params1 = new std::vector<double>(12);
            std::vector<double>* params2 = new std::vector<double>(12);
            std::vector<double>* params3 = new std::vector<double>(12);
            std::vector<double>* params4 = new std::vector<double>(12);
            std::vector<double>* params5 = new std::vector<double>(12);

            std::vector<double>* params6 = new std::vector<double>(4);
            std::vector<double>* params7 = new std::vector<double>(22);
            std::vector<double>* params8 = new std::vector<double>(2);

            TTree *t = (TTree*)f->Get("fitParamsTree");
            TTree *t1 = (TTree*)f1->Get("fitParamsTree");
            TTree *t2 = (TTree*)f2->Get("fitParamsTree");
            TTree* t3 = (TTree*)f3->Get("fitParamsTree");
            // TTree* t4 = (TTree*)f4->Get("fitParamsTree");
            // TTree* t5 = (TTree*)f5->Get("fitParamsTree");
        
            t->SetBranchAddress("nsig", &params1);
            t1->SetBranchAddress("nsig", &params2);
            t2->SetBranchAddress("nsig", &params3);
            t3->SetBranchAddress("nsig", &params4);
        for(int j = 0; j < 3; ++j) {
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



                int nEntries = t->GetEntries();
                for (int k = 0; k < nEntries; ++k) {
                    t->GetEntry(k);
                    t1->GetEntry(k);
                    t2->GetEntry(k);
                    t3->GetEntry(k);
                    t4->GetEntry(k);
                    t5->GetEntry(k);

                    bin30_roocomb->Fill((*params1)[j] * 400);
                    bin30_TFF->Fill((*params2)[j] * 400);
                    bin30_SPD->Fill((*params3)[j] * 400);
                    bin30_IWLS->Fill((*params4)[j] * 400);

                }

                // Fit histograms with Gaussian
                TF1* gaussFit1 = new TF1("gaussFit1", "gaus", 100, 300);
                TF1* gaussFit2 = new TF1("gaussFit2", "gaus", 100, 300);
                TF1* gaussFit3 = new TF1("gaussFit3", "gaus", 100, 300);
                TF1* gaussFit4 = new TF1("gaussFit4", "gaus", 100, 300);
                // TF1* gaussFit5 = new TF1("gaussFit5", "gaus", 100, 300);
                // TF1* gaussFit6 = new TF1("gaussFit6", "gaus", 100, 300);
                // TF1* gaussFit7 = new TF1("gaussFit7", "gaus", 100, 300);
                // TF1* gaussFit8 = new TF1("gaussFit8", "gaus", 100, 300);
                // TF1* gaussFit9 = new TF1("gaussFit9", "gaus", 100, 300);
                // TF1* gaussFit10 = new TF1("gaussFit10", "gaus", 100, 300);

                bin30_roocomb->Fit(gaussFit1, "Q");
                bin30_TFF->Fit(gaussFit2, "Q");
                bin30_SPD->Fit(gaussFit3, "Q");
                bin30_IWLS->Fit(gaussFit4, "Q");

                sig_yield[j][0][i] = gaussFit1->GetParameter(1) / 200;
                sig_yield_error[j][0][i] = gaussFit1->GetParError(1) / 200;

                sig_yield[j][1][i] = gaussFit2->GetParameter(1) / 200;
                sig_yield_error[j][1][i] = gaussFit2->GetParError(1) / 200;
                
                sig_yield[j][2][i] = gaussFit3->GetParameter(1) / 200;
                sig_yield_error[j][2][i] = gaussFit3->GetParError(1) / 200;
                
                sig_yield[j][3][i] = gaussFit4->GetParameter(1) / 200;
                sig_yield_error[j][3][i] = gaussFit4->GetParError(1) / 200;
    
            }
        }
            gStyle->SetOptStat(0);
    TCanvas *c2 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c2->SetLogx();
    c2->SetGridy();
    const char *titles[11] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS", "U_RooComb", "U_RooFit", "U_RooKeysPdf","bin100_roocomb", "bin100_TFF", "bin100_Dembinski", "bin100_IWLS"};
    int colors[11] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange,kViolet, kRed+2, kBlue+2, kGreen+2, kMagenta+2};
    std::vector<TGraphErrors*> graphs;

    for (size_t j = 0; j < 7; ++j) {
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
            graph->GetYaxis()->SetRangeUser(0.98, 1.06);
            graph->Draw("APE");
            graph->SetTitle("Bias vs Template Background N;Template Background N;Bias");
            graph->GetXaxis()->SetTitleSize(0.05);
            graph->GetXaxis()->SetTitleOffset(0.85);
            graph->GetYaxis()->SetTitleSize(0.05);
            graph->GetYaxis()->SetTitleOffset(0.9);
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
    c2->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");


    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TCanvas *c7 = new TCanvas("c7", "Signal Yields Sigma vs Nmc", 800, 600);
    c7->SetLogx();
    c7->SetGridy();
    const char *titles1[11] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS", "U_RooComb", "U_RooFit","U_RooKeysPdf", "bin100_roocomb", "bin100_TFF", "bin100_Dembinski", "bin100_IWLS"};
    int colors1[11] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange,kViolet, kRed+2, kBlue+2, kGreen+2, kMagenta+2};
    std::vector<TGraphErrors*> graphs1;

    for (size_t j = 0; j < 7; ++j) {
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
            graph1->GetYaxis()->SetRangeUser(6, 12);
            graph1->Draw("APE");
            graph1->SetTitle("Signal Yields Sigma vs Template Background N;Template Background N;Signal Yields Sigma");
            graph1->GetXaxis()->SetTitleSize(0.05);
            graph1->GetXaxis()->SetTitleOffset(0.85);
            graph1->GetYaxis()->SetTitleSize(0.05);
            graph1->GetYaxis()->SetTitleOffset(0.9);
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
    c7->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");





    TCanvas *c8 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c8->SetLogx();
    c8->SetGridy();
    const char *titles2[4] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS" };
    int colors2[4] = {kRed, kBlue, kGreen, kMagenta};
    std::vector<TGraphErrors*> graphs2;

    for (size_t j = 0; j < 4; ++j) {
        TGraphErrors *graph2 = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            if(i==0){cout<<all_bias[i][j]<<endl;}
            graph2->SetPoint(i, size[i], all_bias[i][j]);
            graph2->SetPointError(i, 0, all_bias_error[i][j]);
        }

        graph2->SetMarkerStyle(20 + j);
        graph2->SetMarkerColor(colors2[j]);
        graph2->SetTitle(titles2[j]);

        if (j == 0) {
            graph2->GetYaxis()->SetRangeUser(0.98, 1.15);
            graph2->Draw("APE");
            graph2->SetTitle("Bias vs Template Background N;Template Background N;Bias");
            graph2->GetXaxis()->SetTitleSize(0.05);
            graph2->GetXaxis()->SetTitleOffset(0.85);
            graph2->GetYaxis()->SetTitleSize(0.05);
            graph2->GetYaxis()->SetTitleOffset(0.9);
        } else {
            graph2->Draw("P E SAME");
        }
        graphs2.push_back(graph2);
    }

    auto legend2 = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs2.size(); ++j) {
        if(j==graphs2.size()-1){legend2->AddEntry(graphs2[0], titles2[0], "p");continue;}
        legend2->AddEntry(graphs[j+1], titles[j+1], "p");
    }
    legend2->Draw();
    c8->Update();
    c8->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");


    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TCanvas *c9 = new TCanvas("c9", "Signal Yields Sigma vs Nmc", 800, 600);
    c9->SetLogx();
    c9->SetGridy();
    const char *titles3[4] = {"RooCombineFit", "TFF", "ALM-TF", "IWLS",};
    int colors3[4] = {kRed, kBlue, kGreen, kMagenta};
    std::vector<TGraphErrors*> graphs3;

    for (size_t j = 0; j < 4; ++j) {
        TGraphErrors *graph3 = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            graph3->SetPoint(i, size[i], all_fit[i][j]->GetParameter(2));
            graph3->SetPointError(i, 0, all_fit[i][j]->GetParError(2));
            //cout<<all_bias_error[j][i]<<endl;
        }

        graph3->SetMarkerStyle(20 + j);
        graph3->SetMarkerColor(colors1[j]);
        graph3->SetTitle(titles1[j]);

        if (j == 0) {
            graph3->GetYaxis()->SetRangeUser(6, 12);
            graph3->Draw("APE");
            graph3->SetTitle("Signal Yields Sigma vs Template Background N;Template Background N;Signal Yields Sigma");
            graph3->GetXaxis()->SetTitleSize(0.05);
            graph3->GetXaxis()->SetTitleOffset(0.85);
            graph3->GetYaxis()->SetTitleSize(0.05);
            graph3->GetYaxis()->SetTitleOffset(0.9);
        } else {
            graph3->Draw("P E SAME");
        }
        graphs3.push_back(graph3);
    }

    auto legend3 = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs3.size(); ++j) {
        if(j==graphs3.size()-1){legend3->AddEntry(graphs3[0], titles3[0], "p");continue;}
        legend3->AddEntry(graphs1[j+1], titles1[j+1], "p");
    }
    legend3->Draw();
    c9->Update();
    c9->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");





    TCanvas *c10 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c10->SetLogx();
    c10->SetGridy();
    const char *titles4[3] = {"U_RooCombineFit", "U_RooFit", "U_RooKeysPdf"};
    int colors4[3] = {kCyan, kOrange,kViolet};
    std::vector<TGraphErrors*> graphs4;

    for (size_t j = 0; j < 3; ++j) {
        TGraphErrors *graph4 = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            if(i==0){cout<<all_bias[i][j+4]<<endl;}
            graph4->SetPoint(i, size[i], all_bias[i][j+4]);
            graph4->SetPointError(i, 0, all_bias_error[i][j+4]);
        }

        graph4->SetMarkerStyle(20 + j);
        graph4->SetMarkerColor(colors4[j]);
        graph4->SetTitle(titles4[j]);

        if (j == 0) {
            graph4->GetYaxis()->SetRangeUser(0.99, 1.02);
            graph4->Draw("APE");
            graph4->SetTitle("Bias vs Template Background N;Template Background N;Bias");
            graph4->GetXaxis()->SetTitleSize(0.05);
            graph4->GetXaxis()->SetTitleOffset(0.85);
            graph4->GetYaxis()->SetTitleSize(0.05);
            graph4->GetYaxis()->SetTitleOffset(0.9);
        } else {
            graph4->Draw("P E SAME");
        }
        graphs4.push_back(graph4);
    }

    auto legend4 = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs4.size(); ++j) {
        //if(j==graphs4.size()-1){legend4->AddEntry(graphs4[0], titles4[0], "p");continue;}
        legend4->AddEntry(graphs4[j], titles4[j], "p");
    }
    legend4->Draw();
    c10->Update();
    c10->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");


    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(0);
    TCanvas *c11 = new TCanvas("c11", "Signal Yields Sigma vs Nmc", 800, 600);
    c11->SetLogx();
    c11->SetGridy();
    const char *titles5[3] = {"U_RooCombineFit", "U_RooFit", "U_RooKeysPdf"};
    int colors5[3] = {kCyan, kOrange,kViolet};
    std::vector<TGraphErrors*> graphs5;

    for (size_t j = 0; j < 3; ++j) {
        TGraphErrors *graph5 = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            graph5->SetPoint(i, size[i], all_fit[i][j+4]->GetParameter(2));
            graph5->SetPointError(i, 0, all_fit[i][j+4]->GetParError(2));
        }

        graph5->SetMarkerStyle(20 + j);
        graph5->SetMarkerColor(colors5[j]);
        graph5->SetTitle(titles5[j]);

        if (j == 0) {
            graph5->GetYaxis()->SetRangeUser(6, 12);
            graph5->Draw("APE");
            graph5->SetTitle("Signal Yields Sigma vs Template Background N;Template Background N;Signal Yields Sigma");
            graph5->GetXaxis()->SetTitleSize(0.05);
            graph5->GetXaxis()->SetTitleOffset(0.85);
            graph5->GetYaxis()->SetTitleSize(0.05);
            graph5->GetYaxis()->SetTitleOffset(0.9);
        } else {
            graph5->Draw("P E SAME");
        }
        graphs5.push_back(graph5);
    }

    auto legend5 = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs5.size(); ++j) {
        //if(j==graphs5.size()-1){legend5->AddEntry(graphs5[0], titles5[0], "p");continue;}
        legend5->AddEntry(graphs5[j], titles5[j], "p");
    }
    legend5->Draw();
    c11->Update();
    c11->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf");
    

    
    
    c1->Print("/data03/tianye/pbar/syscode/pdf/bias_plot8000.pdf]");
        
    

}