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
#include "TF1.h"
#include "TLegend.h"

void kde_diffwidth(){
    std::vector<int> size = {10,15,20,25,30,40,50,75,100,150,200,300,400,500,700,1000,1500,2000,3000,4000,5000,8000};
    std::vector<std::vector<TH1D*>> all_sig_yield;
    std::vector<std::vector<double>> all_bias;
    std::vector<std::vector<double>> all_bias_error;
    std::vector<std::vector<TF1*>> all_fit;

    // Kernel widths to test
    std::vector<double> rho = {0.7, 0.9, 1.0, 1.1, 1.3};

    // Loop over different Nmc values
    for (int N : size){ 
        std::vector<TH1D*> sig_yield;
        std::vector<double> bias;
        std::vector<double> bias_error;
        std::vector<TF1*> fits;

        // Open files for combined RooFit and RooKeysPdf
        TFile *f_comb = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roocomb_seq0001_npi%04d.root", N), "read");
        TFile *f_keys = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_rookeyspdf_diffwidth_seq0001_npi%04d.root", N), "read");

        // Get trees
        TTree *t_comb = (TTree*)f_comb->Get("fitParamsTree");
        std::vector<TTree*> t_keys;
        for (double w : rho) {
            t_keys.push_back((TTree*)f_keys->Get(Form("fitParamsTree_rho%.1f", w)));
        }

        // Set branch addresses
        std::vector<double>* params_comb = nullptr;
        t_comb->SetBranchAddress("fitparams_unbinned", &params_comb);
        std::vector<std::vector<double>*> params_keys(rho.size(), nullptr);
        for (size_t k = 0; k < rho.size(); ++k) {
            t_keys[k]->SetBranchAddress("fitparams_unbinned", &params_keys[k]);
        }

        // Create histograms
        TH1D* h_comb = new TH1D("h_comb", "Unbinned RooComb", 200, 140, 280);
        std::vector<TH1D*> h_keys;
        for (double w : rho) {
            h_keys.push_back(new TH1D(Form("h_key%.1f", w), Form("Unbinned RooKeysPdf %.1f", w), 200, 140, 280));
        }

        // Fill histograms
        int nEntries = t_comb->GetEntries();
        for (int i = 0; i < nEntries; ++i) {
            t_comb->GetEntry(i);
            h_comb->Fill((*params_comb)[10] * 400);
            for (size_t k = 0; k < rho.size(); ++k) {
                t_keys[k]->GetEntry(i);
                h_keys[k]->Fill((*params_keys[k])[0]);
            }
        }

        // Fit RooComb histogram
        TF1* fit_comb = new TF1("fit_comb", "gaus", 100, 300);
        h_comb->Fit(fit_comb, "Q");
        sig_yield.push_back(h_comb);
        bias.push_back(fit_comb->GetParameter(1) / 200);
        bias_error.push_back(fit_comb->GetParError(1) / 200);
        fits.push_back(fit_comb);

        // Fit RooKeysPdf histograms
        for (size_t k = 0; k < rho.size(); ++k) {
            TF1* f = new TF1(Form("fit_key%.1f", rho[k]), "gaus", 100, 300);
            h_keys[k]->Fit(f, "Q");
            sig_yield.push_back(h_keys[k]);
            bias.push_back(f->GetParameter(1) / 200);
            bias_error.push_back(f->GetParError(1) / 200);
            fits.push_back(f);
        }

        all_sig_yield.push_back(sig_yield);
        all_bias.push_back(bias);
        all_bias_error.push_back(bias_error);
        all_fit.push_back(fits);

        // Close files
        f_comb->Close();
        f_keys->Close();
    }

    // Plot Bias vs Nmc
    TCanvas *c2 = new TCanvas("c2", "Bias vs Nmc", 800, 600);
    c2->SetLogx(); c2->SetGridy();
    const char *titles_bias[] = {"U_RooComb", "U_RooKeysPdf 0.7", "U_RooKeysPdf 0.9", "U_RooKeysPdf 1.0", "U_RooKeysPdf 1.1", "U_RooKeysPdf 1.3"};
    int markers_bias[] = {24,21,22,23,20,25};
    std::vector<TGraphErrors*> graphs_bias;

    for (int j = 0; j < (int)titles_bias->length(); ++j) {
        TGraphErrors *g = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            g->SetPoint(i, size[i], all_bias[i][j]);
            g->SetPointError(i, 0, all_bias_error[i][j]);
        }
        g->SetMarkerStyle(markers_bias[j]);
        if (j == 0) {
            g->Draw("APE");
            g->GetYaxis()->SetRangeUser(0.96, 1.03);
        } else {
            g->Draw("PE SAME");
        }
        graphs_bias.push_back(g);
    }
    auto legend_bias = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs_bias.size(); ++j) {
        legend_bias->AddEntry(graphs_bias[j], titles_bias[j], "p");
    }
    legend_bias->Draw();
    c2->Update();
    c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/kde_diffwidth_bias.pdf");

    // Plot Sigma vs Nmc
    TCanvas *c3 = new TCanvas("c3", "Sigma vs Nmc", 800, 600);
    c3->SetLogx(); c3->SetGridy();
    const char *titles_sigma[] = {"U_RooComb", "U_RooKeysPdf 0.7", "U_RooKeysPdf 0.9", "U_RooKeysPdf 1.0", "U_RooKeysPdf 1.1", "U_RooKeysPdf 1.3"};
    int markers_sigma[] = {24,21,22,23,20,25};
    std::vector<TGraphErrors*> graphs_sigma;

    for (int j = 0; j < (int)titles_sigma->length(); ++j) {
        TGraphErrors *g = new TGraphErrors(size.size());
        for (size_t i = 0; i < size.size(); ++i) {
            double sigma = all_fit[i][j]->GetParameter(2);
            double sigma_err = all_fit[i][j]->GetParError(2);
            g->SetPoint(i, size[i], sigma);
            g->SetPointError(i, 0, sigma_err);
        }
        g->SetMarkerStyle(markers_sigma[j]);
        if (j == 0) {
            g->Draw("APE");
            g->GetYaxis()->SetRangeUser(6, 11.5);
        } else {
            g->Draw("PE SAME");
        }
        graphs_sigma.push_back(g);
    }
    auto legend_sigma = new TLegend(0.5, 0.6, 0.8, 0.9);
    for (size_t j = 0; j < graphs_sigma.size(); ++j) {
        legend_sigma->AddEntry(graphs_sigma[j], titles_sigma[j], "p");
    }
    legend_sigma->Draw();
    c3->Update();
    c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/kde_diffwidth_sigma.pdf");
}
