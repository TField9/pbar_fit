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

void bias_plot(){
    gStyle->SetLabelSize(0.1, "X");
        TFile *f=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/binned_roocomb_seq0001_npi0100.root","read");
        TFile *f1=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/binned_TFF_seq0001_npi0100.root","read");
        //TFile *f2=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/binned_SPD_seq0001_npi0100.root","read");
        TFile *f2=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/binned_SPD_seq0001_npi0100.root","read");
        TFile *f3=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/binned_IWLS1_seq0001_npi0100.root","read");
        TFile *f4=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/unbinned_roocomb_seq0001_npi0100.root","read");
        TFile *f5=new TFile("/data03/tianye/pbar/root/bin25_Nmc100/unbinned_roofit_seq0001_npi0100.root","read");


        std::vector<double>* params1 = new std::vector<double>(22);
        std::vector<double>* params2 = new std::vector<double>(22);
        std::vector<double>* params3 = new std::vector<double>(4);
        std::vector<double>* params4 = new std::vector<double>(4);
        std::vector<double>* params5 = new std::vector<double>(4);
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
                bin20_SPD->Fill((*params4)[0] * 400);
                bin100_SPD->Fill((*params5)[0] * 400);
                bin20_IWLS->Fill((*params6)[0] * 400);
                bin100_IWLS->Fill((*params6)[1] * 400);
                unbinned_comb->Fill((*params7)[10] * 400);
                unbinned_roofit->Fill((*params8)[0]);



                bin20_roocomb_error->Fill((*params1)[21] * 400);
                bin100_roocomb_error->Fill((*params2)[21] * 400);
                bin20_TFF_error->Fill((*params3)[2] * 400);
                bin100_TFF_error->Fill((*params3)[3] * 400);
                bin20_SPD_error->Fill((*params4)[2] * 400);
                bin100_SPD_error->Fill((*params5)[2] * 400);
                bin20_IWLS_error->Fill((*params6)[2] * 400);
                bin100_IWLS_error->Fill((*params6)[3] * 400);
                unbinned_comb_error->Fill((*params7)[21] * 400);
                unbinned_roofit_error->Fill((*params8)[1]);

                bin20_roocomb_pull->Fill(((*params1)[10] - 0.5)/(*params1)[21] );
                bin100_roocomb_pull->Fill(((*params2)[10] - 0.5)/(*params2)[21] );
                bin20_TFF_pull->Fill(((*params3)[0] - 0.5)/(*params3)[2] );
                bin100_TFF_pull->Fill(((*params3)[1] - 0.5)/(*params3)[3] );
                bin20_SPD_pull->Fill(((*params4)[0] - 0.5)/(*params4)[2] );
                bin100_SPD_pull->Fill(((*params5)[0] - 0.5)/(*params5)[2] );
                bin20_IWLS_pull->Fill(((*params6)[0] - 0.5)/(*params6)[2] );
                bin100_IWLS_pull->Fill(((*params6)[1] - 0.5)/(*params6)[3] );
                unbinned_comb_pull->Fill(((*params7)[10] - 0.5)/(*params7)[21] );
                unbinned_roofit_pull->Fill(((*params8)[0] - 200)/(*params8)[1] );
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

            // Draw histograms
            gStyle->SetOptStat(0);
            TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
            c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/bias_plot.pdf[");
            c1->Divide(3, 4);
            c1->cd(1);
            bin20_roocomb->SetLineColor(kRed);
            bin20_roocomb->GetXaxis()->SetNdivisions(504, kTRUE); // Set main divisions and labels
            bin20_roocomb->GetYaxis()->SetNdivisions(505, kTRUE); // Set main divisions and labels
            bin20_roocomb->GetXaxis()->SetLabelSize(0.06); // Set X-axis label size
            bin20_roocomb->GetYaxis()->SetLabelSize(0.06); // Set Y-axis label size
            bin20_roocomb->GetXaxis()->SetTitle("Signal Yield");
            bin20_roocomb->GetYaxis()->SetTitle("Count");
            bin20_roocomb->GetYaxis()->SetTitleOffset(0.85); // Set Y-axis title offset
            bin20_roocomb->GetXaxis()->SetTitleSize(0.05); // Increase X-axis title size
            bin20_roocomb->GetYaxis()->SetTitleSize(0.05); // Increase Y-axis title size
            bin20_roocomb->Draw();
            gaussFit1->SetLineColor(kBlue);
            gaussFit1->Draw("same");
            TLegend *leg1 = new TLegend(0.55, 0.5, 0.88, 0.8);
            leg1->AddEntry(bin20_roocomb, "bin20_roocomb", "l");
            leg1->AddEntry(
                gaussFit1,
                Form("#splitline{mean=%.2f #pm %.2f}{sigma=%.2f}",
                        gaussFit1->GetParameter(1),
                        gaussFit1->GetParError(1),
                        gaussFit1->GetParameter(2)),
                    "l"
                );  
            leg1->SetTextSize(0.05);
            leg1->SetBorderSize(0); // Remove the border
            leg1->Draw();

            c1->cd(2);
            bin20_TFF->SetLineColor(kGreen);
            bin20_TFF->GetXaxis()->SetNdivisions(5);
            bin20_TFF->GetYaxis()->SetNdivisions(5);
            bin20_TFF->Draw();
            gaussFit3->SetLineColor(kMagenta);
            gaussFit3->Draw("same");
            TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg3->AddEntry(bin20_TFF, "bin20_TFF", "l");
            leg3->AddEntry(gaussFit3, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit3->GetParameter(1), gaussFit3->GetParError(1), gaussFit3->GetParameter(2)), "l");
            leg3->SetTextSize(0.04);
            leg3->Draw();

            c1->cd(3);
            bin20_SPD->SetLineColor(kCyan);
            bin20_SPD->GetXaxis()->SetNdivisions(5);
            bin20_SPD->GetYaxis()->SetNdivisions(5);
            bin20_SPD->Draw();
            gaussFit5->SetLineColor(kOrange);
            gaussFit5->Draw("same");
            TLegend *leg5 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg5->AddEntry(bin20_SPD, "bin20_paper", "l");
            leg5->AddEntry(gaussFit5, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit5->GetParameter(1), gaussFit5->GetParError(1), gaussFit5->GetParameter(2)), "l");
            leg5->SetTextSize(0.04);
            leg5->Draw();

            c1->cd(4);
            bin20_IWLS->SetLineColor(kYellow + 1);
            bin20_IWLS->GetXaxis()->SetNdivisions(5);
            bin20_IWLS->GetYaxis()->SetNdivisions(5);
            bin20_IWLS->Draw();
            gaussFit7->SetLineColor(kViolet);
            gaussFit7->Draw("same");
            TLegend *leg7 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg7->AddEntry(bin20_IWLS, "bin20_IWLS", "l");
            leg7->AddEntry(gaussFit7, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit7->GetParameter(1), gaussFit7->GetParError(1), gaussFit7->GetParameter(2)), "l");
            leg7->SetTextSize(0.04);
            leg7->Draw();

            c1->cd(5);
            unbinned_comb->SetLineColor(kPink);
            unbinned_comb->GetXaxis()->SetNdivisions(5);
            unbinned_comb->GetYaxis()->SetNdivisions(5);
            unbinned_comb->Draw();
            gaussFit9->SetLineColor(kTeal);
            gaussFit9->Draw("same");
            TLegend *leg9 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg9->AddEntry(unbinned_comb, "unbinned_comb", "l");
            leg9->AddEntry(gaussFit9, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit9->GetParameter(1), gaussFit9->GetParError(1), gaussFit9->GetParameter(2)), "l");
            leg9->SetTextSize(0.04);
            leg9->Draw();

            c1->cd(6);
            unbinned_roofit->SetLineColor(kAzure);
            unbinned_roofit->GetXaxis()->SetNdivisions(5);
            unbinned_roofit->GetYaxis()->SetNdivisions(5);
            unbinned_roofit->Draw();
            gaussFit10->SetLineColor(kSpring);
            gaussFit10->Draw("same");
            TLegend *leg10 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg10->AddEntry(unbinned_roofit, "unbinned_roofit", "l");
            leg10->AddEntry(gaussFit10, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit10->GetParameter(1), gaussFit10->GetParError(1), gaussFit10->GetParameter(2)), "l");
            leg10->SetTextSize(0.04);
            leg10->Draw();

            c1->cd(7);
            bin100_roocomb->SetLineColor(kRed + 2);
            bin100_roocomb->GetXaxis()->SetNdivisions(5);
            bin100_roocomb->GetYaxis()->SetNdivisions(5);
            bin100_roocomb->Draw();
            gaussFit2->SetLineColor(kBlue + 2);
            gaussFit2->Draw("same");
            TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg2->AddEntry(bin100_roocomb, "bin100_roocomb", "l");
            leg2->AddEntry(gaussFit2, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit2->GetParameter(1), gaussFit2->GetParError(1), gaussFit2->GetParameter(2)), "l");
            leg2->SetTextSize(0.04);
            leg2->Draw();

            c1->cd(8);
            bin100_TFF->SetLineColor(kGreen + 2);
            bin100_TFF->GetXaxis()->SetNdivisions(5);
            bin100_TFF->GetYaxis()->SetNdivisions(5);
            bin100_TFF->Draw();
            gaussFit4->SetLineColor(kMagenta + 2);
            gaussFit4->Draw("same");
            TLegend *leg4 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg4->AddEntry(bin100_TFF, "bin100_TFF", "l");
            leg4->AddEntry(gaussFit4, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit4->GetParameter(1), gaussFit4->GetParError(1), gaussFit4->GetParameter(2)), "l");
            leg4->SetTextSize(0.04);
            leg4->Draw();

            c1->cd(9);
            bin100_SPD->SetLineColor(kCyan + 2);
            bin100_SPD->GetXaxis()->SetNdivisions(5);
            bin100_SPD->GetYaxis()->SetNdivisions(5);
            bin100_SPD->Draw();
            gaussFit6->SetLineColor(kOrange + 2);
            gaussFit6->Draw("same");
            TLegend *leg6 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg6->AddEntry(bin100_SPD, "bin100_paper", "l");
            leg6->AddEntry(gaussFit6, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit6->GetParameter(1), gaussFit6->GetParError(1), gaussFit6->GetParameter(2)), "l");
            leg6->SetTextSize(0.04);
            leg6->Draw();

            c1->cd(10);
            bin100_IWLS->SetLineColor(kYellow + 2);
            bin100_IWLS->GetXaxis()->SetNdivisions(5);
            bin100_IWLS->GetYaxis()->SetNdivisions(5);
            bin100_IWLS->Draw();
            gaussFit8->SetLineColor(kViolet + 2);
            gaussFit8->Draw("same");
            TLegend *leg8 = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg8->AddEntry(bin100_IWLS, "bin100_IWLS", "l");
            leg8->AddEntry(gaussFit8, Form("mean=%.2f #pm %.2f, sigma=%.2f", gaussFit8->GetParameter(1), gaussFit8->GetParError(1), gaussFit8->GetParameter(2)), "l");
            leg8->SetTextSize(0.04);
            leg8->Draw();
            
            for (int i = 1; i <= 12; i++) {
                TPad *pad = (TPad*)c1->GetPad(i);

                // 获取原始位置
                double xlow, ylow, xup, yup;
                pad->GetPadPar(xlow, ylow, xup, yup);

                // 计算缩小后的中心和新尺寸
                double xcenter = (xlow + xup) / 2.0;
                double ycenter = (ylow + yup) / 2.0;
                double xwidth = (xup - xlow) * 0.9; // 宽度缩小到90%
                double yheight = (yup - ylow) * 0.9; // 高度缩小到90%

                // 设置新的Pad位置
                pad->SetPad(xcenter - xwidth / 2.0, ycenter - yheight / 2.0-0.02,
                            xcenter + xwidth / 2.0, ycenter + yheight / 2.0-0.02);

            }

            c1->cd();
            //c1->SetTopMargin(0.15); // Increase the top margin to leave more space for the title
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.02);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.98, "Signal Yields Distribution for Different Binning Methods(N_bkg template =100)");
            c1->Update();
            c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/bias_plot.pdf");

            TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
            c2->Divide(3, 4);

            c2->cd(1);
            bin20_roocomb_error->Draw();

            c2->cd(2);
            bin20_TFF_error->Draw();

            c2->cd(3);
            bin20_SPD_error->Draw();

            c2->cd(4);
            bin20_IWLS_error->Draw();

            c2->cd(5);
            unbinned_comb_error->Draw();

            c2->cd(6);
            unbinned_roofit_error->Draw();

            c2->cd(7);
            bin100_roocomb_error->Draw();

            c2->cd(8);
            bin100_TFF_error->Draw();

            c2->cd(9);
            bin100_SPD_error->Draw();

            c2->cd(10);
            bin100_IWLS_error->Draw();

            c2->cd();
            TLatex latex1;
            latex1.SetNDC();
            latex1.SetTextSize(0.02);
            latex1.SetTextAlign(22); // Center alignment
            latex1.DrawLatex(0.5, 0.98, "Signal Yields error for Different Binning Methods");

            c2->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/bias_plot.pdf");
            
            TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
            c3->Divide(3, 4);
            c3->cd(1);
            bin20_roocomb_pull->Draw();
            TLegend *leg1_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg1_pull->AddEntry(bin20_roocomb_pull, "bin20_roocomb_pull", "l");
            leg1_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_roocomb_pull->GetMean()), "");
            leg1_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_roocomb_pull->GetStdDev()), "");
            leg1_pull->SetTextSize(0.05);
            leg1_pull->Draw();

            c3->cd(2);
            bin20_TFF_pull->Draw();
            TLegend *leg3_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg3_pull->AddEntry(bin20_TFF_pull, "bin20_TFF_pull", "l");
            leg3_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_TFF_pull->GetMean()), "");
            leg3_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_TFF_pull->GetStdDev()), "");
            leg3_pull->SetTextSize(0.05);
            leg3_pull->Draw();

            c3->cd(3);
            bin20_SPD_pull->Draw();
            TLegend *leg5_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg5_pull->AddEntry(bin20_SPD_pull, "bin20_paper_pull", "l");
            leg5_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_SPD_pull->GetMean()), "");
            leg5_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_SPD_pull->GetStdDev()), "");
            leg5_pull->SetTextSize(0.05);
            leg5_pull->Draw();

            c3->cd(4);
            bin20_IWLS_pull->Draw();
            TLegend *leg7_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg7_pull->AddEntry(bin20_IWLS_pull, "bin20_IWLS_pull", "l");
            leg7_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_IWLS_pull->GetMean()), "");
            leg7_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_IWLS_pull->GetStdDev()), "");
            leg7_pull->SetTextSize(0.05);
            leg7_pull->Draw();

            c3->cd(5);
            unbinned_comb_pull->Draw();
            TLegend *leg9_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg9_pull->AddEntry(unbinned_comb_pull, "unbinned_comb_pull", "l");
            leg9_pull->AddEntry((TObject*)0, Form("Mean = %.2f", unbinned_comb_pull->GetMean()), "");
            leg9_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", unbinned_comb_pull->GetStdDev()), "");
            leg9_pull->SetTextSize(0.05);
            leg9_pull->Draw();

            c3->cd(6);
            unbinned_roofit_pull->Draw();
            TLegend *leg10_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg10_pull->AddEntry(unbinned_roofit_pull, "unbinned_roofit_pull", "l");
            leg10_pull->AddEntry((TObject*)0, Form("Mean = %.2f", unbinned_roofit_pull->GetMean()), "");
            leg10_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", unbinned_roofit_pull->GetStdDev()), "");
            leg10_pull->SetTextSize(0.05);
            leg10_pull->Draw();

            c3->cd(7);
            bin100_roocomb_pull->Draw();
            TLegend *leg2_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg2_pull->AddEntry(bin100_roocomb_pull, "bin100_roocomb_pull", "l");
            leg2_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_roocomb_pull->GetMean()), "");
            leg2_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_roocomb_pull->GetStdDev()), "");
            leg2_pull->SetTextSize(0.05);
            leg2_pull->Draw();

            c3->cd(8);
            bin100_TFF_pull->Draw();
            TLegend *leg4_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg4_pull->AddEntry(bin100_TFF_pull, "bin100_TFF_pull", "l");
            leg4_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_TFF_pull->GetMean()), "");
            leg4_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_TFF_pull->GetStdDev()), "");
            leg4_pull->SetTextSize(0.05);
            leg4_pull->Draw();

            c3->cd(9);
            bin100_SPD_pull->Draw();
            TLegend *leg6_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg6_pull->AddEntry(bin100_SPD_pull, "bin100_paper_pull", "l");
            leg6_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_SPD_pull->GetMean()), "");
            leg6_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_SPD_pull->GetStdDev()), "");
            leg6_pull->SetTextSize(0.05);
            leg6_pull->Draw();

            c3->cd(10);
            bin100_IWLS_pull->Draw();
            TLegend *leg8_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            leg8_pull->AddEntry(bin100_IWLS_pull, "bin100_IWLS_pull", "l");
            leg8_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_IWLS_pull->GetMean()), "");
            leg8_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_IWLS_pull->GetStdDev()), "");
            leg8_pull->SetTextSize(0.05);
            leg8_pull->Draw();

            // Display mean and standard deviation for pull distributions
            TLatex latex3;
            latex3.SetNDC();
            latex3.SetTextSize(0.02);
            latex3.SetTextAlign(22); // Center alignment

            // c3->cd(1);
            // bin20_roocomb_pull->Draw();
            // TLegend *leg1_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg1_pull->AddEntry(bin20_roocomb_pull, "bin20_roocomb_pull", "l");
            // leg1_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_roocomb_pull->GetMean()), "");
            // leg1_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_roocomb_pull->GetStdDev()), "");
            // leg1_pull->SetTextSize(0.05);
            // leg1_pull->Draw();

            // c3->cd(2);
            // bin20_TFF_pull->Draw();
            // TLegend *leg3_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg3_pull->AddEntry(bin20_TFF_pull, "bin20_TFF_pull", "l");
            // leg3_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_TFF_pull->GetMean()), "");
            // leg3_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_TFF_pull->GetStdDev()), "");
            // leg3_pull->SetTextSize(0.05);
            // leg3_pull->Draw();

            // c3->cd(3);
            // bin20_SPD_pull->Draw();
            // TLegend *leg5_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg5_pull->AddEntry(bin20_SPD_pull, "bin20_paper_pull", "l");
            // leg5_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_SPD_pull->GetMean()), "");
            // leg5_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_SPD_pull->GetStdDev()), "");
            // leg5_pull->SetTextSize(0.05);
            // leg5_pull->Draw();

            // c3->cd(4);
            // bin20_IWLS_pull->Draw();
            // TLegend *leg7_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg7_pull->AddEntry(bin20_IWLS_pull, "bin20_IWLS_pull", "l");
            // leg7_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin20_IWLS_pull->GetMean()), "");
            // leg7_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin20_IWLS_pull->GetStdDev()), "");
            // leg7_pull->SetTextSize(0.05);
            // leg7_pull->Draw();

            // c3->cd(5);
            // unbinned_comb_pull->Draw();
            // TLegend *leg9_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg9_pull->AddEntry(unbinned_comb_pull, "unbinned_comb_pull", "l");
            // leg9_pull->AddEntry((TObject*)0, Form("Mean = %.2f", unbinned_comb_pull->GetMean()), "");
            // leg9_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", unbinned_comb_pull->GetStdDev()), "");
            // leg9_pull->SetTextSize(0.05);
            // leg9_pull->Draw();

            // c3->cd(6);
            // unbinned_roofit_pull->Draw();
            // TLegend *leg10_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg10_pull->AddEntry(unbinned_roofit_pull, "unbinned_roofit_pull", "l");
            // leg10_pull->AddEntry((TObject*)0, Form("Mean = %.2f", unbinned_roofit_pull->GetMean()), "");
            // leg10_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", unbinned_roofit_pull->GetStdDev()), "");
            // leg10_pull->SetTextSize(0.05);
            // leg10_pull->Draw();

            // c3->cd(7);
            // bin100_roocomb_pull->Draw();
            // TLegend *leg2_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg2_pull->AddEntry(bin100_roocomb_pull, "bin100_roocomb_pull", "l");
            // leg2_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_roocomb_pull->GetMean()), "");
            // leg2_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_roocomb_pull->GetStdDev()), "");
            // leg2_pull->SetTextSize(0.05);
            // leg2_pull->Draw();

            // c3->cd(8);
            // bin100_TFF_pull->Draw();
            // TLegend *leg4_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg4_pull->AddEntry(bin100_TFF_pull, "bin100_TFF_pull", "l");
            // leg4_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_TFF_pull->GetMean()), "");
            // leg4_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_TFF_pull->GetStdDev()), "");
            // leg4_pull->SetTextSize(0.05);
            // leg4_pull->Draw();

            // c3->cd(9);
            // bin100_SPD_pull->Draw();
            // TLegend *leg6_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg6_pull->AddEntry(bin100_SPD_pull, "bin100_paper_pull", "l");
            // leg6_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_SPD_pull->GetMean()), "");
            // leg6_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_SPD_pull->GetStdDev()), "");
            // leg6_pull->SetTextSize(0.05);
            // leg6_pull->Draw();

            // c3->cd(10);
            // bin100_IWLS_pull->Draw();
            // TLegend *leg8_pull = new TLegend(0.6, 0.7, 0.9, 0.9);
            // leg8_pull->AddEntry(bin100_IWLS_pull, "bin100_IWLS_pull", "l");
            // leg8_pull->AddEntry((TObject*)0, Form("Mean = %.2f", bin100_IWLS_pull->GetMean()), "");
            // leg8_pull->AddEntry((TObject*)0, Form("StdDev = %.2f", bin100_IWLS_pull->GetStdDev()), "");
            // leg8_pull->SetTextSize(0.05);
            // leg8_pull->Draw();

            c3->cd();
            TLatex latex2;
            latex2.SetNDC();
            latex2.SetTextSize(0.02);
            latex2.SetTextAlign(22); // Center alignment
            latex2.DrawLatex(0.5, 0.98, "Pull Distribution for Different Binning Methods");

            c3->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/bias_plot.pdf");
            
            c1->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/pdf/bias_plot.pdf]");

}