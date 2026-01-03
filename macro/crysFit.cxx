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
#include <cstdlib>
#include <vector>
#include <string>
#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TH3F.h"
#include "TMath.h"
#include "RooCrystalBall.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "RooAddPdf.h"
#include <iostream>
#include <map>
#include <stdio.h>
#include "advfit_unbinned.h"
#include "Math/MinimizerOptions.h"
#include "RooMsgService.h"
#include "RooFit.h"

// double combine_crys(double *x, double *p) {
//     TF1 *fFit1=new TF1("fit1",CrystalBall,-1,3,10);
//     TF1 *fFit2=new TF1("fit2",CrystalBall,-1,3,10);
//    //  fFit1->SetParameters(p);
//    //  fFit2->SetParameters(&p[7]);
//     double fsig = fFit1->EvalPar(x, &p[10]);
//     double fbkg = fFit2->EvalPar(x, p);
//     return p[20]*fsig + p[21]*fbkg;
// }


void crysFit(){
    TFile *f=new TFile("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq0001_npi2000.root","read");
    TH1D *hpi = (TH1D*)f->Get("data_pi25_0000");
    TH1D *hpr = (TH1D*)f->Get("data_pr25_0000");
    TH1D *hneg = (TH1D*)f->Get("data_neg25_0000");
    

    TH1D *hpi_clone = (TH1D*)hpi->Clone();
    //hpi_clone->Scale(1/hpi_clone->GetSumOfWeights()/hpi_clone->GetBinWidth(1));

    TH1D *hpr_clone = (TH1D*)hpr->Clone();
    //hpr_clone->Scale(1/hpr_clone->GetSumOfWeights()/hpr_clone->GetBinWidth(1));

    TH1D *hneg_clone = (TH1D*)hneg->Clone();
    //hneg_clone->Scale(1/hneg_clone->GetSumOfWeights()/hneg_clone->GetBinWidth(1));

    RooRealVar x("x", "x", -1, 3);
        int nbins = hpi->GetNbinsX();
    x.setBins(nbins);
    RooDataHist data1("data1", "dataset with x", x, RooFit::Import(*hpi_clone));
    RooDataHist data2("data2", "dataset with x", x, RooFit::Import(*hpr_clone));
    RooDataHist data3("data3", "dataset with x", x, RooFit::Import(*hneg_clone));




    RooRealVar xs_0("xs_0", "xs_0",hpr_clone->GetMean(), 0.5, 1.3);
    RooRealVar xs_1("xs_1", "sigmaL",hpr_clone->GetRMS(), 0.05, 1.1);
    RooRealVar xs_2("xs_2", "sigmaR",hpr_clone->GetRMS(), 0.05, 1.1);
    RooRealVar xs_3("xs_3", "alphaL",1, 0.01, 5.0);
    RooRealVar xs_4("xs_4", "nL",2,1.0001, 100);
    RooRealVar xs_5("xs_5", "alphaR",1, 0.01, 5.0);
    RooRealVar xs_6("xs_6", "nR",2,1.0001, 100);



    RooRealVar xb_0("xb_0", "xb_0",0.001, 0.001, 0.35);
    RooRealVar xb_1("xb_1", "sigmaL",hpi_clone->GetRMS(), 0.05, 1.1);
    RooRealVar xb_2("xb_2", "sigmaR",hpi_clone->GetRMS(), 0.05, 1.1);
    RooRealVar xb_3("xb_3", "alphaL",1, 0.01, 5.0);
    RooRealVar xb_4("xb_4", "nL",2,1.0001, 100);
    RooRealVar xb_5("xb_5", "alphaR",1, 0.01, 5.0);
    RooRealVar xb_6("xb_6", "nR",2,1.0001, 100);

    RooRealVar scal_sig("scal_sig", "scal_sig", 100, 10, 1000000);
    RooRealVar scal_bkg("scal_bkg", "scal_bkg", 100, 10, 1000000);

    RooCrystalBall cb1("cb1", "cb1", x, xs_0, xs_1,xs_2, xs_3, xs_4, xs_5, xs_6);
    RooCrystalBall cb2("cb2", "cb2", x, xb_0, xb_1,xb_2, xb_3, xb_4, xb_5, xb_6);

    RooAddPdf model("model", "model", RooArgList(cb1, cb2), RooArgList(scal_sig, scal_bkg));
    
    cb1.fitTo(data2);
    cb2.fitTo(data1);


    xs_0.setConstant(kTRUE);
    xs_1.setConstant(kTRUE);
    xs_2.setConstant(kTRUE);
    xs_3.setConstant(kTRUE);
    xs_4.setConstant(kTRUE);
    xs_5.setConstant(kTRUE);
    xs_6.setConstant(kTRUE);

    xb_0.setConstant(kTRUE);
    xb_1.setConstant(kTRUE);
    xb_2.setConstant(kTRUE);
    xb_3.setConstant(kTRUE);
    xb_4.setConstant(kTRUE);
    xb_5.setConstant(kTRUE);
    xb_6.setConstant(kTRUE);

    model.fitTo(data3);

    RooPlot* frame1 = x.frame();
    data1.plotOn(frame1);
    cb2.plotOn(frame1);

    RooPlot* frame2 = x.frame();
    data2.plotOn(frame2);
    cb1.plotOn(frame2);

    RooPlot* frame3 = x.frame();
    data3.plotOn(frame3);
    model.plotOn(frame3);

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(2, 2);

    c->cd(1);
    frame1->Draw();

    c->cd(2);
    frame2->Draw();

    c->cd(3);
    frame3->Draw();

    cout<<"nsig: "<<scal_sig.getVal()<<endl;
    cout<<"nbkg: "<<scal_bkg.getVal()<<endl;

    // c->cd(1);
    // frame1->Draw();

    // c->cd(2);
    // frame2->Draw();

    c->Print("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/pdf/crys_results.pdf");
}