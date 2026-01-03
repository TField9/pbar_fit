#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include <vector>
#include <iostream>

void RooKPdf(RooDataSet &sigData, RooDataSet &bkgData,RooDataSet &negData,
             double kde_wide, double &nsig, double &nsig_error) {
    RooRealVar x("x","x",-1.0,3.0);
    // 构建 KDE
    RooKeysPdf f_sig("f_sig","", x, sigData, RooKeysPdf::MirrorBoth, kde_wide);
    RooKeysPdf f_bkg("f_bkg","", x, bkgData, RooKeysPdf::MirrorBoth, kde_wide);

    RooRealVar scal_sig("scal_sig","",0.5,0.0,1.0);
    RooAddPdf combinedPdf("combinedPdf","", RooArgList(f_sig,f_bkg), RooArgList(scal_sig));


    // 拟合并计算 nsig
    RooFitResult* res = combinedPdf.fitTo(
        negData,
        RooFit::Save(true),
        RooFit::PrintLevel(-1),
        RooFit::Verbose(false),
        RooFit::Minimizer("Minuit2","migrad")
    );
    Int_t Ntot = negData.sumEntries();
    if (res && res->status()==0) {
        nsig       = scal_sig.getVal()  * Ntot;
        nsig_error = scal_sig.getError()* Ntot;
    } else {
        nsig = nsig_error = -1;
    }
}

void RKP(int npr,int npi){
    TStopwatch timer; timer.Start();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().silentMode();

    TFile *f = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npr),"read");
    TFile *f1 = new TFile(Form("/data03/tianye/pbar/root/both_vary/toymc_bin1_seq0001_npi%04d.root",npi),"read");
    TFile *fOut = TFile::Open(
        Form("/data03/tianye/pbar/root/2d_result/"
             "unbinned_rookeyspdf_diffwidth_seq0001_npr%04d_npi%04d.root",
           npr,npi),
      "RECREATE");


    std::vector<double> rho = {0.7, 0.9, 1.0, 1.1, 1.3};
    std::vector<std::vector<double>> nsig(rho.size(), std::vector<double>(2));
    TTree *treeOut = new TTree("fitParamsTree","Fit Results");
    treeOut->Branch("nsig",&nsig);

    for (int j = 0; j < 10000; ++j) {
        // --- 填充两套模板的 sigData/bkgData ---
        RooRealVar x("x","x",-1.0,3.0);
        RooDataSet sigData0("sig0","",RooArgSet(x)), bkgData0("bkg0","",RooArgSet(x));
        RooDataSet negData0("neg0","",RooArgSet(x));
        Double_t v;

        // 第一套模板
        TTree *s0 = (TTree*)f->Get(Form("sigtree_%04d",j));
        TTree *b0 = (TTree*)f1->Get(Form("bkgtree_%04d",j));
        s0->SetBranchAddress("pr_mass2",&v);
        for (Long64_t i=0,n=s0->GetEntries(); i<n; ++i) {
            s0->GetEntry(i); x.setVal(v); sigData0.add(RooArgSet(x));
        }
        b0->SetBranchAddress("pi_mass2",&v);
        for (Long64_t i=0,n=b0->GetEntries(); i<n; ++i) {
            b0->GetEntry(i); x.setVal(v); bkgData0.add(RooArgSet(x));
        }


        // --- 对每个 rho 和每种模板、每个 negtree 做拟合 ---
        TTree *n0 = (TTree*)f1->Get(Form("negtree_%04d",j));
        n0->SetBranchAddress("neg_mass2",&v);
        for (Long64_t i=0,n=n0->GetEntries(); i<n; ++i) {
            n0->GetEntry(i); x.setVal(v); negData0.add(RooArgSet(x));
        }

        for (size_t i = 0; i < rho.size(); ++i) {
            RooKPdf(sigData0,bkgData0,negData0,rho[i],
                    nsig[i][0], nsig[i][1]);
        }

        treeOut->Fill();
        if(j%100==0){
            std::cout<<"j = "<<j<<std::endl;
            std::cout<<"sig_yields = "<<nsig[2][0]<<" "<<nsig[2][1]<<std::endl;
        }
        delete s0;
        delete b0;
        delete n0;

    }

    fOut->cd();
    treeOut->Write();
    fOut->Close();
    std::cout<<"CPU="<<timer.CpuTime()
             <<"s, Real="<<timer.RealTime()<<"s\n";
}

int main(int argc, char* argv[]){
    RKP(atoi(argv[1]),atoi(argv[2]));
    return 0;
}
