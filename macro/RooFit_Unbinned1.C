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
#include <fstream>
#include <chrono>
#include <thread>
#include <RooFitResult.h>

 class expgausexp : public RooAbsPdf {
    public:
        expgausexp() {}
        expgausexp(const char *name, const char *title,
                   RooAbsReal& _x,
                   RooAbsReal& _x0,
                   RooAbsReal& _sigmaL,
                   RooAbsReal& _alphaL,
                   RooAbsReal& _sigmaR,
                   RooAbsReal& _alphaR) :
            RooAbsPdf(name, title),
            x("x", "x", this, _x),
            x0("x0", "x0", this, _x0),
            sigmaL("sigmaL", "sigmaL", this, _sigmaL),
            alphaL("alphaL", "alphaL", this, _alphaL),
            sigmaR("sigmaR", "sigmaR", this, _sigmaR),
            alphaR("alphaR", "alphaR", this, _alphaR) {}

        expgausexp(const expgausexp& other, const char* name = 0) :
            RooAbsPdf(other, name),
            x("x", this, other.x),
            x0("x0", this, other.x0),
            sigmaL("sigmaL", this, other.sigmaL),
            alphaL("alphaL", this, other.alphaL),
            sigmaR("sigmaR", this, other.sigmaR),
            alphaR("alphaR", this, other.alphaR) {}

        virtual TObject* clone(const char* newname) const { return new expgausexp(*this, newname); }
        inline virtual ~expgausexp() {}

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
            if (matchArgs(allVars, analVars, x)) return 1;
            return 0;
        }

        Double_t analyticalIntegral(Int_t code, const char* rangeName) const {
            assert(code == 1);
            constexpr double sqrtPiOver2 = 1.2533141373;
            constexpr double sqrt2 = 1.4142135624;
            double xmin = x.min(rangeName);
            double xmax = x.max(rangeName);
            double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
            double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
            double sum = 0;
            if (tmin < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double lv = tmin;
                double uv = std::min(tmax, -alphaL->getVal());
                sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
            }
            if (tmax > alphaR) {
                double a = 0.5 * alphaR * alphaR;
                double lv = std::max(tmin, alphaR->getVal());
                double uv = tmax;
                sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
            }
            if (tmin < alphaR && tmax > -alphaL) {
                double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
                double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
                sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR->getVal()) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
            }
            return sum;
        }

    protected:
        RooRealProxy x;
        RooRealProxy x0;
        RooRealProxy sigmaL;
        RooRealProxy alphaL;
        RooRealProxy sigmaR;
        RooRealProxy alphaR;

        Double_t evaluate() const {
            constexpr double sqrtPiOver2 = 1.2533141373;
            constexpr double sqrt2 = 1.4142135624;

            double t = (x - x0) / (x < x0 ? sigmaL : sigmaR);
            double v = 0;
            if (t < -alphaL) {
                double a = 0.5 * alphaL * alphaL;
                double b = alphaL * t;
                v = std::exp(a + b);
            } else if (t <= alphaR) {
                v = std::exp(-0.5 * t * t);
            } else {
                double a = 0.5 * alphaR * alphaR;
                double b = alphaR * (-t);
                v = std::exp(a + b);
            }

            return v;
        }
    };



void RooFit_(TTree *sigtree, TTree *bkgtree, TTree *negtree, double &nsig, double &nsig_error, int max_iter = 10, double tol = 0.01) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooRealVar x("x", "x", -1.0, 3.0);
    int Ntot = negtree->GetEntries();
    // 信号与背景参数初始化
    RooRealVar sig_mean("sig_mean", "sig_mean", 0.8, 0.5, 1.1);
    RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 1, 0.05, 1.1);
    RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1, 0.01, 5.0);
    RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 1, 0.05, 1.1);
    RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1, 0.01, 6.0);
    RooRealVar bkg_mean("bkg_mean", "bkg_mean", 0.1, 0.001, 0.35);
    RooRealVar bkg_sigmaL("bkg_sigmaL", "bkg_sigmaL", 1, 0.05, 1.1);
    RooRealVar bkg_alphaL("bkg_alphaL", "bkg_alphaL", 1, 0.01, 5.0);
    RooRealVar bkg_sigmaR("bkg_sigmaR", "bkg_sigmaR", 1, 0.05, 1.1);
    RooRealVar bkg_alphaR("bkg_alphaR", "bkg_alphaR", 1, 0.01, 5.0);
    
    // 修改点1: 使用分数而不是绝对比例
    RooRealVar frac_sig("frac_sig", "signal fraction", 0.3, 0.0, 1.0);  // 信号比例参数

    // 数据读取（保持不变）
    Double_t pr_mass2_value, pi_mass2_value, neg_mass2_value;
    if (sigtree->GetBranch("pr_mass2")) {
        sigtree->SetBranchAddress("pr_mass2", &pr_mass2_value);
    }
    else{
        sigtree->SetBranchAddress("pr1_mass2", &pr_mass2_value);
    }
    if (bkgtree->GetBranch("pi_mass2")) {
        bkgtree->SetBranchAddress("pi_mass2", &pi_mass2_value);
    }
    else{
        bkgtree->SetBranchAddress("pi1_mass2", &pi_mass2_value);
    }
    if (negtree->GetBranch("neg_mass2")) {
        negtree->SetBranchAddress("neg_mass2", &neg_mass2_value);
    } else if (negtree->GetBranch("neg1_mass2")) {
        negtree->SetBranchAddress("neg1_mass2", &neg_mass2_value);
    }

    RooDataSet sigData("sigData", "Signal Data", RooArgSet(x));
    RooDataSet bkgData("bkgData", "Background Data", RooArgSet(x));
    RooDataSet negData("negData", "Negative Data", RooArgSet(x));

    for (Long64_t i = 0; i < sigtree->GetEntries(); i++) {
        sigtree->GetEntry(i);
        x.setVal(pr_mass2_value);
        sigData.add(RooArgSet(x));
    }
    for (Long64_t i = 0; i < bkgtree->GetEntries(); i++) {
        bkgtree->GetEntry(i);
        x.setVal(pi_mass2_value);
        bkgData.add(RooArgSet(x));
    }
    for (Long64_t i = 0; i < negtree->GetEntries(); i++) {
        negtree->GetEntry(i);
        x.setVal(neg_mass2_value);
        negData.add(RooArgSet(x));
    }

    // 构建模型
    expgausexp f_sig("sigPdf", "Signal PDF", x, sig_mean, sig_sigmaL, sig_alphaL, sig_sigmaR, sig_alphaR);
    expgausexp f_bkg("bkgPdf", "Background PDF", x, bkg_mean, bkg_sigmaL, bkg_alphaL, bkg_sigmaR, bkg_alphaR);

    // 第一步：分别拟合 signal 和 background（保持不变）
    f_sig.fitTo(sigData, 
        RooFit::Save(), 
        RooFit::PrintLevel(-1), 
        RooFit::Verbose(kFALSE),
        RooFit::Minimizer("Minuit2", "migrad")
    );

    f_bkg.fitTo(bkgData, 
        RooFit::Save(), 
        RooFit::PrintLevel(-1), 
        RooFit::Verbose(kFALSE),
        RooFit::Minimizer("Minuit2", "migrad")
    );
    
    // 固定信号/背景形状参数
    for (auto* var : {&sig_mean, &sig_sigmaL, &sig_alphaL, &sig_sigmaR, &sig_alphaR,
                      &bkg_mean, &bkg_sigmaL, &bkg_alphaL, &bkg_sigmaR, &bkg_alphaR}) {
        var->setConstant(true);
    }

    // 修改点2: 使用分数构建组合PDF
    // 信号比例：frac_sig
    // 背景比例：1 - frac_sig
    RooAddPdf combinedPdf("combinedPdf", "combinedPdf", 
                          RooArgList(f_sig, f_bkg), 
                          RooArgList(frac_sig)); // 只有信号比例是自由参数

    // 第二步：对 negData 做组合拟合



    RooFitResult* fitResult = combinedPdf.fitTo(negData, 
        RooFit::Save(true), 
        RooFit::PrintLevel(-1), 
        RooFit::Verbose(kFALSE),
        RooFit::Minimizer("Minuit2", "migrad")
    );

    // 检查拟合状态
    if (fitResult && fitResult->status() == 0) {
        // 修改点3: 获取信号比例及其误差
        nsig = frac_sig.getVal()*Ntot;
        nsig_error = frac_sig.getError()*Ntot;
    } else {
        std::cerr << "Fit failed!" << std::endl;
        nsig = -1;
        nsig_error = -1;
    }
    
    // 删除动态分配的对象（如果之前创建了）
    // 注意：frac_sig 是自动变量，不需要手动删除
}

void expgausexp_roofit(int seq, int npi, int low, int high){
        // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
    // TFile *f = new TFile(Form("/data03/tianye/pbar/root/template_bin25/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/10percent/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/sys_result/unbinned_roofit_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");

    std::vector<double> nsig(8);
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    fitParamsTree->Branch("nsig", &nsig);
    for(int j=low;j<high;j++){
        TTree *tsigtree=(TTree *)f->Get(Form("sigtree_%04d",j));
        TTree *tbkgtree=(TTree *)f->Get(Form("bkgtree_%04d",j));
        TTree *tsigtree1=(TTree *)f->Get(Form("sigtree1_%04d",j));
        TTree *tbkgtree1=(TTree *)f->Get(Form("bkgtree1_%04d",j));
        // TTree *bkgtree1=(TTree *)f->Get(Form("bkgtree1_%04d",j));
        // TTree *bkgtree2=(TTree *)f->Get(Form("bkgtree2_%04d",j));
        // TTree *bkgtree3=(TTree *)f->Get(Form("bkgtree3_%04d",j));
        // TTree *bkgtree4=(TTree *)f->Get(Form("bkgtree4_%04d",j));
        TTree *tnegtree=(TTree *)f->Get(Form("negtree_%04d",j));
        TTree *tnegtree1=(TTree *)f->Get(Form("negtree1_%04d",j));

        RooFit_(tsigtree, tbkgtree, tnegtree, nsig[0], nsig[4]);
        RooFit_(tsigtree, tbkgtree, tnegtree1, nsig[1], nsig[5]);
        RooFit_(tsigtree1, tbkgtree1, tnegtree, nsig[2], nsig[6]);
        RooFit_(tsigtree1, tbkgtree1, tnegtree1, nsig[3], nsig[7]);

    // RooFit_(sigtree,bkgtree1,negtree,nsig[1],nsig[6]);
    // RooFit_(sigtree,bkgtree2,negtree,nsig[2],nsig[7]);
    // RooFit_(sigtree,bkgtree3,negtree,nsig[3],nsig[8]);
    // RooFit_(sigtree,bkgtree4,negtree,nsig[4],nsig[9]);


    fitParamsTree->Fill();


    fOutput->cd();

   std::cout <<"j="<< j << std::endl;
        std::cout <<"sig_yields"<< nsig[0] << "sig_yields"<<nsig[1]<<std::endl;
        std::cout <<"sig_yields"<< nsig[2] << "sig_yields"<<nsig[3]<<std::endl;
        //std::cout <<"sig_error"<< nsig[3] << std::endl;

    }
    fOutput->cd();
    fitParamsTree->Write();
    fitParamsTree->Delete();
    fOutput->Close();
    
    std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/RooFit_Unbinned_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "RooFit_Unbinned - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "RooFit_Unbinned - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
        timeFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

int main(int argc,char* argv[]){
    expgausexp_roofit(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    std::cout<<"Finished!"<<std::endl;
    return 0;
}