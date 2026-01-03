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
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


class QMinimizer : public ROOT::Math::IBaseFunctionMultiDim {
public:
    // 构造时仍然传入三幅直方图
    QMinimizer(TH1* signal, TH1* background, TH1* observed)
        : sig(dynamic_cast<TH1*>(signal->Clone("signal_clone"))),
          bkg(dynamic_cast<TH1*>(background->Clone("background_clone"))),
          data(observed),
          a_opt(0), a_err(0) ,Ntot(data->GetEntries())
    {}

    ~QMinimizer() override {
        delete sig;
        delete bkg;
    }

    // 克隆接口
    QMinimizer* Clone() const override {
        return new QMinimizer(sig, bkg, data);
    }

    // 维度改为 1（只优化参数 a）
    unsigned int NDim() const override {
        return 1;
    }

    // 评价函数：p[0] 对应 a，y1 = a, y2 = 1 - a
    double DoEval(const double* p) const override {
        double a = p[0];
        double y1 = a;
        double y2 = 1.0 - a;
        // 为了与原来逻辑保持一致，需要把 y1, y2 转回“计数”尺度：
        // 原代码里，y1, y2 是 count*400，只有在 GetY1()、GetY2() 时才除以 400。
        // 这里我们直接把比例 a 映射到 count 空间：令 N_total = 400，
        // 则 count1 = a*N_total, count2 = (1-a)*N_total。
        
        double y1_count = a * Ntot;
        double y2_count = (1.0 - a) * Ntot;
        return calculateQ(y1_count, y2_count);
    }

    // 执行最小化
    void Minimize() {
        // 用 Minuit2 Migrad
        auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetFunction(*this);
        minimizer->SetMaxFunctionCalls(1000000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.001);

        // 只设置一个游程参数 a，范围 [0,1]
        if(Ntot==400){
            minimizer->SetLimitedVariable(0, "a", 0.4, 0.0001, 0.0, 1.0);
        }
        else{
            minimizer->SetLimitedVariable(0, "a", 0.2, 0.0001, 0.0, 1.0);
        }
    

        minimizer->Minimize();

        const double* xopt = minimizer->X();
        const double* xerr = minimizer->Errors();
        a_opt = xopt[0];
        a_err = xerr[0];
        delete minimizer;
    }

    // 外部获取 y1, y2 及误差（都除以400，恢复原来 count 单位）：
    double GetY1() const { return a_opt * Ntot ;}
    double GetY2() const { return (1.0 - a_opt) * Ntot ; }
    double GetY1Error() const { return a_err * Ntot ; }
    double GetY2Error() const { return a_err * Ntot ; }

private:
    TH1* sig;
    TH1* bkg;
    TH1* data;
    double a_opt, a_err;
    double Ntot;

    // 原来的 calculateQ 函数不变，只是输入参数不再是真正的 y1/y2，
    // 而是我们在 DoEval 中传入的 “计数值”。
    double calculateQ(double y1, double y2) const {
        double Q = 0;
        double sumSig = sig->GetSumOfWeights();
        double sumBkg = bkg->GetSumOfWeights();

        for(int i = 1; i <= sig->GetNbinsX(); ++i) {
            double n = data->GetBinContent(i);
            double t = 1.0;
            // 计算 mu0 = y1 * sig_i/sumSig + y2 * bkg_i/sumBkg
            double mu0 = 0;
            if(sumSig > 0)  mu0 += y1 * (sig->GetBinContent(i) / sumSig);
            if(sumBkg > 0) mu0 += y2 * (bkg->GetBinContent(i) / sumBkg);

            if(mu0 <= 0) continue;

            // 计算 Vmu = y1^2 * sig_i/sumSig^2 + y2^2 * bkg_i/sumBkg^2
            double Vmu = 0;
            if(sumSig > 0)  Vmu += (y1*y1) * (sig->GetBinContent(i) / (sumSig*sumSig));
            if(sumBkg > 0) Vmu += (y2*y2) * (bkg->GetBinContent(i) / (sumBkg*sumBkg));

            double s = (Vmu > 0 ? mu0 / Vmu : 0.0);
            double beta = (t * n + s * mu0) / (t * mu0 + s * mu0);

            double c1 = cash(t * n,    beta * t * mu0);
            double c2 = cash(s * mu0,  beta * s * mu0);
            Q += c1 + c2;
        }
        return Q;
    }

    double cash(double observed, double expected) const {
        if(expected <= 0 || observed <= 0) return 0.0;
        return expected - observed + observed * log(observed / expected);
    }
};

// 外部调用接口保持不变，只是拿到 y1,y1Error 时可以直接呼叫 GetY1()/GetY1Error()
void RunQOptimization(TH1* signal, TH1* background, TH1* data, double &nsig, double &nsig_error) {
    QMinimizer qMinimizer(signal, background, data);
    qMinimizer.Minimize();
    nsig = qMinimizer.GetY1();
    nsig_error = qMinimizer.GetY1Error();
}

void SPD(int seq, int npi, int low, int high){
        // 创建 TStopwatch 对象
    TStopwatch timer;
    // 启动计时器
    timer.Start();
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
    TFile *f = new TFile(Form("/data03/tianye/pbar/root/10percent/toymc_bin1_seq%04d_npi%04d.root",seq,npi),"read");
    TFile *fOutput = new TFile(Form("/data03/tianye/pbar/root/sys_result/binned_SPD_seq%04d_npi%04d_%04d.root",seq,npi,low),"recreate");
    
    TTree *fitParamsTree = new TTree("fitParamsTree", "Fit Parameters Tree");
    std::vector<double> nsig(4);
    fitParamsTree->Branch("nsig", &nsig);
 

    for(int j=low;j<high;j++){
    TH1D *data_pr30=(TH1D*)f->Get(Form("data_pr30_%04d",j));
    TH1D *data_pr100=(TH1D*)f->Get(Form("data_pr100_%04d",j));
    TH1D *data_pi30=(TH1D*)f->Get(Form("data_pi30_%04d",j));
    TH1D *data_pi100=(TH1D*)f->Get(Form("data_pi100_%04d",j));
    TH1D *data_neg30=(TH1D*)f->Get(Form("data_neg30_%04d",j));
    TH1D *data1_neg30=(TH1D*)f->Get(Form("data1_neg30_%04d",j));
    TH1D *data_neg100=(TH1D*)f->Get(Form("data_neg100_%04d",j));

    RunQOptimization(data_pr30, data_pi30, data_neg30, nsig[0], nsig[2]);
    RunQOptimization(data_pr30, data_pi30, data1_neg30, nsig[1], nsig[3]);
    // RunQOptimization(data_pr30, data_pi1_30, data_neg30, nsig[1], nsig[6]);
    // RunQOptimization(data_pr30, data_pi2_30, data_neg30, nsig[2], nsig[7]);
    // RunQOptimization(data_pr30, data_pi3_30, data_neg30, nsig[3], nsig[8]);
    // RunQOptimization(data_pr30, data_pi4_30, data_neg30, nsig[4], nsig[9]);

    fitParamsTree->Fill();

    std::cout<<"j = "<<j<<std::endl;
    }
    fOutput->cd();
    fitParamsTree->Write();
    //fitParamsTree->Delete();
    fOutput->Close();
    std::ofstream timeFile("/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/time/Root_SPD_time.txt", std::ios::app); // 以追加模式打开文件
    if (timeFile.is_open()) {
        timeFile << "Root_SPD - Real time: " << timer.RealTime() << " seconds" << std::endl;
        timeFile << "Root_SPD - CPU time: " << timer.CpuTime() << " seconds" << std::endl;
        timeFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
   
}

int main(int argc,char* argv[]){
    SPD(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    return 0;
}



