#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include <vector>
#include <iostream>

class QMinimizer : public ROOT::Math::IBaseFunctionMultiDim {
public:
    QMinimizer(TH1* signal, TH1* background, TH1* observed)
        : sig(dynamic_cast<TH1*>(signal->Clone("signal_clone"))),
          bkg(dynamic_cast<TH1*>(background->Clone("background_clone"))),
          data(observed),
          y1_opt(0), y2_opt(0), y1_err(0), y2_err(0) {}

    ~QMinimizer() override {
        delete sig;
        delete bkg;
    }

    // 克隆函数
    QMinimizer* Clone() const override {
        return new QMinimizer(sig, bkg, data);
    }

    // 维度数量：2
    unsigned int NDim() const override {
        return 2;
    }

    // 核心评估函数
    double DoEval(const double* p) const override {
        double y1 = p[0];
        double y2 = p[1];
        return calculateQ(y1, y2);
    }

    // 优化接口：运行优化并返回结果
    void Minimize() {
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimizer->SetFunction(*this);
        minimizer->SetMaxFunctionCalls(10000);
        minimizer->SetMaxIterations(1000);
        minimizer->SetTolerance(0.001);

        // 设置参数初始值及边界
        minimizer->SetLimitedVariable(0, "y1", 200, 0.01, 0.0, 500); // 参数 y1
        minimizer->SetLimitedVariable(1, "y2", 700, 0.01, 500, 1000); // 参数 y2

        // 执行最小化
        minimizer->Minimize();

        // 提取优化结果
        const double* results = minimizer->X();
        const double* errors = minimizer->Errors();

        y1_opt = results[0];
        y2_opt = results[1];
        y1_err = errors[0];
        y2_err = errors[1];

        delete minimizer;
    }

    // 获取结果
    double GetY1() const { return y1_opt; }
    double GetY2() const { return y2_opt; }
    double GetY1Error() const { return y1_err; }
    double GetY2Error() const { return y2_err; }

private:
    TH1* sig;
    TH1* bkg;
    TH1* data;
    double y1_opt, y2_opt; // 最优参数
    double y1_err, y2_err; // 参数误差

    // Q 计算逻辑
    double calculateQ(double y1, double y2) const {
           double Q = 0;
        for (int i = 0; i < sig->GetNbinsX(); i++) {
            double n = data->GetBinContent(i + 1);
            double Vn = data->GetBinContent(i + 1);
            double t = 1;
            double mu0 = y1 * sig->GetBinContent(i + 1) / sig->GetSumOfWeights() +
                         y2 * bkg->GetBinContent(i + 1) / bkg->GetSumOfWeights();

            if(mu0!=0){

                double Vmu = y1 * y1 * sig->GetBinContent(i + 1) / (sig->GetSumOfWeights()*sig->GetSumOfWeights()) +
                            y2 * y2 * bkg->GetBinContent(i + 1) / (bkg->GetSumOfWeights()*bkg->GetSumOfWeights());
                double s = (Vmu != 0) ? mu0 / Vmu : 0;
                double beta=(t * n + s * mu0) / (t * mu0 + s * mu0);
                double cash1 = cash(t * n, beta * t * mu0);
                double cash2 = cash(s * mu0, beta * s * mu0);
                Q += cash1 + cash2;
            }
            else{
                Q+=0;
            }
        } 
        return Q;
    }

    // Cash 统计量计算
    double cash(double observed, double expected) const {
        if (expected <= 0) return 0;
        if (observed <= 0) return 0;
        return expected - observed + observed * log(observed / expected);
    }
};

std::vector <double> RunQOptimization(TH1* signal, TH1* background, TH1* data) {
    QMinimizer qMinimizer(signal, background, data);
    qMinimizer.Minimize();

    // // 打印优化结果
    // std::cout << "Optimal y1: " << qMinimizer.GetY1()
    //           << " ± " << qMinimizer.GetY1Error() << std::endl;
    // std::cout << "Optimal y2: " << qMinimizer.GetY2()
    //           << " ± " << qMinimizer.GetY2Error() << std::endl;
    std::vector<double> par(4);
    par[0] = qMinimizer.GetY1();
    par[1] = qMinimizer.GetY2();
    par[2] = qMinimizer.GetY1Error();
    par[3] = qMinimizer.GetY2Error();
    return par;
}

void paper_fit(){

TFile *f=new TFile("/data03/tianye/pbar/root/paper_data/paper_data.root","read");
TFile *f1=new TFile("/data03/tianye/pbar/root/paper_data/paper_pull.root","recreate");
TTree *t = new TTree("fitParamsTree", "fitParamsTree");
TH1D* pull_hist= new TH1D("pull_hist", "pull_hist", 30, -6, 6);
std::vector<double> params(4);
double pull;
t->Branch("fitparams", &params);
t->Branch("pull", &pull);

for(int i=0;i<2000;i++){
    TH1D* sig_hist = (TH1D*)f->Get(Form("sig_hist_100_%d",i));
    TH1D* bkg_hist = (TH1D*)f->Get(Form("bkg_hist_100_%d",i));
    TH1D* data_hist = (TH1D*)f->Get(Form("data_hist_100_%d",i));
    params = RunQOptimization(sig_hist, bkg_hist, data_hist);
    pull=(params[0]-250)/params[2];
    pull_hist->Fill(pull);
    t->Fill();
    delete sig_hist;
    delete bkg_hist;
    delete data_hist;
    if(i%100==0){
        std::cout<<i<<std::endl;
    }
}

TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
pull_hist->Draw();
c1->Print("/data03/tianye/pbar/root/paper_data/paper_pull.pdf");
f1->cd();
t->Write();
pull_hist->Write();
f1->Close();

}