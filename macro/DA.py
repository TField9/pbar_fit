import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import Template
import uproot
import awkward as ak
import ROOT  # 添加这一行
from ROOT import Form  # 添加这一行

def RunQOptimization(sig_hist, bkg_hist, data_hist,methods="da"):

    sig_array = sig_hist.values()
    bkg_array = bkg_hist.values()
    data_array = data_hist.values()
    

    xe = sig_hist.axes[0].edges()
    
    cost = Template(data_array, xe, [sig_array, bkg_array],method=methods)
    m = Minuit(cost, 200, 200)
    m.limits = [(0, 500), (500, 1000)]
    m.migrad()
    
    # 获取拟合结果
    params = m.values
    errors = m.errors
    
    # 获取拟合模板和标准差
    y, yerr = cost.prediction(params)
    
    # 获取拉距值（pulls）
    pulls = cost.pulls(params)
    
    return params, errors, y, yerr, pulls

# 主程序
def main():
    # 打开 ROOT 文件
    f = uproot.open("/data03/tianye/pbar/root/paper_data/paper_data.root")
    
    # 创建 TTree 和 TH1D
    f1 = ROOT.TFile("/data03/tianye/pbar/root/paper_data/fit_results.root", "RECREATE")
    c2 = ROOT.TCanvas("c2", "c2", 1200, 900)
    c2.Print("/data03/tianye/pbar/root/paper_data/paper_pull_comparison.pdf[")
    for k in range(5):
        t = ROOT.TTree(f"fit_results_{k}", "Fit Results")
        pull_hist = ROOT.TH1D("pull_hist", "Pull Distribution1", 30, -6, 6)
        pull_hist1 = ROOT.TH1D("pull_hist1", "Pull Distribution2", 30, -6, 6)
        pull_hist2 = ROOT.TH1D("pull_hist2", "Pull Distribution3", 30, -6, 6)
        
        params = np.zeros(2)
        params1= np.zeros(2)
        params2= np.zeros(2)
        pull = np.zeros(1)
        pull1 = np.zeros(1)
        pull2 = np.zeros(1)
        t.Branch("fitparams", params, "fitparams[2]/D")
        t.Branch("fitparams1", params1, "fitparams1[2]/D")
        t.Branch("fitparams2", params2, "fitparams2[2]/D")
        # t.Branch("pull", pull, "pull/D")
        # t.Branch("pull1", pull1, "pull1/D")
        # t.Branch("pull2", pull2, "pull2/D")

        
        for i in range(2000):
            sig_hist = f[f"sig_hist_seed{k}_nmc100_{i}"]
            bkg_hist = f[f"bkg_hist_seed{k}_nmc100_{i}"]
            data_hist = f[f"data_hist_seed{k}_nmc100_{i}"]
            
            params, errors, y, yerr, pulls = RunQOptimization(sig_hist, bkg_hist, data_hist,"da")
            params1, errors1, y1, yerr1, pulls1 = RunQOptimization(sig_hist, bkg_hist, data_hist,methods="jsc")
            params2, errors2, y2, yerr2, pulls2 = RunQOptimization(sig_hist, bkg_hist, data_hist,methods="asy")

            t.Fill()
            pull_hist.Fill((params[0]-250)/errors[0])
            pull_hist1.Fill((params1[0]-250)/errors1[0])
            pull_hist2.Fill((params2[0]-250)/errors2[0])
            
            if i % 100 == 0:
                print(f"Processed {i} histograms")
                print(params)

        
            # 绘制 Pull 分布
            # c1 = ROOT.TCanvas("c1", "c1", 800, 600)
            # pull_hist.Draw()
            # c1.Print("/data03/tianye/pbar/root/paper_data/paper_pull.pdf")
            # 创建一个画布并分割成四个pad
        c2 = ROOT.TCanvas("c2", "c2", 1200, 900)
        c2.Divide(2, 2)
        
        # 第一个pad画三个method的pull
        # ROOT.gStyle.SetOptStat(0)
        pull_hist.GetYaxis().SetRangeUser(0, 700)
        pull_hist1.GetYaxis().SetRangeUser(0, 700)
        pull_hist2.GetYaxis().SetRangeUser(0, 700)
        c2.cd(1)
        pull_hist.SetLineColor(ROOT.kRed)
        pull_hist1.SetLineColor(ROOT.kBlue)
        pull_hist2.SetLineColor(ROOT.kGreen)
        pull_hist.Draw()
        pull_hist1.Draw("SAME")
        pull_hist2.Draw("SAME")
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        legend.AddEntry(pull_hist, "DA Method", "l")
        legend.AddEntry(pull_hist1, "JSC Method", "l")
        legend.AddEntry(pull_hist2, "ASY Method", "l")
        legend.Draw()
        
        # 第二个pad画DA method的pull
        c2.cd(2)
        pull_hist.Draw()
        pull_hist.SetTitle("DA Method Pull Distribution")
        
        # 第三个pad画JSC method的pull
        c2.cd(3)
        pull_hist1.Draw()
        pull_hist1.SetTitle("JSC Method Pull Distribution")
        
        # 第四个pad画ASY method的pull
        c2.cd(4)
        pull_hist2.Draw()
        pull_hist2.SetTitle("ASY Method Pull Distribution")
        
        # 保存画布
        c2.Print("/data03/tianye/pbar/root/paper_data/paper_pull_comparison.pdf")
        
        # 保存结果
        f1.cd()
        t.Write()
        pull_hist.Delete()
        pull_hist1.Delete()
        pull_hist2.Delete()
        print(k)

    c2.Print("/data03/tianye/pbar/root/paper_data/paper_pull_comparison.pdf]")
    f1.Close()
    

if __name__ == "__main__":
    main()