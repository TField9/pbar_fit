import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import Template
import uproot
import awkward as ak
import ROOT  # 添加这一行
from ROOT import Form  # 添加这一行
import sys
import ROOT
from ROOT import TFile, TTree

from iminuit import Minuit



def RunQOptimization(sig_hist, bkg_hist, data_hist, methods="da"):


    sig_array = sig_hist.values()
    bkg_array = bkg_hist.values()
    data_array = data_hist.values()
    

    xe = sig_hist.axes[0].edges()

    # 构造 Template 对象
    cost = Template(
        data_array,
        xe,
        [sig_array, bkg_array],
        method=methods
    )

    # 初始化 Minuit (假设两个可调参数初值为 200, 200)
    m = Minuit(cost, 200, 200)
    # 设置两个参数的取值范围
    m.limits = [(0, 500), (0,500)]
    # 运行拟合
    m.migrad()

    # 获取拟合结果
    params = m.values
    errors = m.errors
    # （如需查看 pulls、y, yerr 等，也可在此获取）
    return np.array([params[0], params[1],errors[0], errors[1]])


def SPD(seq, npi, low, high):
    # 计时器
    timer = ROOT.TStopwatch()
    timer.Start()

    # 日志等级
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    # 打开输入 ROOT 文件 (READ 模式)
    # 注意这里加了 'f' 以正确插入 seq、npi
    # f_input = uproot.open(f"/data03/tianye/pbar/root/template_multi_bin/toymc_bin1_seq{seq:04d}_npi{npi:04d}.root")
    f_input = uproot.open(f"/data03/tianye/pbar/root/template_morebin/toymc_bin1_seq{seq:04d}_npi{npi:04d}.root")

    # 打开输出 ROOT 文件 (RECREATE)
    f_output = ROOT.TFile(f"/data03/tianye/pbar/root/multi_bin/binned_SPD_seq{seq:04d}_npi{npi:04d}_{low:04d}.root", "RECREATE")

    fitParamsTree = TTree("fitParamsTree", "Fit Parameters Tree")

    from array import array
    params1 = array('f', [0., 0., 0., 0.])
    params2 = array('f', [0., 0., 0., 0.])

    fitParamsTree.Branch("bin30_SPD_pararr", params1, "bin30_SPD_pararr[4]/F")
    fitParamsTree.Branch("bin100_SPD_pararr", params2, "bin100_SPD_pararr[4]/F")

    for j in range(low, high):
        data_pr25 = f_input[f"data_pr30_{j:04d}"]
        data_pr100= f_input[f"data_pr100_{j:04d}"]
        data_pi25 = f_input[f"data_pi30_{j:04d}"]
        data_pi100= f_input[f"data_pi100_{j:04d}"]
        data_neg25= f_input[f"data_neg30_{j:04d}"]
        data_neg100= f_input[f"data_neg100_{j:04d}"]

        res1 = RunQOptimization(data_pr25, data_pi25, data_neg25)
        res2 = RunQOptimization(data_pr100, data_pi100, data_neg100)

        for idx, val in enumerate(res1):
            params1[idx] = val
        for idx, val in enumerate(res2):
            params2[idx] = val

        fitParamsTree.Fill()

        print(f"j = {j}, res2 = {res2}")

    f_output.cd()
    fitParamsTree.Write()
    f_output.Close()

    timer.Stop()
    # 如需记录时间，再自行追加
    # ...



def main():
    """
    Python main() 函数，对应原来的 C++ main(int argc, char* argv[])
    """
    seq = int(sys.argv[1])
    npi = int(sys.argv[2])
    low = int(sys.argv[3])
    high = int(sys.argv[4])

    SPD(seq, npi, low, high)


if __name__ == "__main__":
    main()
