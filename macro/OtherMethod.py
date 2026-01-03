# #!/usr/bin/env python3
# import sys, time
# import numpy as np
# from array import array
# import ROOT
# from iminuit import Minuit
# from iminuit.cost import template_chi2_da, template_chi2_jsc, template_nll_asy

# def hist_to_array(h):
#     """将 TH1D 转为 numpy 数组并返回原始计数及其 Poisson 方差"""
#     nb   = h.GetNbinsX()
#     vals = np.array([h.GetBinContent(i+1) for i in range(nb)], dtype=float)
#     var  = vals.copy()
#     return vals, var

# def make_cost(n_obs, T_sig, var_sig, T_bkg, var_bkg, cost_func):
#     """构造 iminuit cost 函数，参数名为 y_sig，背景事件数自动计算为 Ntot - y_sig"""
#     Ntot = n_obs.sum()
    
#     def cost(y_sig):
#         y_bkg = Ntot - y_sig
#         mu = y_sig * T_sig + y_bkg * T_bkg
#         if np.any(mu <= 0):
#             return np.inf
#         mu_var = y_sig**2 * var_sig + y_bkg**2 * var_bkg
#         return cost_func(n_obs, mu, mu_var)
    
#     return cost

# def fit_nbs(n_obs, T_sig, var_sig, T_bkg, var_bkg):
#     """
#     对给定观测数据执行三种拟合方法，返回每种方法的最终 y_sig 值列表。
#     """
#     # 定义三种拟合方法及其对应的误差定义
#     cost_fns = [
#         (template_chi2_da, Minuit.LEAST_SQUARES),
#         (template_chi2_jsc, Minuit.LEAST_SQUARES),
#         (template_nll_asy, Minuit.LIKELIHOOD),
#     ]
    
#     results = []  # 存储每种方法的拟合结果
#     Ntot = n_obs.sum()  # 总事件数
    

    
#     # 遍历每种拟合方法
#     for cost_fn, errdef in cost_fns:
#         # 创建成本函数 - 注意这里只使用一个参数 y_sig
#         cost_func = make_cost(n_obs, T_sig, var_sig, T_bkg, var_bkg, cost_fn)
        
#         # 根据总事件数设置初始值 - 只设置信号初始值
#         if Ntot < 300:
#             init_sig = Ntot * 0.5
#         else:
#             init_sig = Ntot * 0.5
        
#         # 创建并配置 Minuit 优化器 - 只有 y_sig 一个参数
#         m = Minuit(cost_func, y_sig=init_sig)
#         m.limits['y_sig'] = (0, Ntot)
#         m.strategy = 2
#         m.tol = 0.01
#         m.errordef = errdef
#         m.errors['y_sig'] = max(0.01 * init_sig, 0.5)  # 确保步长不过小
        
#         # 执行优化和误差计算
#         m.migrad(ncall=10000, iterate=5)
#         m.hesse()
        
#         # 获取信号事件数
#         y_sig_cur = m.values['y_sig']
#         results.append(y_sig_cur)
        
#         # 输出当前方法的结果
#         method_name = cost_fn.__name__
#         # 计算背景事件数
#         y_bkg_cur = Ntot - y_sig_cur
#         #print(f"{method_name} | y_sig = {y_sig_cur:.2f}, y_bkg = {y_bkg_cur:.2f}, 总事件数 = {y_sig_cur + y_bkg_cur:.2f}")
    
#     # 返回所有方法的结果
#     return results

# def fit_start(npr,npi):
#     t0 = time.time()

#     fin_name  = f"/data03/tianye/pbar/root/both_vary//toymc_bin1_seq0001_npi{npr:04d}.root"
#     fin_name1 = f"/data03/tianye/pbar/root/both_vary//toymc_bin1_seq0001_npi{npi:04d}.root"
#     fout_name = f"/data03/tianye/pbar/root/2d_result/binned_OM_seq0001_npr{npr:04d}_npi{npi:04d}.root"

#     fin  = ROOT.TFile.Open(fin_name,  "READ")
#     fin1 = ROOT.TFile.Open(fin_name1,  "READ")
#     fout = ROOT.TFile.Open(fout_name, "RECREATE")
#     tree = ROOT.TTree("fitParamsTree", "Fit Parameters Tree")

#     # 定义数组：第一套模板结果
#     arr_DA  = array('d', [0.0, 0.0])
#     arr_JSC = array('d', [0.0, 0.0])
#     arr_asy = array('d', [0.0, 0.0])


#     # 创建分支
#     tree.Branch("Ysig_DA",   arr_DA,  "Ysig_DA[2]/D")
#     tree.Branch("Ysig_JSC",  arr_JSC, "Ysig_JSC[2]/D")
#     tree.Branch("Ysig_asy",  arr_asy, "Ysig_asy[2]/D")


#     binname = "30"
#     total = 10000

#     for j in range(0, 10000):
#         # 读直方图
#         hsig  = fin.Get(f"data_pr{binname}_{j:04d}")
#         hbkg  = fin1.Get(f"data_pi{binname}_{j:04d}")
#         hneg  = fin1.Get(f"data_neg{binname}_{j:04d}")

#         if not all([hsig, hbkg, hneg]):
#             print(f"[Warning] j={j} 有缺失直方图，跳过")
#             continue

#         # 转 numpy
#         raw_sig, var_sig_c = hist_to_array(hsig)
#         raw_bkg, var_bkg_c = hist_to_array(hbkg)
#         n_obs, _   = hist_to_array(hneg)
#         # 构造模板 0
#         eps = 1e-8
#         rs   = np.maximum(raw_sig, eps);   rb   = np.maximum(raw_bkg, eps)
#         ds   = rs.sum();                    db   = rb.sum()
#         T_sig = rs / ds;                    T_bkg = rb / db
#         var_sig = var_sig_c / (ds**2);      var_bkg = var_bkg_c / (db**2)

#         # 防止方差为0
#         epsv = 1e-8
#         var_sig  = np.maximum(var_sig,  epsv)
#         var_bkg  = np.maximum(var_bkg,  epsv)


#         # 四次拟合
#         r01 = fit_nbs(n_obs,  T_sig,  var_sig,  T_bkg,  var_bkg)
 

#         # 填入数组
#         arr_DA[0], arr_JSC[0], arr_asy[0]   = r01


#         tree.Fill()
#         if j%100==0:
#             print(f"[nsig={r01}]")



#     fout.cd()
#     tree.Write()
#     fout.Close()
#     fin.Close()

#     total_t = time.time() - t0
#     print(f"全部完成: {total} fits 用时 {total_t:.1f}s, 平均 {total_t/total:.3f}s/fit")

# if __name__ == "__main__":
#     npr,npi = map(int, sys.argv[1:])
#     fit_start(npr,npi)

#!/usr/bin/env python3
import sys, time
import numpy as np
from array import array
import ROOT
from iminuit import Minuit
from iminuit.cost import template_chi2_da, template_chi2_jsc, template_nll_asy

def hist_to_array(h):
    """将 TH1D 转为 numpy 数组并返回原始计数及其 Poisson 方差"""
    nb   = h.GetNbinsX()
    vals = np.array([h.GetBinContent(i+1) for i in range(nb)], dtype=float)
    var  = vals.copy()
    return vals, var

def make_cost(n_obs, T_sig, var_sig, T_bkg, var_bkg, cost_func):
    """构造 iminuit cost 函数，参数名为 y_sig，背景事件数自动计算为 Ntot - y_sig"""
    Ntot = n_obs.sum()
    def cost(y_sig):
        y_bkg = Ntot - y_sig
        mu = y_sig * T_sig + y_bkg * T_bkg
        if np.any(mu <= 0):
            return np.inf
        mu_var = y_sig**2 * var_sig + y_bkg**2 * var_bkg
        return cost_func(n_obs, mu, mu_var)
    return cost

def fit_nbs(n_obs, T_sig, var_sig, T_bkg, var_bkg):
    """
    对给定观测数据执行三种拟合方法，返回每种方法的最终 y_sig 值列表。
    """
    cost_fns = [
        (template_chi2_da, Minuit.LEAST_SQUARES),
        (template_chi2_jsc, Minuit.LEAST_SQUARES),
        (template_nll_asy, Minuit.LIKELIHOOD),
    ]
    results = []
    Ntot = n_obs.sum()

    for cost_fn, errdef in cost_fns:
        cost_func = make_cost(n_obs, T_sig, var_sig, T_bkg, var_bkg, cost_fn)
        init_sig = Ntot * 0.5
        m = Minuit(cost_func, y_sig=init_sig)
        m.limits['y_sig'] = (0, Ntot)
        m.strategy = 2
        m.tol = 0.1
        m.errordef = errdef
        m.errors['y_sig'] = 0.1
        m.migrad(ncall=100000, iterate=5)
        m.hesse()
        results.append(m.values['y_sig'])
    return results

def fit_start(npr, npi):
    t0 = time.time()

    fin_name  = f"/data03/tianye/pbar/root/both_vary//toymc_bin1_seq0001_npi{npr:04d}.root"
    fin_name1 = f"/data03/tianye/pbar/root/both_vary//toymc_bin1_seq0001_npi{npi:04d}.root"
    fout_name = f"/data03/tianye/pbar/root/2d_result/binned_OM_seq0001_npr{npr:04d}_npi{npi:04d}.root"

    fin  = ROOT.TFile.Open(fin_name,  "READ")
    fin1 = ROOT.TFile.Open(fin_name1, "READ")
    fout = ROOT.TFile.Open(fout_name, "RECREATE")
    tree = ROOT.TTree("fitParamsTree", "Fit Parameters Tree")

    arr_DA  = array('d', [0.0, 0.0])
    arr_JSC = array('d', [0.0, 0.0])
    arr_asy = array('d', [0.0, 0.0])

    tree.Branch("Ysig_DA",   arr_DA,  "Ysig_DA[2]/D")
    tree.Branch("Ysig_JSC",  arr_JSC, "Ysig_JSC[2]/D")
    tree.Branch("Ysig_asy",  arr_asy, "Ysig_asy[2]/D")

    total = 10000
    for j in range(total):
        hsig  = fin.Get(f"data_pr30_{j:04d}")
        hbkg  = fin1.Get(f"data_pi30_{j:04d}")
        hneg  = fin1.Get(f"data_neg30_{j:04d}")
        if not all([hsig, hbkg, hneg]):
            print(f"[Warning] j={j} 有缺失直方图，跳过")
            continue

        raw_sig, var_sig_c = hist_to_array(hsig)
        raw_bkg, var_bkg_c = hist_to_array(hbkg)
        n_obs, _          = hist_to_array(hneg)

        # 构造模板并用二项分布计算误差
        eps = 1e-12
        rs = np.maximum(raw_sig, eps)
        rb = np.maximum(raw_bkg, eps)
        ds = rs.sum()
        db = rb.sum()
        T_sig = rs / ds
        T_bkg = rb / db

        # 二项分布方差：Var(T_i) = T_i*(1-T_i)/N
        var_sig = T_sig * (1.0 - T_sig) / ds
        var_bkg = T_bkg * (1.0 - T_bkg) / db

        # 防止方差为 0
        epsv = 1e-12
        var_sig = np.maximum(var_sig, epsv)
        var_bkg = np.maximum(var_bkg, epsv)

        # 执行拟合
        r01 = fit_nbs(n_obs, T_sig, var_sig, T_bkg, var_bkg)

        arr_DA[0], arr_JSC[0], arr_asy[0] = r01
        tree.Fill()
        if j % 100 == 0:
            print(f"[nsig={r01}]", flush=True)

    fout.cd()
    tree.Write()
    fout.Close()
    fin.Close()
    fin1.Close()

    total_t = time.time() - t0
    print(f"全部完成: {total} fits 用时 {total_t:.1f}s, 平均 {total_t/total:.3f}s/fit")

if __name__ == "__main__":
    npr, npi = map(int, sys.argv[1:])
    fit_start(npr, npi)
