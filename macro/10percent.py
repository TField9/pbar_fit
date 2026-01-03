import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from matplotlib.backends.backend_pdf import PdfPages  # For multi-page PDF saving

# Define the expgausexp function
def expgausexp(x, x0, sigmaL, alphaL, sigmaR, alphaR):
    is_scalar = np.isscalar(x)
    if is_scalar:
        x = np.array([x])
    else:
        x = np.asarray(x)
    
    t = np.zeros_like(x)
    v = np.zeros_like(x)
    
    left_mask  = (x < x0)
    right_mask = (x >= x0)
    t[left_mask]  = (x[left_mask]  - x0) / sigmaL
    t[right_mask] = (x[right_mask] - x0) / sigmaR
    
    mask_left_exp  = (t < -alphaL)
    mask_center    = (~mask_left_exp) & (t <= alphaR)
    mask_right_exp = (t > alphaR)
    
    if np.any(mask_left_exp):
        a_left = 0.5 * alphaL * alphaL
        v[mask_left_exp] = np.exp(a_left + alphaL * t[mask_left_exp])
    if np.any(mask_center):
        v[mask_center] = np.exp(-0.5 * t[mask_center] ** 2)
    if np.any(mask_right_exp):
        a_right = 0.5 * alphaR * alphaR
        v[mask_right_exp] = np.exp(a_right - alphaR * t[mask_right_exp])
    
    return v[0] if is_scalar else v

# Original background parameters
orig_bg_params = [0.001, 0.090, 0.498, 0.115, 0.530]

# Four sets of modified parameters
mod_params_list = [
    [0.001, 0.090, 0.498, 0.115, 0.400],  # func_bkg1
    [0.001, 0.090, 0.498, 0.115, 0.660],  # func_bkg2
    [0.001, 0.090, 0.498, 0.095,  0.530],  # func_bkg3
    [0.001, 0.090, 0.498, 0.140,  0.530]   # func_bkg4
]

# Labels for the parameter changes
labels = [
    r'$\alpha_R:0.530\rightarrow0.400$',
    r'$\alpha_R:0.530\rightarrow0.660$',
    r'$\sigma_R:0.115\rightarrow0.095$',
    r'$\sigma_R:0.115\rightarrow0.140$'
]

# x range for plotting
x = np.linspace(-1, 3, 1000)

# Normalize to total area = 1 (PDF)
def normalize_to_pdf(values, x):
    integral = trapezoid(values, x)
    if integral == 0:
        return values
    return values / integral

# Compute normalized PDF for the original background
orig_vals = expgausexp(x, *orig_bg_params)
orig_pdf  = normalize_to_pdf(orig_vals, x)

# Open a PdfPages object to save four plots in a single PDF
save_pdf_path = '/data03/tianye/pbar/syscode/pdf/background_shadow_area_all.pdf'
with PdfPages(save_pdf_path) as pdf:
    for idx, mod_params in enumerate(mod_params_list):
        mod_vals = expgausexp(x, *mod_params)
        mod_pdf  = normalize_to_pdf(mod_vals, x)

        # 差值计算
        diff = orig_pdf - mod_pdf

        # 找到交点（第一个从正变负的位置）
        crossing_idx = np.where(np.diff(np.sign(diff)))[0]
        if len(crossing_idx) == 0:
            # 无交点，跳过绘图
            continue
        right_crossing = crossing_idx[0] + 1  # 向右取一点以避免边界问题

        # 仅在交点右侧且 orig > mod 时填充
        mask =  (np.arange(len(x)) >= right_crossing)

        area = trapezoid(diff[mask], x[mask])

        # 绘图
        plt.figure(figsize=(8, 5))
        plt.plot(x, orig_pdf, 'r-', linewidth=2, label='Original Background (PDF)')
        plt.plot(x, mod_pdf,  'r--', linewidth=2, label='Modified Background (PDF)')
        plt.fill_between(x, orig_pdf, mod_pdf, where=mask,
                         color='gray', alpha=0.3, interpolate=True,
                         label='Shadow Region (right of crossing)')

        plt.title(f'Background Comparison ({labels[idx]})\nShadow Area = {abs(area):.4f}', fontsize=12)
        plt.xlabel('x', fontsize=10)
        plt.ylabel('Probability Density (PDF)', fontsize=10)
        plt.legend(fontsize=9)
        plt.grid(alpha=0.3)
        plt.xlim(-1, 3)
        y_max = max(np.max(orig_pdf), np.max(mod_pdf)) * 1.2
        plt.ylim(0, y_max)

        plt.tight_layout()
        pdf.savefig()
        plt.close()


print(f"Saved all four plots in one PDF: {save_pdf_path}")



# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import trapezoid  # 使用scipy的trapezoid替代np.trapz

# # 定义expgausexp函数（基于您提供的evaluate函数）
# def expgausexp(x, x0, sigmaL, alphaL, sigmaR, alphaR):
#     # 判断输入x是否为标量
#     is_scalar = np.isscalar(x)
#     # 转为数组以便统一处理
#     if is_scalar:
#         x = np.array([x])
#     else:
#         x = np.asarray(x)
    
#     # 初始化t和v数组
#     t = np.zeros_like(x)
#     v = np.zeros_like(x)
    
#     # 创建左右掩码
#     left_mask = (x < x0)
#     right_mask = (x >= x0)
    
#     # 计算左侧t
#     t[left_mask] = (x[left_mask] - x0) / sigmaL
#     t[right_mask] = (x[right_mask] - x0) / sigmaR
    
#     # 定义三个区域的掩码
#     mask_left_exp = (t < -alphaL)
#     mask_center = (~mask_left_exp) & (t <= alphaR)
#     mask_right_exp = (t > alphaR)
    
#     # 计算左侧指数区域
#     if np.any(mask_left_exp):
#         a_left = 0.5 * alphaL * alphaL
#         v[mask_left_exp] = np.exp(a_left + alphaL * t[mask_left_exp])
    
#     # 计算中心高斯区域
#     if np.any(mask_center):
#         v[mask_center] = np.exp(-0.5 * t[mask_center] ** 2)
    
#     # 计算右侧指数区域
#     if np.any(mask_right_exp):
#         a_right = 0.5 * alphaR * alphaR
#         v[mask_right_exp] = np.exp(a_right - alphaR * t[mask_right_exp])
    
#     # 如果是标量输入，返回标量
#     if is_scalar:
#         return v[0]
#     else:
#         return v

# # 设置参数
# signal_params = [0.869, 0.336, 1.143, 0.374, 0.940]  # [x0, sigmaL, alphaL, sigmaR, alphaR]
# background_params = [0.001, 0.090, 0.498, 0.115, 0.530]  # 原始背景参数
# modified_bg_params = [0.001, 0.090, 0.498, 0.115, 0.4]  # 修改后的背景参数

# # 创建x值范围
# x = np.linspace(-1, 3, 1000)
# dx = x[1] - x[0]  # 计算步长用于归一化

# # 计算原始函数值
# signal_vals = expgausexp(x, *signal_params)
# bg_vals_original = expgausexp(x, *background_params)
# bg_vals_modified = expgausexp(x, *modified_bg_params)

# # 归一化到200个事例的函数
# def normalize_to_events(values, x, num_events=200):
#     integral = trapezoid(values, x)  # 使用scipy.integrate.trapezoid
#     return values * num_events / integral

# # 归一化信号和背景曲线
# signal_vals_norm = normalize_to_events(signal_vals, x, 200)
# bg_vals_original_norm = normalize_to_events(bg_vals_original, x, 200)
# bg_vals_modified_norm = normalize_to_events(bg_vals_modified, x, 200)

# # 计算总曲线（信号+背景）
# total_vals_original = signal_vals_norm + bg_vals_original_norm
# total_vals_modified = signal_vals_norm + bg_vals_modified_norm

# # 创建图形
# plt.figure(figsize=(12, 8))

# # 绘制信号曲线 - 归一化到200个事例
# plt.plot(x, signal_vals_norm, 'b-', linewidth=2.5, label='Signal (200 events)')

# # 绘制背景曲线 - 归一化到200个事例
# plt.plot(x, bg_vals_original_norm, 'r-', linewidth=2, label='Background (αR=0.530, 200 events)')
# plt.plot(x, bg_vals_modified_norm, 'r--', linewidth=2.5, label='Background (αR=0.350, 200 events)')

# # 绘制总曲线 - 归一化到400个事例
# plt.plot(x, total_vals_original, 'g-', linewidth=1.5, alpha=0.7, label='Total (original αR, 400 events)')
# plt.plot(x, total_vals_modified, 'g--', linewidth=2, alpha=0.9, label='Total (modified αR, 400 events)')

# # 添加标记和标签
# plt.axvline(x=signal_params[0], color='blue', linestyle=':', alpha=0.7, label=f'Signal x0={signal_params[0]}')
# plt.axvline(x=background_params[0], color='red', linestyle=':', alpha=0.7, label=f'Background x0={background_params[0]}')

# # 添加修改影响的箭头
# arrow_x = 3.5
# # 找到最接近arrow_x的索引
# arrow_idx = np.argmin(np.abs(x - arrow_x))
# arrow_y_norm = bg_vals_modified_norm[arrow_idx]

# plt.annotate('αR reduced', 
#              xy=(arrow_x, arrow_y_norm), 
#              xytext=(arrow_x-0.5, arrow_y_norm + 0.2),
#              arrowprops=dict(arrowstyle='->', color='black'),
#              fontsize=10,
#              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1, alpha=0.8))

# # 设置标题和标签
# plt.title('Signal and Background Distributions Normalized to Event Counts', fontsize=14)
# plt.xlabel('x', fontsize=12)
# plt.ylabel('Event Density', fontsize=12)
# plt.legend(fontsize=10, loc='upper right')
# plt.grid(True, alpha=0.3)

# # 设置坐标范围
# plt.xlim(-1, 3.0)
# plt.ylim(0, max(np.max(total_vals_original), np.max(total_vals_modified)) * 1.2)

# # 添加参数说明
# plt.text(0.02, 0.95, 'Signal Parameters:', transform=plt.gca().transAxes, fontsize=10)
# plt.text(0.02, 0.90, f'x0={signal_params[0]}, σL={signal_params[1]}, αL={signal_params[2]}', transform=plt.gca().transAxes, fontsize=9)
# plt.text(0.02, 0.85, f'σR={signal_params[3]}, αR={signal_params[4]}', transform=plt.gca().transAxes, fontsize=9)

# plt.text(0.02, 0.75, 'Background Parameters:', transform=plt.gca().transAxes, fontsize=10)
# plt.text(0.02, 0.70, f'x0={background_params[0]}, σL={background_params[1]}, αL={background_params[2]}', transform=plt.gca().transAxes, fontsize=9)
# plt.text(0.02, 0.65, f'σR={background_params[3]}, αR={background_params[4]} (original)', color='red', transform=plt.gca().transAxes, fontsize=9)
# plt.text(0.02, 0.60, f'σR={modified_bg_params[3]}, αR={modified_bg_params[4]} (modified)', color='red', transform=plt.gca().transAxes, fontsize=9)

# # 添加事件计数说明
# plt.text(0.02, 0.45, 'Normalized Event Counts:', transform=plt.gca().transAxes, fontsize=10, weight='bold')
# plt.text(0.02, 0.40, '- Signal: 200 events (blue line)', transform=plt.gca().transAxes, fontsize=9)
# plt.text(0.02, 0.35, '- Background: 200 events (red lines)', transform=plt.gca().transAxes, fontsize=9)
# plt.text(0.02, 0.30, '- Total (signal + background): 400 events (green lines)', transform=plt.gca().transAxes, fontsize=9)

# # 保存图像
# plt.savefig('/data03/tianye/pbar/syscode/pdf/expgausexp_plot.pdf', bbox_inches='tight')
# print("归一化后包含总曲线的图像已保存为指定路径的PDF文件")

# # 关闭图形以释放内存
# plt.close()