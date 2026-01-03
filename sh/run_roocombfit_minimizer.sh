#!/bin/bash

# 设置资源限制
ulimit -m unlimited      # 设置内存限制为无限制
#ulimit -n 100000       # 设置最大打开文件数为10000
ulimit -t unlimited      # 设置CPU时间为无限制
ulimit -v unlimited      # 设置虚拟内存为无限制

# 运行你的程序
/data03/tianye/pbar/macro/unbinned_fit/bin1/macro/fixneg/macro/Nmc100/bin/RooCombFit_bin25_minimizer "$@"

