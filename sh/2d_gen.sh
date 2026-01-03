#!/bin/bash

# 清空文件
> /data03/tianye/pbar/syscode/sh/2d_gen.txt

# 定义数组，也可以直接写在循环中
pr_values="10 15 20 25 30 40 50 75 100 200 400"
pi_values="10 15 20 25 30 40 50 75 100 200 400"

# 或者直接循环
for pr in $pr_values; do
    for pi in $pi_values; do
        echo "$pr $pi" >> /data03/tianye/pbar/syscode/sh/2d_gen.txt
    done
done