#!/bin/bash

# 设置资源限制
ulimit -m unlimited
ulimit -t unlimited
ulimit -v unlimited

# 定义所有参数组合（9x9=81种组合）
params=(
    10 10
    10 15
    10 20
    10 25
    10 30
    10 40
    10 50
    10 75
    10 100
    15 10
    15 15
    15 20
    15 25
    15 30
    15 40
    15 50
    15 75
    15 100
    20 10
    20 15
    20 20
    20 25
    20 30
    20 40
    20 50
    20 75
    20 100
    25 10
    25 15
    25 20
    25 25
    25 30
    25 40
    25 50
    25 75
    25 100
    30 10
    30 15
    30 20
    30 25
    30 30
    30 40
    30 50
    30 75
    30 100
    40 10
    40 15
    40 20
    40 25
    40 30
    40 40
    40 50
    40 75
    40 100
    50 10
    50 15
    50 20
    50 25
    50 30
    50 40
    50 50
    50 75
    50 100
    75 10
    75 15
    75 20
    75 25
    75 30
    75 40
    75 50
    75 75
    75 100
    100 10
    100 15
    100 20
    100 25
    100 30
    100 40
    100 50
    100 75
    100 100
)

# 程序路径
program="/data03/tianye/pbar/syscode/bin/RooCombFit_bin25_unbinned"

# 日志目录
log_dir="/data03/tianye/pbar/syscode/log"
mkdir -p "$log_dir"

# 并行任务计数器
parallel_jobs=30  # 修改为30个并行任务
active_jobs=0
job_counter=0

# 运行任务函数
run_job() {
    local param1=$1
    local param2=$2
    local log_file="$log_dir/${param1}_${param2}.log"
    
    ((job_counter++))
    
    echo "[$(date +'%T')] 启动任务 [$job_counter/${#params[@]/2}]：$param1 $param2" | tee -a "$log_file"
    
    # 运行程序（后台执行）
    "$program" "$param1" "$param2" >> "$log_file" 2>&1 &
    
    local pid=$!
    echo "任务 [$job_counter] PID: $pid" | tee -a "$log_file"
    
    # 存储PID到数组，用于错误检测
    pids[$pid]=${param1}_${param2}
    echo $pid
}

# 创建关联数组存储PID和任务信息
declare -A pids
total_tasks=$(( ${#params[@]} / 2 ))
start_time=$(date +%s)

# 开始处理所有任务
for ((i=0; i<${#params[@]}; i+=2)); do
    param1=${params[i]}
    param2=${params[i+1]}
    
    # 如果并行任务已达上限，等待一个完成
    if (( active_jobs >= parallel_jobs )); then
        # 等待任意一个后台任务完成
        wait -n
        ((active_jobs--))
        completed=$?
        
        # 检查任务退出状态
        if [ $completed -ne 0 ]; then
            failed_pid=$!
            echo "警告: 任务 ${pids[$failed_pid]} (PID: $failed_pid) 失败，退出码: $completed"
        fi
    fi
    
    # 启动新任务
    pid=$(run_job $param1 $param2)
    unset pids[$pid]  # 从pids中移除已启动的任务
    ((active_jobs++))
    
    # 显示进度
    elapsed=$(( $(date +%s) - start_time ))
    completed_percent=$(( job_counter * 100 / total_tasks ))
    echo "进度: $job_counter/$total_tasks (${completed_percent}%) | 活动任务: $active_jobs | 耗时: ${elapsed}秒"
done

# 等待所有剩余任务完成
echo "等待最后 $active_jobs 个任务完成..."
wait

end_time=$(date +%s)
total_time=$((end_time - start_time))

# 显示摘要
echo -e "\n===== 任务完成摘要 ====="
echo "总任务数: $total_tasks"
echo "成功完成: $(( total_tasks - ${#pids[@]} ))"
echo "失败任务: ${#pids[@]}"
[ ${#pids[@]} -gt 0 ] && echo "失败任务详情: ${pids[@]}"
echo "总耗时: ${total_time}秒 ($((total_time/60))分$((total_time%60))秒)"
echo "平均任务时间: $((total_time/total_tasks))秒"

# 检查失败任务
if [ ${#pids[@]} -gt 0 ]; then
    echo -e "\n===== 失败任务日志 ====="
    for pid in "${!pids[@]}"; do
        echo "任务 ${pids[$pid]} 日志: $log_dir/${pids[$pid]}.log"
    done
    exit 1
else
    echo "所有任务成功完成!"
    exit 0
fi