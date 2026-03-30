#!/usr/bin/env bash
# pull_results.sh — 从远程拉取模拟结果到对应 results_seed{SEED}/ 目录
#
# 用法: bash pull_results.sh
# 依赖: ssh/rsync，~/.ssh/config 中已配置 cufe.zbw.uk 的 ProxyCommand

set -euo pipefail

REMOTE="u2023312303@cufe.zbw.uk"
REMOTE_MAIN="~/simulation/main.R"
REMOTE_RESULTS="~/simulation/results/"
LOCAL_BASE="$(cd "$(dirname "$0")" && pwd)"  # simulation/ 目录

# 1. 检测远程 seed
echo ">>> 检测远程 seed..."
SEED=$(ssh "$REMOTE" "grep -oP '(?<=set\.seed\()[0-9]+' $REMOTE_MAIN")
if [[ -z "$SEED" ]]; then
  echo "错误: 无法从远程 main.R 中读取 seed 值" >&2
  exit 1
fi
echo ">>> 远程 seed = $SEED"

# 2. 确定本地目标目录
LOCAL_DIR="$LOCAL_BASE/results_seed${SEED}"
echo ">>> 本地目标目录: $LOCAL_DIR"
mkdir -p "$LOCAL_DIR"

# 3. 备份（如已存在同名目录且非空）
if [[ -n "$(ls -A "$LOCAL_DIR" 2>/dev/null)" ]]; then
  BACKUP_DIR="${LOCAL_DIR}_backup_$(date +%Y%m%d_%H%M%S)"
  echo ">>> 备份已有数据到 $BACKUP_DIR ..."
  cp -r "$LOCAL_DIR" "$BACKUP_DIR"
fi

# 4. 拉取（先清空目标目录，再 scp）
echo ">>> 清空目标目录..."
rm -rf "${LOCAL_DIR:?}"/*

echo ">>> 开始拉取（scp）..."
scp -r "$REMOTE:$REMOTE_RESULTS*" "$LOCAL_DIR/"

echo ""
echo ">>> 完成！结果已同步到: $LOCAL_DIR"
echo ">>> 文件列表:"
ls -lh "$LOCAL_DIR"
