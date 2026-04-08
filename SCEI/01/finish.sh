#!/bin/bash

# 获取当前脚本所在的路径（B文件夹的路径）
B_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 获取A文件夹的路径，即B文件夹的父目录
A_DIR="$(dirname "$B_DIR")"

# 获取B文件夹的名字（假设是数字）
B_NAME=$(basename "$B_DIR")

# 去掉B_NAME中的前导零
B_NAME=$(echo "$B_NAME" | sed 's/^0*//')

# 检查B_NAME是否为数字
if ! [[ "$B_NAME" =~ ^[0-9]+$ ]]; then
    echo "Error: B folder name '$B_NAME' is not a valid number."
    exit 1
fi

# 计算C文件夹的名字（B文件夹名字+1）
C_NAME=$((B_NAME + 1))

# 如果C_NAME小于10，前面加一个0
if [ "$C_NAME" -lt 10 ]; then
    C_NAME="0$C_NAME"
fi

# 复制B文件夹及其内容到A文件夹中的C文件夹
cp -r "$B_DIR" "$A_DIR/$C_NAME"

# 进入C文件夹
cd "$A_DIR/$C_NAME"
