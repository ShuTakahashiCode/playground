FROM rust:1.86.0-bullseye

# 開発に便利なツールをインストール（任意）
RUN apt-get update && apt-get install -y \
    pkg-config \
    libssl-dev \
    build-essential \
    && apt-get clean

# Rustのツールチェインをインストール
RUN rustup target add x86_64-pc-windows-gnu 
RUN apt install -y gcc-mingw-w64

# 一式をコピー
COPY . .
