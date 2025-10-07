# infer_pcfg_julia

## 概要 / Overview

このリポジトリは、ノンパラメトリックベイズ手法を用いた(k,l)-文脈依存確率文法の学習アルゴリズムの実装を提供します。

This repository provides an implementation of learning algorithms for (k,l)-context-sensitive probabilistic grammars using nonparametric Bayesian approaches.

## 論文 / Reference Paper

このプログラムは以下の論文で使用されました：

This program was used in the following research paper:

Chihiro Shibata:  
**Learning (k,l)-context-sensitive probabilistic grammars with nonparametric Bayesian approach.**  
Machine Learning, 113(5): 3267-3301 (2024)

## 説明 / Description

### 日本語

(k,l)-文脈依存確率文法は、形式言語理論における文法の拡張で、文脈を考慮した確率的な生成規則を持ちます。本実装では、ノンパラメトリックベイズ手法（特にディリクレ過程を用いた手法）を使用して、データから自動的にこれらの文法を学習します。

#### 主な特徴
- **ノンパラメトリックベイズアプローチ**: モデルの複雑さをデータから自動的に決定
- **(k,l)-文脈依存性**: 左側k個、右側l個のシンボルを文脈として考慮
- **確率文法の推定**: 文字列データから確率的生成規則を学習

### English

(k,l)-context-sensitive probabilistic grammars are extensions of formal grammars that incorporate probabilistic generation rules with context awareness. This implementation uses nonparametric Bayesian methods (specifically Dirichlet process-based approaches) to automatically learn these grammars from data.

#### Key Features
- **Nonparametric Bayesian Approach**: Automatically determines model complexity from data
- **(k,l)-Context Sensitivity**: Considers k symbols on the left and l symbols on the right as context
- **Probabilistic Grammar Inference**: Learns probabilistic generation rules from string data

## 使用方法 / Usage

### 必要要件 / Requirements
- Julia (推奨バージョン: 1.0以降 / Recommended version: 1.0 or later)

### インストール / Installation

```bash
git clone https://github.com/ChihiroShibata/infer_pcfg_julia.git
cd infer_pcfg_julia
```

### 実行例 / Example Usage

Juliaコードが追加された後、詳細な使用方法をここに記載します。

Detailed usage instructions will be provided when Julia code is added to the repository.

## ライセンス / License

MIT License - 詳細は [LICENSE](LICENSE) ファイルを参照してください。

MIT License - See [LICENSE](LICENSE) file for details.

## 著者 / Author

C. S.

## 引用 / Citation

この実装を研究で使用する場合は、以下の論文を引用してください：

If you use this implementation in your research, please cite:

```bibtex
@article{Shibata2024,
  author    = {Chihiro Shibata},
  title     = {Learning (k,l)-context-sensitive probabilistic grammars with nonparametric Bayesian approach},
  journal   = {Machine Learning},
  volume    = {113},
  number    = {5},
  pages     = {3267--3301},
  year      = {2024}
}
```