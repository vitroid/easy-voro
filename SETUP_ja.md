# easy-voro (Voronoi 解析ツール)

このディレクトリは、`voro++` を Python から簡単に利用するためのラッパー `easy-voro` と、それを利用した解析ツール群を含んでいます。
GROMACS の構造ファイル（`.gro`）の読み込みにも対応しています。

## 準備するもの (Mac + Homebrew)

Mac環境で Homebrew を使用して必要なツールをインストールします。

### 1. C++ コンパイラとライブラリ

- **Xcode Command Line Tools**: `g++` などの開発ツール一式が必要です。インストールしていない場合は、以下のコマンドを実行してください。
  ```bash
  xcode-select --install
  ```
- **voro++**: Voronoi解析のコアライブラリです。
  ```bash
  brew install voro++
  ```

### 2. Python 環境

- **Python 3.10**: このプロジェクトは Python 3.10 系を想定しています。
- **Poetry**: Pythonのパッケージ管理ツールです。
  ```bash
  brew install poetry
  ```

---

## セットアップ手順

### ステップ1: C++ ラッパーのビルド

Voronoi解析を高速に行うための C++ プログラム `voro_volumes` をコンパイルします。

```bash
make
```

コンパイルされたバイナリは `src/easy_voro/voro_volumes` に生成されます。

### ステップ2: Python 依存関係のインストール

Poetryを使用して、必要なPythonライブラリ（numpy, pyvista, matplotlibなど）をインストールします。

```bash
poetry install
```

---

## 実行方法

### 1. Voronoi解析の結果（体積・面数）を表示する

`examples/gromacs_to_voronoi.py` は、標準入力から `.gro` ファイルを受け取り、解析結果を標準出力に表示します。

```bash
poetry run python examples/gromacs_to_voronoi.py < your_structure.gro
```

### 2. Voronoiセルの体積偏差を可視化・アニメーション化する

`examples/animate_voronoi_volumes.py` を使用して、3D可視化や動画作成が可能です。

```bash
# 動画（MP4）として保存する場合
poetry run python examples/animate_voronoi_volumes.py < your_structure.gro

# インタラクティブなウィンドウで確認する場合
poetry run python examples/animate_voronoi_volumes.py --interactive < your_structure.gro
```

---

## 注意事項

- **セルの形状**: 現在のスクリプトは直方体（Orthorhombic）セルのみに対応しています。
- **入力フォーマット**: GROMACSの固定形式 `.gro` ファイルを想定しています。
