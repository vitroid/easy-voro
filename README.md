# easy-voro

`easy-voro` is a simple and fast Python wrapper for [voro++](math.lbl.gov/voro++/), a powerful C++ library for Voronoi analysis.

Unlike other wrappers that use complex bindings, `easy-voro` uses a clean JSON-based interface between a C++ executable and Python, making it robust, easy to build, and highly portable.

## Directory Structure

- `easy_voro/`: Python ラッパー。ビルド後は `easy_voro/voro_volumes` がここにできる。
- `src/`: C++ ソース（voro++ 呼び出し）。
- `examples/`: GROMACS 用ユーティリティ（`gromacs_to_voronoi.py`, `gromacs_to_vvol.py`, `vvol_to_yaplot.py` は macOS/Ubuntu 共通。`animate_voronoi_volumes.py` は macOS のみ想定）。
- `Makefile`: C++ ラッパーのビルド（macOS / Linux 対応）。

## Features

- **Fast**: Leverages the full speed of `voro++`.
- **Easy to Build**: No complex Cython or SWIG dependencies. Just a simple C++ compiler and `voro++`.
- **Pythonic API**: Returns data in familiar Python dictionaries and NumPy arrays, compatible with `pyvoro`.
- **Periodic Boundary Conditions**: Fully supports 3D periodic Voronoi tessellation.
- **Visualizations**: Includes tools for 3D visualization using `pyvista`.

## Prerequisites

### macOS（Homebrew）

1. **C++ コンパイラと voro++**
   - Xcode Command Line Tools: `xcode-select --install`
   - voro++: `brew install voro++`

2. **Python 環境**
   - Python 3.10+
   - Poetry（推奨）: `brew install poetry`

### Ubuntu Linux（gromacs_to_vvol.py まで利用する場合）

1. **C++ コンパイラと voro++**

   ```bash
   sudo apt update
   sudo apt install build-essential libvoro++-dev
   ```

   - `libvoro++-dev`: ヘッダとライブラリ（開発用）。実行のみなら `libvoro++1` でも可だが、本プロジェクトでは自前ビルドのため dev が必要。

2. **Python 環境**
   ```bash
   sudo apt install python3 python3-pip
   pip install --user poetry
   # または: pipx install poetry
   ```

   - Python 3.10 以上を推奨。必要なら `deadsnakes` PPA で 3.10 を入れる。

## Quick Start

macOS でも Ubuntu でも、前提を入れたあとは同じ手順です。

### Step 1: C++ ラッパーをビルド

```bash
make
```

（Makefile が `uname` で OS を判定し、macOS なら Homebrew、Linux なら `/usr` の voro++ を参照します。）

### Step 2: Python 依存関係のインストール

```bash
poetry install
```

### Step 3: ラッパーを使う

```python
from easy_voro.voro_wrapper import compute_voronoi_cells_voro
import numpy as np

# Define positions and box size
positions = np.random.rand(100, 3) * 10.0
box = np.diag([10.0, 10.0, 10.0])

# Compute Voronoi cells
cells = compute_voronoi_cells_voro(positions, box)

for i, cell in enumerate(cells):
    print(f"Cell {i}: volume = {cell['volume']}, vertices = {len(cell['vertices'])}")
```

## GROMACS Analysis Utilities

GROMACS の `.gro` ファイルから Voronoi 解析を行うユーティリティです。

**macOS / Ubuntu 共通で利用できるもの:**

```bash
# 体積と面数を表示
poetry run python examples/gromacs_to_voronoi.py < system.gro

# .gro を .vvol に変換（可視化用に事前計算しておく）
poetry run python examples/gromacs_to_vvol.py system.gro output.vvol

# .vvol を Yaplot 形式（.yap）に変換（体積偏差を RdBu_r で彩色）
poetry run python examples/vvol_to_yaplot.py output.vvol scene.yap
```

**Yaplot 形式について**（[Yaplot](https://github.com/vitroid/Yaplot) で再生）:

- サイズ（`r`）やパレット（`@`）などの状態は **そのフレーム内だけ有効** で、空行でフレームが終わるとリセットされます。
- パレット番号は通常 0〜255 で、**0〜2 はシステム用（黒・灰・白）** のため、書き換えないほうが無難です。本ツールは 3 番以降を使用します。

**macOS のみ（可視化）:**

3D アニメーションは macOS 上でのみ実行を想定しています。

```bash
# 体積偏差の 3D アニメーション（macOS のみ）
poetry run python examples/animate_voronoi_volumes.py < system.gro
```

## License

MIT
