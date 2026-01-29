#!/usr/bin/env python
"""
.vvol を Yaplot 形式（.yap）に変換する。

体積偏差（平均からのずれ）を animate_voronoi_volumes.py と同様に RdBu_r で彩色し、
各サイトを円（Circle）で描画します。

yaplotlib: https://vitroid.github.io/yaplotlib/yaplotlib.html
Yaplot の仕様: https://github.com/vitroid/Yaplot
- サイズ（r）やパレット（@）はそのフレーム内だけ有効で、空行でフレームが終わるとリセットされる。
- パレット番号 0〜2 はシステム用（黒・灰・白）のため、本スクリプトでは 3 番以降を使用する。
"""

import sys
import numpy as np
from pathlib import Path

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.cm as cm
except ImportError:
    cm = None

try:
    from yaplotlib import YaplotDocument
except ImportError:
    YaplotDocument = None

try:
    from .vvol_indexer import VvolIndexer
except ImportError:
    from vvol_indexer import VvolIndexer


# パレットは 0〜2 を避け、3 番から使用（最大 255 まで。余裕をみて 3〜66 の 64 段階）
PALETTE_START = 3
NCOLORS = 64


def build_rdbu_r_palette(ncolors, cmap_name="RdBu_r"):
    """RdBu_r の離散パレットを [ (r,g,b), ... ] で返す。インデックス 0 が青側、ncolors-1 が赤側。"""
    if cm is None:
        raise RuntimeError("matplotlib が必要です: pip install matplotlib")
    cmap = cm.get_cmap(cmap_name)
    palette = []
    for i in range(ncolors):
        rgba = cmap(i / max(ncolors - 1, 1))
        r, g, b = (
            int(round(rgba[0] * 255)),
            int(round(rgba[1] * 255)),
            int(round(rgba[2] * 255)),
        )
        palette.append((r, g, b))
    return palette


def dev_to_palette_index(dev, vmin, vmax, ncolors):
    """dev を [vmin, vmax] で正規化し、パレット番号 0..ncolors-1 にマップ。"""
    if vmax <= vmin:
        return (ncolors - 1) // 2
    norm = (dev - vmin) / (vmax - vmin)
    norm = max(0.0, min(1.0, norm))
    idx = int(norm * (ncolors - 1) + 0.5)
    return min(idx, ncolors - 1)


def add_frame_to_doc(doc, data, radius, palette_rgb, palette_start, ncolors):
    """
    1 フレーム分のデータを doc の現在フレームに追加する。
    彩色は animate_voronoi_volumes.py と同様: 平均体積に対する±5%の範囲でグラデーション。
    """
    pos = data["pos"]
    dev = data["dev"]
    mean_v = data["mean_v"]
    dev_pct = dev / mean_v * 100.0
    vmin, vmax = -5.0, 5.0

    frame = doc.current
    frame.Size(radius)
    for i in range(ncolors):
        r, g, b = palette_rgb[i]
        frame.SetPalette(palette_start + i, r, g, b, maxval=255.0)

    for j in range(len(pos)):
        idx = dev_to_palette_index(dev_pct[j], vmin, vmax, ncolors)
        frame.Color(palette_start + idx)
        frame.Circle(pos[j].tolist())


def vvol_to_yaplot(
    vvol_path, output_path, radius=0.07, ncolors=NCOLORS, palette_start=PALETTE_START
):
    """
    .vvol を Yaplot 形式で output_path に書き出す。

    Parameters
    ----------
    vvol_path : str or Path
        入力 .vvol ファイル
    output_path : str or Path
        出力 .yap ファイル（または - で stdout）
    radius : float
        円マークの半径（Yaplot の r コマンド）
    ncolors : int
        体積偏差の段階数（パレット 3 番から ncolors 個使用）
    palette_start : int
        使用開始パレット番号（0〜2 は避けるためデフォルト 3）
    """
    if YaplotDocument is None:
        raise RuntimeError("yaplotlib が必要です: pip install yaplotlib")

    if palette_start < 3 or palette_start + ncolors > 256:
        raise ValueError(
            "パレットは 3〜255 の範囲で指定してください（0〜2 はシステム予約）"
        )

    palette_rgb = build_rdbu_r_palette(ncolors)
    indexer = VvolIndexer(vvol_path)
    nframes = len(indexer)

    if tqdm is not None:
        frame_iter = tqdm(range(nframes), desc="vvol→yaplot", unit="frame")
    else:
        frame_iter = range(nframes)

    doc = YaplotDocument()
    first_frame = True
    written_count = 0
    for frame_idx in frame_iter:
        data = indexer.get_frame_data(frame_idx)
        if data is None:
            continue
        if not first_frame:
            doc.new_frame()
        first_frame = False
        add_frame_to_doc(doc, data, radius, palette_rgb, palette_start, ncolors)
        written_count += 1

    if str(output_path) == "-":
        sys.stdout.write(doc.dumps())
    else:
        doc.save(output_path)
        if written_count:
            print(f"Wrote {written_count} frames to {output_path}", file=sys.stderr)


def main():
    import argparse

    p = argparse.ArgumentParser(
        description=".vvol を Yaplot 形式（.yap）に変換（体積偏差を RdBu_r で彩色）"
    )
    p.add_argument("input_file", help="入力 .vvol ファイル")
    p.add_argument(
        "output_file",
        nargs="?",
        default="-",
        help="出力 .yap ファイル（省略時は標準出力）",
    )
    p.add_argument(
        "-r",
        "--radius",
        type=float,
        default=0.07,
        help="円の半径 (default: 0.07)",
    )
    p.add_argument(
        "-n",
        "--ncolors",
        type=int,
        default=NCOLORS,
        help=f"パレット段階数 (default: {NCOLORS})",
    )
    args = p.parse_args()

    vvol_to_yaplot(
        args.input_file,
        args.output_file,
        radius=args.radius,
        ncolors=args.ncolors,
    )


if __name__ == "__main__":
    main()
