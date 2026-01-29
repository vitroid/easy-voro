#!/usr/bin/env python
"""
Voronoiセルの体積偏差を色で表した3Dアニメーションを作成
"""

import sys
import os
import warnings
import numpy as np
import pyvista as pv
from pathlib import Path
from tqdm import tqdm
from easy_voro.voro_wrapper import compute_voronoi_cells_voro

# 相対インポートまたは絶対インポートに対応
try:
    from .gromacs_to_voronoi import read_gro, is_diagonal
except ImportError:
    from gromacs_to_voronoi import read_gro, is_diagonal

try:
    from .vvol_indexer import VvolIndexer
except ImportError:
    from vvol_indexer import VvolIndexer

warnings.filterwarnings("ignore", category=UserWarning)


class GroIndexer:
    """Gromacsファイルからフレームを読み込むクラス（従来の実装）"""

    def __init__(self, filename):
        self.filename = filename
        self.offsets = []
        self._scan_frames()

    def _scan_frames(self):
        file_size = Path(self.filename).stat().st_size
        with open(self.filename, "rb") as f, tqdm(
            total=file_size, unit="B", unit_scale=True, desc="Indexing", file=sys.stderr
        ) as pbar:
            while True:
                offset = f.tell()
                line = f.readline()
                if not line:
                    break
                try:
                    line2 = f.readline()
                    if not line2:
                        break
                    n_atoms = int(line2.strip())
                    for _ in range(n_atoms + 1):
                        f.readline()
                    self.offsets.append(offset)
                    pbar.n = f.tell()
                    pbar.refresh()
                except:
                    break

    def get_frame(self, frame_idx):
        if frame_idx < 0 or frame_idx >= len(self.offsets):
            return None
        with open(self.filename, "rb") as f:
            f.seek(self.offsets[frame_idx])

            class Decoder:
                def __init__(self, f):
                    self.f = f

                def readline(self):
                    l = self.f.readline()
                    return l.decode("utf-8") if l else ""

            try:
                return next(read_gro(Decoder(f)))
            except:
                return None

    def __len__(self):
        return len(self.offsets)


def create_animation_interactive(indexer, atom_type="C1", use_vvol=False):
    if len(indexer) == 0:
        return

    plotter = pv.Plotter()
    cache = {}
    state = {"idx": 0, "vmin": -0.01, "vmax": 0.01}

    def get_processed_frame(idx):
        if idx in cache:
            return cache[idx]

        # `.vvol`ファイルを使用する場合
        if use_vvol and isinstance(indexer, VvolIndexer):
            data = indexer.get_frame_data(idx)
            if data is not None:
                cache[idx] = data
            return data

        # 従来の方法（Gromacsファイルから直接計算）
        frame = indexer.get_frame(idx)
        if frame is None:
            return None
        box = frame["cell"]
        atom_mask = frame["atom"] == atom_type
        if not np.any(atom_mask):
            return None
        positions = frame["position"][atom_mask, :]
        Lx, Ly, Lz = box[0, 0], box[1, 1], box[2, 2]
        if Lx > 0:
            positions %= [Lx, Ly, Lz]
        try:
            cells = compute_voronoi_cells_voro(positions, box)
            volumes = np.array([c["volume"] for c in cells])
            devs = volumes - np.mean(volumes)
            data = {
                "pos": positions,
                "dev": devs,
                "mean_v": np.mean(volumes),
                "box": box,
            }
            cache[idx] = data
            return data
        except:
            return None

    def update_scale(idx):
        d = get_processed_frame(idx)
        if d is not None and d["dev"].size > 0:
            # 平均体積に対する±5%の範囲でグラデーション
            state["vmin"], state["vmax"] = -5.0, 5.0

    def show_frame(idx):
        idx = int(max(0, min(idx, len(indexer) - 1)))
        state["idx"] = idx

        # デバッグ情報の表示
        if isinstance(indexer, VvolIndexer):
            print(f"DEBUG: Frame {idx+1}/{len(indexer)}", file=sys.stderr)
        else:
            with open(indexer.filename, "rb") as f:
                f.seek(indexer.offsets[idx])
                title = f.readline().decode("utf-8").strip()
                print(f"DEBUG: Frame {idx+1}: {title}", file=sys.stderr)

        data = get_processed_frame(idx)
        plotter.clear()

        if data is None:
            plotter.add_text(f"Frame {idx+1}: Error", color="red")
        else:
            pc = pv.PolyData(data["pos"])
            # 平均体積に対する%偏差に変換
            dev_pct = data["dev"] / data["mean_v"] * 100.0
            pc["d"] = dev_pct
            glyph = pc.glyph(geom=pv.Sphere(radius=0.07), orient=False, scale=False)
            plotter.add_mesh(
                glyph, scalars="d", cmap="RdBu_r", clim=[state["vmin"], state["vmax"]]
            )
            box = data["box"]
            plotter.add_mesh(
                pv.Box(bounds=[0, box[0, 0], 0, box[1, 1], 0, box[2, 2]]),
                style="wireframe",
                opacity=0.3,
            )

        info = (
            f"Frame {idx+1}/{len(indexer)}\n"
            "Arrows: 1, Shift: 10, Opt+Shift: 100, Cmd: 1000\n"
            "[s]: Autoscale, [q]: Quit"
        )
        plotter.add_text(info, font_size=10, name="info")

    def on_key_press(obj, event):
        key = obj.GetKeySym()
        shift = obj.GetShiftKey()
        alt = obj.GetAltKey()
        ctrl = obj.GetControlKey()  # macOS Command or Control

        step = 0
        if key == "Right":
            step = 1
        elif key == "Left":
            step = -1

        if step != 0:
            if ctrl:  # Command or Ctrl
                step *= 1000
            elif shift and alt:
                step *= 100
            elif shift:
                step *= 10

            show_frame(state["idx"] + step)
            obj.SetAbortFlag(1)

        elif key == "s":
            update_scale(state["idx"])
            show_frame(state["idx"])
            obj.SetAbortFlag(1)
        elif key == "q":
            plotter.close()
            obj.SetAbortFlag(1)

    plotter.iren.add_observer("KeyPressEvent", on_key_press, 1.0)

    update_scale(0)
    show_frame(0)
    plotter.show()


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("input_file", help="入力ファイル (.gro または .vvol)")
    p.add_argument("--interactive", action="store_true")
    p.add_argument(
        "--atom-type", default="C1", help="原子タイプ（.groファイル使用時のみ）"
    )
    args = p.parse_args()

    if args.interactive:
        input_path = Path(args.input_file)

        # `.vvol`ファイルの場合はVvolIndexerを使用（prescan + seek）
        if input_path.suffix == ".vvol":
            indexer = VvolIndexer(args.input_file)
            create_animation_interactive(indexer, args.atom_type, use_vvol=True)
        else:
            # Gromacsファイルの場合は従来の方法
            indexer = GroIndexer(args.input_file)
            create_animation_interactive(indexer, args.atom_type, use_vvol=False)


if __name__ == "__main__":
    main()
