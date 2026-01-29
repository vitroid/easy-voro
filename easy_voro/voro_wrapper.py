import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


def _find_voro_volumes_exe():
    """voro_volumes 実行ファイルのパスを探す"""
    # 1. 環境変数
    exe_env = os.environ.get("VORO_VOLUMES_EXE") or os.environ.get("EASY_VORO_EXE")
    if exe_env and os.path.isfile(exe_env):
        return exe_env
    # 2. このモジュールと同じディレクトリ（インストール先 or ソースツリー）
    pkg_dir = Path(__file__).resolve().parent
    exe_in_pkg = pkg_dir / "voro_volumes"
    if exe_in_pkg.is_file():
        return str(exe_in_pkg)
    # 3. カレントディレクトリからの easy_voro/voro_volumes（poetry run でプロジェクトルートから実行時）
    for cwd in (Path.cwd(), Path.cwd().parent):
        cand = cwd / "easy_voro" / "voro_volumes"
        if cand.is_file():
            return str(cand)
    # 4. PATH
    exe_in_path = shutil.which("voro_volumes")
    if exe_in_path:
        return exe_in_path
    raise FileNotFoundError(
        "voro_volumes 実行ファイルが見つかりません。"
        " プロジェクトルートで `make` を実行してビルドするか、"
        "環境変数 VORO_VOLUMES_EXE にバイナリのパスを設定してください。"
    )


def compute_voronoi_cells_voro(positions: np.ndarray, box: np.ndarray, exe: str = None):
    """
    Voro++ の C++ ラッパー (`voro_volumes`) を叩いて、pyvoro 互換の cell 構造のリストを返す。

    Parameters
    ----------
    positions : (N, 3) ndarray
        カルテシアン座標の点群（[0, Lx]×[0, Ly]×[0, Lz] 内にあることを前提）。
    box : (3, 3) ndarray
        直方体セルの行列（対角成分に Lx, Ly, Lz を持つことを前提）。
    exe : str, optional
        Voro++ ラッパーバイナリへのパス。
        指定されない場合は、このスクリプトと同じディレクトリにある `voro_volumes` を使用します。

    Returns
    -------
    cells : list[dict]
        各セルについて、
        {
            "original": [x, y, z],
            "volume": float,
            "vertices": [[x, y, z], ...],
            "faces": [
                {"adjacent_cell": int, "vertices": [int, ...]},
                ...
            ],
            "adjacency": [[int, ...], ...],  # 各頂点ごとの隣接頂点インデックス
        }
        のような辞書を要素とするリスト。
    """
    if exe is None:
        exe = _find_voro_volumes_exe()
    positions = np.asarray(positions, float)
    Lx, Ly, Lz = float(box[0, 0]), float(box[1, 1]), float(box[2, 2])

    # 一時ファイルに書き出し
    # delete=False にして、Windows等での競合を避ける（Pathlibで管理）
    fd, path = tempfile.mkstemp(suffix=".txt", prefix="voro_pos_")
    pos_path = Path(path)
    try:
        with open(pos_path, "w") as f:
            np.savetxt(f, positions)

        # ファイルを閉じてから別プロセスで開く
        os.close(fd)

        cmd = [exe, str(Lx), str(Ly), str(Lz)]
        with open(pos_path, "r") as f_in:
            proc = subprocess.run(
                cmd,
                stdin=f_in,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
    finally:
        pos_path.unlink(missing_ok=True)

    # JSON Lines をパースしてセル構造を組み立てる
    cells = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        cid = rec["id"]

        cell = {
            "original": rec["position"],
            "volume": rec["volume"],
            "vertices": rec.get("vertices", []),
            "faces": rec.get("faces", []),
        }

        # adjacency（頂点グラフ）は faces から構成する
        verts = cell["vertices"]
        n_vertices = len(verts)
        adjacency = [set() for _ in range(n_vertices)]

        for face in cell["faces"]:
            vs = face.get("vertices", [])
            m = len(vs)
            if m < 2:
                continue
            for k in range(m):
                a = vs[k]
                b = vs[(k + 1) % m]
                if 0 <= a < n_vertices and 0 <= b < n_vertices:
                    adjacency[a].add(b)
                    adjacency[b].add(a)

        cell["adjacency"] = [sorted(list(nei)) for nei in adjacency]
        cells[cid] = cell

    # ID順に並べたリストにする（pyvoro と似た挙動にしておく）
    ordered_cells = [cells[i] for i in sorted(cells.keys())]
    return ordered_cells
