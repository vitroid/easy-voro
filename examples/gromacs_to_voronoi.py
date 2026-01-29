#!/usr/bin/env python

import sys
import numpy as np
from easy_voro.voro_wrapper import compute_voronoi_cells_voro


def read_gro(file):
    """
    gromacsの.groファイルを読みこむ。

    あとで出力する場合にそなえ、できるだけデータをそのままの形で保持する。
    """

    # 無限ループ
    while True:
        frame = {
            "resi_id": [],
            "residue": [],
            "atom": [],
            "atom_id": [],
            "position": [],
        }

        title = file.readline()
        # 終了判定。1文字も読めない時はファイルの終わり。
        if len(title) == 0:
            return
        n_atom = int(file.readline())
        for i in range(n_atom):
            line = file.readline()
            residue_id = int(line[0:5])
            residue = line[5:10].strip()
            atom = line[10:15].strip()
            atom_id = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            # 速度は省略

            frame["resi_id"].append(residue_id)
            frame["residue"].append(residue)
            frame["atom"].append(atom)
            frame["atom_id"].append(atom_id)
            frame["position"].append([x, y, z])

        cell = [float(x) for x in file.readline().split()]

        # numpy形式に変換しておく。
        frame["resi_id"] = np.array(frame["resi_id"])
        frame["residue"] = np.array(frame["residue"])
        frame["atom"] = np.array(frame["atom"])
        frame["atom_id"] = np.array(frame["atom_id"])
        frame["position"] = np.array(frame["position"])

        # cellは行列の形にしておく。
        if len(cell) == 3:
            # 直方体セルの場合
            cell = np.diag(cell)
        else:
            # 9パラメータで指定される場合は、順番がややこしい。
            # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            x = [cell[0], cell[5], cell[7]]
            y = [cell[3], cell[1], cell[8]]
            z = [cell[4], cell[6], cell[2]]
            cell = np.array([x, y, z])

        frame["cell"] = cell
        # returnの代わりにyieldを使うと、繰り返し(iterator)にできる。
        yield frame


def read_one_frame(file):
    """
    現在のファイル位置から1フレームだけ読み込む。
    並列処理でオフセット指定して読むときに使う。
    読み取れなかった場合（EOFなど）は None を返す。
    """
    title = file.readline()
    if len(title) == 0:
        return None
    n_atom_line = file.readline()
    if len(n_atom_line) == 0:
        return None
    n_atom = int(n_atom_line)
    frame = {
        "resi_id": [],
        "residue": [],
        "atom": [],
        "atom_id": [],
        "position": [],
    }
    for i in range(n_atom):
        line = file.readline()
        if len(line) == 0:
            return None
        residue_id = int(line[0:5])
        residue = line[5:10].strip()
        atom = line[10:15].strip()
        atom_id = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        frame["resi_id"].append(residue_id)
        frame["residue"].append(residue)
        frame["atom"].append(atom)
        frame["atom_id"].append(atom_id)
        frame["position"].append([x, y, z])
    cell_line = file.readline()
    if len(cell_line) == 0:
        return None
    cell = [float(x) for x in cell_line.split()]
    frame["resi_id"] = np.array(frame["resi_id"])
    frame["residue"] = np.array(frame["residue"])
    frame["atom"] = np.array(frame["atom"])
    frame["atom_id"] = np.array(frame["atom_id"])
    frame["position"] = np.array(frame["position"])
    if len(cell) == 3:
        cell = np.diag(cell)
    else:
        x = [cell[0], cell[5], cell[7]]
        y = [cell[3], cell[1], cell[8]]
        z = [cell[4], cell[6], cell[2]]
        cell = np.array([x, y, z])
    frame["cell"] = cell
    return frame


"""
ユーティリティ関数
"""


def is_diagonal(matrix, tol=1e-10):
    """
    3x3行列が対角行列かどうかを判定する

    Parameters
    ----------
    matrix : array_like
        判定する3x3行列
    tol : float, optional
        非対角要素が0とみなす許容誤差（デフォルト: 1e-10）

    Returns
    -------
    bool
        対角行列の場合True、そうでない場合False

    Examples
    --------
    >>> import numpy as np
    >>> from codes.utils import is_diagonal
    >>> box = np.diag([10.0, 10.0, 10.0])
    >>> is_diagonal(box)
    True
    >>> box[0, 1] = 0.1
    >>> is_diagonal(box)
    False
    """
    matrix = np.asarray(matrix)

    # 3x3行列であることを確認
    if matrix.shape != (3, 3):
        raise ValueError(f"行列は3x3である必要があります。現在の形状: {matrix.shape}")

    # 非対角要素を取得
    # 対角要素以外の要素をチェック
    mask = ~np.eye(3, dtype=bool)
    off_diagonal = matrix[mask]

    # 非対角要素がすべて許容誤差以内で0かどうかを判定
    return np.allclose(off_diagonal, 0.0, atol=tol)


def main():
    for frame in read_gro(sys.stdin):
        # シミュレーションセル（対角行列を仮定）
        box = frame["cell"]
        assert is_diagonal(box), "シミュレーションセルは直方体である必要があります"

        # 炭素原子の座標
        atoms = frame["position"][frame["atom"] == "C1", :]

        # OW原子が存在しないフレームは計算をスキップ
        if atoms.shape[0] == 0:
            print(
                "Warning: this frame has no atoms for Voronoi. Skipping.",
                file=sys.stderr,
            )
            continue

        # 直方体セルのサイズ
        Lx, Ly, Lz = box[0, 0], box[1, 1], box[2, 2]

        # voro_wrapper による周期 Voronoi 解析
        cells = compute_voronoi_cells_voro(atoms, box)

        # 各 Voronoi セルの面数と体積を出力
        for i, cell in enumerate(cells):
            nfaces = len(cell["faces"])
            volume = cell["volume"]
            print(f"Cell {i} has {nfaces} faces and a volume of {volume}")


if __name__ == "__main__":
    main()
