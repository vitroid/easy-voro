#!/usr/bin/env python
"""
Gromacsファイルを読み込み、Voronoi計算結果を.vvolファイルに保存するツール

各フレームのVoronoiセルの中心座標と体積をCSV形式で保存します。
これにより、アニメーション表示時に毎回Voronoi計算を行う必要がなくなります。
"""

import sys
import numpy as np
from pathlib import Path
from tqdm import tqdm
from easy_voro.voro_wrapper import compute_voronoi_cells_voro

# 相対インポートまたは絶対インポートに対応
try:
    from .gromacs_to_voronoi import read_gro, read_one_frame, is_diagonal
except ImportError:
    from gromacs_to_voronoi import read_gro, read_one_frame, is_diagonal


def format_number(x):
    """数値を小数点以下4桁の文字列に変換"""
    return f"{float(x):.4f}"


def process_frame(frame, atom_type="C1"):
    """
    1フレームを処理してVoronoi計算を実行
    
    Parameters
    ----------
    frame : dict
        read_gro()で読み込んだフレームデータ
    atom_type : str
        処理する原子タイプ（デフォルト: "C1"）
    
    Returns
    -------
    dict or None
        処理結果の辞書、またはエラー時はNone
    """
    box = frame["cell"]
    
    # 対角行列でない場合はエラー
    if not is_diagonal(box):
        return None
    
    # 指定された原子タイプのマスク
    atom_mask = frame["atom"] == atom_type
    if not np.any(atom_mask):
        return None
    
    positions = frame["position"][atom_mask, :]
    Lx, Ly, Lz = box[0, 0], box[1, 1], box[2, 2]
    
    # 周期境界条件を適用
    if Lx > 0:
        positions %= [Lx, Ly, Lz]
    
    try:
        # Voronoi計算を実行
        cells = compute_voronoi_cells_voro(positions, box)
        
        # セルの中心座標と体積を取得（小数点以下4桁）
        centers = [c["original"] for c in cells]
        volumes = [c["volume"] for c in cells]
        
        return {
            "box": box,
            "centers": centers,
            "volumes": volumes,
        }
    except Exception as e:
        print(f"Error processing frame: {e}", file=sys.stderr)
        return None


def convert_gromacs_to_vvol(input_file, output_file, atom_type="C1"):
    """
    Gromacsファイルを.vvolファイルに変換
    
    Parameters
    ----------
    input_file : str or Path
        入力Gromacsファイルのパス
    output_file : str or Path
        出力.vvolファイルのパス
    atom_type : str
        処理する原子タイプ（デフォルト: "C1"）
    """
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    if not input_path.exists():
        raise FileNotFoundError(f"入力ファイルが見つかりません: {input_file}")
    
    # ファイルを開いて処理（CSV形式で保存）
    with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
        for frame_idx, frame in enumerate(tqdm(read_gro(f_in), desc="処理中", unit="frame")):
            result = process_frame(frame, atom_type)
            
            if result is None:
                # エラーの場合は空のフレームを書き込む
                f_out.write(f"# FRAME {frame_idx}\n")
                f_out.write("BOX 0.0 0.0 0.0\n")
                f_out.write("\n")
                continue
            
            box = result["box"]
            centers = result["centers"]
            volumes = result["volumes"]
            n_cells = len(centers)
            
            # フレームヘッダー
            f_out.write(f"# FRAME {frame_idx} {n_cells}\n")
            
            # ボックス情報（対角成分のみ、小数点以下4桁）
            Lx, Ly, Lz = box[0, 0], box[1, 1], box[2, 2]
            f_out.write(f"BOX {format_number(Lx)} {format_number(Ly)} {format_number(Lz)}\n")
            
            # 各セルの中心座標と体積（x y z volume の形式、小数点以下4桁）
            for center, volume in zip(centers, volumes):
                x, y, z = center
                f_out.write(f"{format_number(x)} {format_number(y)} {format_number(z)} {format_number(volume)}\n")
            
            # フレーム区切り（空行）
            f_out.write("\n")
    
    print(f"完了: {output_path}", file=sys.stderr)


def _index_gro_offsets(input_path):
    """groファイルを1パス走査し、各フレームの先頭バイトオフセットのリストを返す"""
    offsets = []
    with open(input_path, "rb") as f:
        while True:
            off = f.tell()
            if not f.readline():
                break
            line2 = f.readline()
            if not line2:
                break
            try:
                n_atoms = int(line2.strip())
            except ValueError:
                break
            for _ in range(n_atoms + 1):
                if not f.readline():
                    break
            offsets.append(off)
    return offsets


def _worker_one_frame(args):
    """
    1フレームを処理するワーカー（ProcessPoolExecutor用・トップレベルで定義）。
    args: (input_path, offset, frame_idx, atom_type)
    returns: (frame_idx, result_dict or None)
    """
    import os
    import sys
    input_path, offset, frame_idx, atom_type = args
    # 子プロセスで examples から import できるようにする
    _examples_dir = Path(__file__).resolve().parent
    if str(_examples_dir) not in sys.path:
        sys.path.insert(0, str(_examples_dir))
    from easy_voro.voro_wrapper import compute_voronoi_cells_voro
    from gromacs_to_voronoi import read_one_frame, is_diagonal
    import numpy as np

    with open(input_path, "r") as f:
        f.seek(offset)
        frame = read_one_frame(f)
    if frame is None:
        return (frame_idx, None)
    box = frame["cell"]
    if not is_diagonal(box):
        return (frame_idx, None)
    atom_mask = frame["atom"] == atom_type
    if not np.any(atom_mask):
        return (frame_idx, None)
    positions = np.asarray(frame["position"], dtype=np.float64)[atom_mask, :]
    Lx, Ly, Lz = float(box[0, 0]), float(box[1, 1]), float(box[2, 2])
    if Lx > 0:
        positions = positions % np.array([Lx, Ly, Lz])
    try:
        cells = compute_voronoi_cells_voro(positions, box)
    except Exception:
        return (frame_idx, None)
    centers = [[round(float(x), 4) for x in c["original"]] for c in cells]
    volumes = [round(float(c["volume"]), 4) for c in cells]
    box_list = [[round(float(box[i, j]), 4) for j in range(3)] for i in range(3)]
    return (frame_idx, {"box": box_list, "centers": centers, "volumes": volumes})


def convert_gromacs_to_vvol_parallel(input_file, output_file, atom_type="C1", workers=None):
    """
    groをvvolに変換（並列版）。manycore環境で高速化する。
    直列版 convert_gromacs_to_vvol は変更しない。
    """
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    input_path = Path(input_file)
    output_path = Path(output_file)
    if not input_path.exists():
        raise FileNotFoundError(f"入力ファイルが見つかりません: {input_file}")

    n_workers = workers if workers is not None else (os.cpu_count() or 4)
    n_workers = max(1, int(n_workers))

    print(f"オフセットをスキャン中: {input_path}", file=sys.stderr)
    offsets = _index_gro_offsets(input_path)
    n_frames = len(offsets)
    print(f"フレーム数: {n_frames}, ワーカー数: {n_workers}", file=sys.stderr)

    input_path_str = str(input_path.resolve())
    task_args = [(input_path_str, off, i, atom_type) for i, off in enumerate(offsets)]

    # 完了順はバラバラでも、書き出しは next_write の昇順のみ行うのでフレーム順序は保証される
    pending = {}
    next_write = 0

    with open(output_path, "w") as f_out:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_worker_one_frame, a): a[2] for a in task_args}
            for fut in tqdm(as_completed(futures), total=n_frames, desc="処理中", unit="frame"):
                frame_idx = futures[fut]
                try:
                    idx, result = fut.result()
                except Exception as e:
                    print(f"Frame {frame_idx} エラー: {e}", file=sys.stderr)
                    result = None
                pending[idx] = result
                while next_write in pending:
                    res = pending.pop(next_write)
                    if res is None:
                        f_out.write(f"# FRAME {next_write}\n")
                        f_out.write("BOX 0.0 0.0 0.0\n\n")
                    else:
                        box = res["box"]
                        centers = res["centers"]
                        volumes = res["volumes"]
                        n_cells = len(centers)
                        f_out.write(f"# FRAME {next_write} {n_cells}\n")
                        Lx, Ly, Lz = box[0][0], box[1][1], box[2][2]
                        f_out.write(f"BOX {format_number(Lx)} {format_number(Ly)} {format_number(Lz)}\n")
                        for (x, y, z), v in zip(centers, volumes):
                            f_out.write(f"{format_number(x)} {format_number(y)} {format_number(z)} {format_number(v)}\n")
                        f_out.write("\n")
                    next_write += 1
    print(f"完了: {output_path}", file=sys.stderr)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Gromacsファイルを.vvolファイルに変換",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
例:
  %(prog)s input.gro output.vvol
  %(prog)s input.gro output.vvol --atom-type OW
        """
    )
    parser.add_argument("input_file", help="入力Gromacsファイル (.gro)")
    parser.add_argument("output_file", help="出力.vvolファイル")
    parser.add_argument(
        "--atom-type",
        default="C1",
        help="処理する原子タイプ (デフォルト: C1)"
    )
    parser.add_argument(
        "--workers",
        "-j",
        type=int,
        default=1,
        metavar="N",
        help="並列ワーカー数 (1=直列, 2以上で並列。省略時は1)"
    )
    
    args = parser.parse_args()
    
    try:
        if args.workers <= 1:
            convert_gromacs_to_vvol(args.input_file, args.output_file, args.atom_type)
        else:
            convert_gromacs_to_vvol_parallel(
                args.input_file, args.output_file, args.atom_type, workers=args.workers
            )
    except Exception as e:
        print(f"エラー: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
