"""
.vvol ファイル用のインデクサ。フレーム先頭オフセットを保持し、get_frame_data で1フレーム分のデータを返す。

vvol_to_yaplot や animate_voronoi_volumes から利用可能。pyvista に依存しない。
"""

import sys
import numpy as np
from pathlib import Path

# バイナリで "# FRAME" を探す用（デコードしないので高速）
_MARKER = b"# FRAME"
_NEWLINE_MARKER = b"\n# FRAME"


class VvolIndexer:
    """`.vvol`からVoronoi結果を読む。インデックスは初回アクセス時にバイナリスキャンで一括取得。"""

    def __init__(self, filename):
        self.filename = filename
        self._offsets = None  # 初回アクセス時にのみ構築

    def _ensure_index(self):
        if self._offsets is not None:
            return
        self._offsets = []
        path = Path(self.filename)
        file_size = path.stat().st_size
        buf_size = 4 * 1024 * 1024  # 4MB
        overlap = 32  # 境界で "# FRAME" が分割されないように
        step = buf_size - overlap
        with open(self.filename, "rb") as f:
            chunk_start = 0
            while chunk_start < file_size:
                f.seek(chunk_start)
                buf = f.read(buf_size)
                if not buf:
                    break
                # 先頭フレーム: ファイル先頭が "# FRAME"
                if (
                    chunk_start == 0
                    and buf.startswith(_MARKER)
                    and (
                        len(buf) <= len(_MARKER)
                        or buf[len(_MARKER) : len(_MARKER) + 1] in b" \t\n"
                    )
                ):
                    self._offsets.append(0)
                # 検索開始位置（重複を避ける）
                search_from = 0 if chunk_start == 0 else overlap
                pos = 0
                while True:
                    idx = buf.find(_NEWLINE_MARKER, search_from + pos)
                    if idx == -1:
                        break
                    self._offsets.append(chunk_start + idx + 1)  # +\n の次
                    pos = idx + 1
                chunk_start += step
        print(
            f".vvol: {len(self._offsets)} フレーム（インデックス済）", file=sys.stderr
        )

    @property
    def offsets(self):
        self._ensure_index()
        return self._offsets

    def get_frame_data(self, frame_idx):
        """指定フレームの先頭にseekしてそのフレームだけ読み、処理済みデータを返す"""
        self._ensure_index()
        if frame_idx < 0 or frame_idx >= len(self._offsets):
            return None

        current_frame = None
        current_box = None
        current_centers = []
        current_volumes = []

        with open(self.filename, "r") as f:
            f.seek(self.offsets[frame_idx])
            for line in f:
                line = line.strip()
                if not line:
                    # 空行でフレーム終端
                    break
                if line.startswith("# FRAME"):
                    parts = line.split()
                    if len(parts) >= 3:
                        current_frame = int(parts[2])
                elif line.startswith("BOX"):
                    parts = line.split()
                    if len(parts) >= 4:
                        Lx, Ly, Lz = float(parts[1]), float(parts[2]), float(parts[3])
                        current_box = np.diag([Lx, Ly, Lz])
                else:
                    parts = line.split()
                    if len(parts) >= 4:
                        x, y, z, volume = (
                            float(parts[0]),
                            float(parts[1]),
                            float(parts[2]),
                            float(parts[3]),
                        )
                        current_centers.append([x, y, z])
                        current_volumes.append(volume)

        if current_box is None or not current_centers:
            return None

        positions = np.array(current_centers, dtype=np.float64)
        volumes = np.array(current_volumes, dtype=np.float64)
        devs = volumes - np.mean(volumes)

        return {
            "pos": positions,
            "dev": devs,
            "mean_v": np.mean(volumes),
            "box": current_box,
        }

    def __len__(self):
        return len(self.offsets)
