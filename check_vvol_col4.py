#!/usr/bin/env python3
"""00002.vvol の第4カラム（体積）の特異値・外れ値をチェックする"""
import sys
from collections import Counter

path = "/Volumes/Caches/Dropbox/00Done/gitbox/easy-voro/00002.vvol"

n = 0
total = 0.0
min_v = float("inf")
max_v = float("-inf")
negatives = []
zeros = 0
non_numeric = []
# 外れ値検出用: 1パスで平均は出ないので、まず全データを読んでから統計（メモリ注意）
# 代わりにヒストグラム用バケットで分布を見る
buckets = [0] * 20  # 0-0.05, 0.05-0.1, ... など
very_small = []   # 0.01未満
very_large = []   # 0.5以上

try:
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("BOX"):
                continue
            parts = line.split()
            if len(parts) != 4:
                non_numeric.append((n, line[:80]))
                continue
            try:
                v = float(parts[3])
            except ValueError:
                non_numeric.append((n, line[:80]))
                continue
            n += 1
            total += v
            if v < min_v:
                min_v = v
            if v > max_v:
                max_v = v
            if v < 0:
                negatives.append((n, v))
            if v == 0:
                zeros += 1
            if v < 0.01 and len(very_small) < 20:
                very_small.append((n, v))
            if v > 0.5 and len(very_large) < 20:
                very_large.append((n, v))
            # バケット 0~0.5 を 0.025 刻みで
            idx = int(v / 0.025)
            if 0 <= idx < 20:
                buckets[idx] += 1
            elif idx >= 20:
                buckets[19] += 1
except KeyboardInterrupt:
    print("中断されました", file=sys.stderr)

if n == 0:
    print("データ行が0件でした")
    sys.exit(1)

mean = total / n
print("=== 第4カラム（体積）基本統計 ===")
print(f"有効データ行数: {n}")
print(f"最小値: {min_v}")
print(f"最大値: {max_v}")
print(f"平均値: {mean}")

print("\n=== 特異値の有無 ===")
if negatives:
    print(f"負の値: {len(negatives)} 件")
    for row, val in negatives[:10]:
        print(f"  行目付近 {row}: {val}")
else:
    print("負の値: なし")

print(f"ゼロの個数: {zeros}")

if non_numeric:
    print(f"非数値・カラム不足: {len(non_numeric)} 件")
    for row, s in non_numeric[:5]:
        print(f"  {row}: {s!r}")
else:
    print("非数値・カラム不足: なし")

if very_small:
    print(f"\n0.01 未満の値の例（最大20件）: {len(very_small)} 件")
    for row, val in very_small[:10]:
        print(f"  {row}: {val}")
if very_large:
    print(f"\n0.5 以上の値の例（最大20件）: {len(very_large)} 件")
    for row, val in very_large[:10]:
        print(f"  {row}: {val}")

print("\n=== 分布（0〜0.5 を 0.025 刻み） ===")
for i in range(20):
    lo, hi = i * 0.025, (i + 1) * 0.025
    print(f"  [{lo:.3f}, {hi:.3f}): {buckets[i]}")
