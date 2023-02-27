# # memo
# C++のboostのLU分解で逆行列を求める関数を作ったが、
# numpy.linalg.invより遅すぎたのでPythonを使ってみようと思う。
# (200*200の行列の逆行列計算で、C++: 603ms, Python: 110ms @atcoder code test)
# boostのublasの線形代数演算が遅すぎるのだろうか。
# 
# ガウス過程回帰を使う。

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import time
import sys
from itertools import product

DEBUG = False
if DEBUG:
  import matplotlib.pyplot as plt

start_time = None
cur_time = None

alphas = [38, 42, 48, 49, 87, 109, 152, 234]
coefs = [
  1.071071836416865,
  1.1103374315829098,
  1.1581850602369728,
  1.2295262864675027,
  1.2309734139372803,
  1.4340920873990537,
  1.4505429558219907,
  1.7479564275184625
]

def get_time():
  cur_time = time.time()
  return start_time - cur_time

def excavate(i: int, j: int, power: int):
  print(f"{i} {j} {power}", flush=True)
  return int(sys.stdin.readline().strip())
  
class UnionFind:
  def __init__(self, n):
    self.sz: int = n
    self.par = [-1] * n

  def same(self, a: int, b: int) -> bool :
    return self.leader(a) == self.leader(b)

  def leader(self, a: int) -> int:
    if self.par[a] == -1:
      return a
    p = self.leader(self.par[a])
    self.par[a] = p
    return p

  def merge(self, a: int, b: int):
    pa = self.leader(a)
    self.par[pa] = self.leader(b)


def solve():
  N, W, K, C = map(int, sys.stdin.readline().split())
  if DEBUG:
    for _ in range(N): input()
  sources = [tuple(map(int, sys.stdin.readline().split())) for _ in range(W)]
  houses = [tuple(map(int, sys.stdin.readline().split())) for _ in range(K)]
  pos = sources + houses
  dists = []
  for i in range(W+K):
    for j in range(i):
      dist = (pos[i][0]-pos[j][0])**2+(pos[i][1]-pos[j][1])**2
      dists.append((dist, j, i))
  dists.sort()

  print("Init Done", file=sys.stderr)

  uf = UnionFind(N)
  boss = [-1] * (W+K)
  edges = []
  for t in dists:
    dist, i, j = t
    if uf.same(i, j) or (i < W and j < W):
      continue
    if i < W and boss[j] != -1:
      continue
    if boss[i] != -1 and boss[j] != -1 and boss[i] != boss[j]:
      continue
    uf.merge(i, j)
    edges.append((i, j))
    for w in range(W):
      for k in range(K):
        if uf.same(w, k+W):
          boss[W+k] = w

  M = len(edges)  
  L = 10**18
  ans_grid = None
  def fill(ps, pt, g):
    if ps[0] == pt[0]:
      for j in range(min(ps[1], pt[1]), max(ps[1], pt[1]) + 1):
        g[ps[0]][j] = 1
    else:
      for i in range(min(ps[0], pt[0]), max(ps[0], pt[0]) + 1):
        g[i][ps[1]] = 1

      
  for i in range(1<<M):
    l = 0
    grid = np.zeros((N, N), dtype=np.int8)
    for j, edge in enumerate(edges):
      s, t = edge
      ps = pos[s]
      pt = pos[t]
      if (1 << j) & i:
        pu = (pt[0], ps[1])
      else:
        pu = (ps[0], pt[1])
        
      fill(ps, pu, grid)
      fill(pu, pt, grid)

    l = np.sum(grid)
    if l < L:
      L = l
      ans_grid = grid

  print("Found MST", file=sys.stderr)

  X_train = []
  y_train = []

  # calc powers
  powers = []
  for i in range(8):
    if 2 ** i == C:
      cur = alphas[i]
      
      while cur < 5000:
        powers.append(int(cur)) 
        cur *= coefs[i]
      powers.append(5000)

  power_list = [alphas[i]]
  for i in range(len(powers) - 1):
    power_list.append(powers[i+1]-powers[i])

  for i in range(N):
    for j in range(N):
      if (i + j) % 3 == 0 and ans_grid[i, j] == 1:
        res = 0
        
        total_power = 0
        for power in power_list:
          res = excavate(i, j, power)
          total_power += power
          if res == 1 or res == 2:
            break
          
        print(f"added: {i} {j} {total_power-power / 2}", file=sys.stderr)
        X_train.append((i, j))
        y_train.append(total_power - power / 2)
  
  print("Got train data", file=sys.stderr)

  X_train = np.array(X_train)
  y_train = np.array(y_train)
  X_all = np.array([p for p in product(range(N), range(N))])
  kernel = RBF()
  gp = GaussianProcessRegressor(kernel=kernel)
  gp.fit(X_train, y_train)
  y_pred = gp.predict(X_all).reshape(N, N)
  
  print(f"sample predict: {y_pred[20, 20]}, {y_pred[70, 29]}, {y_pred[169, 195]}", file=sys.stderr)
  print("train done", file=sys.stderr)

  for i in range(N):
    for j in range(N):
      if (i + j) % 3 == 0 or ans_grid[i, j] == 0:
        continue

      pred = max(min(int(y_pred[i, j] + 10), 5000), 10)
      print(f"{i} {j} {pred}", file=sys.stderr)
      res = excavate(i, j, pred)
      if res == 2:
        print("All done", file=sys.stderr)
        return
      elif res == 1:
        continue
      elif res == -1:
        print("Invalid output", file=sys.stderr)
        exit(1)
      else: # res == 0
        power = C * 5
        while res == 0:
          res = excavate(i, j, power)


if __name__ == "__main__":
  start_time = time.time()
  solve()
