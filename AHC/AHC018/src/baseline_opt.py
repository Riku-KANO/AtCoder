# # memo
# C++のboostのLU分解で逆行列を求める関数を作ったが、
# numpy.linalg.invより遅すぎたのでPythonを使ってみようと思う。
# (200*200の行列の逆行列計算で、C++: 603ms, Python: 110ms @atcoder code test)
# boostのublasの線形代数演算が遅すぎるのだろうか。
# 
# ガウス過程回帰を使う。

import numpy as np
import time
import sys
from itertools import product
import optuna
import bisect
import glob

def load_input(file_name):
  X = []
  y = []
  waters = []
  houses = []
  S = []
  with open(file_name, "r") as f:
    s = iter(f.read().split("\n"))
    N, W, K, C = map(int, next(s).split())
    for i in range(N):
      ss = list(map(int, next(s).split()))
      S.append(ss)
      for j in range(N):
        X.append([i, j])
        y.append(np.log10(ss[j]))
    for i in range(W):
      waters.append(tuple(map(int, next(s).split())))
    for j in range(K):
      houses.append(tuple(map(int, next(s).split())))
  
  return N, W, K, C, S, waters, houses, X, y


DEBUG = False
if DEBUG:
  import matplotlib.pyplot as plt

start_time = None
cur_time = None

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


def solve(file_name, start, coef):
  N, W, K, C, S, sources, houses, X, y = load_input(file_name=file_name)
  pos = sources + houses
  dists = []
  for i in range(W+K):
    for j in range(i):
      dist = (pos[i][0]-pos[j][0])**2+(pos[i][1]-pos[j][1])**2
      dists.append((dist, j, i))
  dists.sort()

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


  powers = [start]
  power = start
  while sum(powers) <= 5000:
    power *= coef
    powers.append(power)

  cum_powers = [0]
  for power in powers:
    cum_powers.append(cum_powers[-1] + power)

  total_cost = 0
  for i in range(N):
    for j in range(N):
      if ans_grid[i][j]:
        idx = bisect.bisect_left(cum_powers, S[i][j])
        
        total_cost += min(cum_powers[idx], 5000) + C * idx
        
  return total_cost / L 


def objective(trial):
  start = trial.suggest_int('start', 10, 300)
  coef = trial.suggest_uniform("coef", 1, 2)
  total_cost = 0

  file_names = glob.glob("./in/C2/*/*/0000.txt")
  file_names2 = glob.glob("./in/C2/*/*/0001.txt")
  file_names = file_names + file_names2
  for file_name in file_names:
    total_cost += solve(file_name, start, coef)

  return total_cost

if __name__ == "__main__":
  start_time = time.time()
  
  study = optuna.create_study()
  study.optimize(objective, n_trials=100, n_jobs=4)
  print("C=2")
  print(study.best_params)