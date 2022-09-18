import sys
from collections import defaultdict

sys.setrecursionlimit(10000)

TEST = False

def test():
  N = 1000
  p = []
  for i in range(50):
    for j in range(20):
      p.append((i, j))
  return N, p

class UnionFind:
  def __init__(self, n):
    self.n = n
    self.parent = [x for x in range(n)]

  def same(self, a, b):
    pass

  def merge(self, a, b):
    par = self.parent[b]
    while self.parent[par] != par:
      par = self.parent[par]
      print(par)
    self.parent[par] = self.leader(a)

  def leader(self, a):
    par = self.parent[a]
    while par != self.parent[par]:
      par = self.parent[par]
    return par

  def get_num_group(self):
    pass

dx = [-1, -1, 0, 0, 1, 1]
dy = [-1, 0, -1, 1, 0, 1]

def main(n, p):
  N = int(input())
  exist = defaultdict(int)
  points = [tuple(map(int, input().split())) for _ in range(N)]
  if n is not None:
    N = n
    points = p

  uf = UnionFind(N + 10)
  for i in range(N):
    exist[points[i]] = i + 1
    x = points[i][0]
    y = points[i][1]
    for d in range(6):
      nx = x + dx[d]
      ny = y + dy[d]
      adi = exist[(nx, ny)]
      if adi > 0:
        uf.merge(i + 1, adi)

  s = set()
  for i in range(1, N+1):
    s.add(uf.leader(i))

  print(len(s))

if __name__ == "__main__":
  N = None
  P = None
  if TEST:
    N, P = test()

  main(N, P)
  # test()