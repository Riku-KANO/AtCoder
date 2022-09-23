from typing import List

class UnionFind:
  def __init__(self, n: int) -> None:
    self.n: int = n
    self.par: List[int] = [-1] * n

  def merge(self, a: int, b: int) -> None:
    a_leader: int = self.leader(a)
    b_leader: int = self.leader(b)
    self.par[a_leader] = b_leader

  def same(self, a: int, b: int) -> bool:
    return self.leader(a) == self.leader(b)

  def leader(self, a: int) -> int:
    while self.par[a] != -1:
      a = par[a]
    return a
  
  def groups(self) -> List[List[int]]:
    ret: List[List[int]] = []
    memo = [[]] * self.n
    for i in range(self.n):
      if p == -1:
        ret.append(i)
      else:
        memo[p].append(i)
    for lst in memo:
      if len(lst) != 0:
        ret.append(lst)
    return ret
