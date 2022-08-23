import numpy as np
N, M = map(int, input().split())
Adj = np.zeros(shape=(N+1, N+1), dtype=int)
for _ in range(M):
    U, V = map(int, input().split())
    Adj[U][V]=1
print(sum(Adj[i][j]*Adj[i][k]*Adj[j][k] for i in range(N+1) for j in range(N+1) for k in range(N+1)))