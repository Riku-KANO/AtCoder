import sys
import numpy as np
import copy
from collections import namedtuple
from itertools import product
from dataclasses import dataclass
from typing import Optional, List

di = [1, -1, 0, 0, 0, 0]
dj = [0, 0, 1, -1, 0, 0]
dk = [0, 0, 0, 0, 1, -1]

def is_out(x, y, z, d):
    return (x<0 or x>=d or y<0 or y>=d or z<0 or z>=d)

@dataclass
class Silhouette:
    f: List[List[int]]
    r: List[List[int]]
    h_f: List[List[int]]
    h_r: List[List[int]]
    exist: List[List[List[int]]]
    adj: List[List[List[int]]]
    vol: int

def load_input():
    D = int(input())
    f = [[] for i in range(2)]
    r = [[] for i in range(2)]
    for i in range(2):
        for k in range(D):
            f[i].append(input())
        for k in range(D):
            r[i].append(input())

    return D, f, r

def naive_solve1():
    D, f, r = load_input()
    b = [[0 for j in range(D * D * D)] for i in range(2)]
    n = 0
    for i in range(2):
        for x in range(D):
            for y in range(D):
                for z in range(D):
                    if f[i][z][x] == '1' and r[i][z][y] == '1':
                        n += 1
                        b[i][x * D * D + y * D + z] = n

    print(n)
    print(' '.join(map(str, b[0])))
    print(' '.join(map(str, b[1])))

    vol1 = 0
    vol2 = 0
    for v in b[0]:
        if v > 0:
            vol1 += 1

    for v in b[1]:
        if v > 0:
            vol2 += 1

    print(f"Volume 1: {vol1}", file=sys.stderr)
    print(f"Volume 2: {vol2}", file=sys.stderr)


def naive_solve2():
    D, f, r = load_input()
    f1 = [[0 for _ in range(D)] for _ in range(D)]
    r1 = [[0 for _ in range(D)] for _ in range(D)]
    f2 = [[0 for _ in range(D)] for _ in range(D)]
    r2 = [[0 for _ in range(D)] for _ in range(D)]
    h_f1 = [[0 for _ in range(D)] for _ in range(D)]
    h_r1 = [[0 for _ in range(D)] for _ in range(D)]
    h_f2 = [[0 for _ in range(D)] for _ in range(D)]
    h_r2 = [[0 for _ in range(D)] for _ in range(D)]
    exist1 = [[[0 for _ in range(D)] for _ in range(D)] for _ in range(D)]
    exist2 = [[[0 for _ in range(D)] for _ in range(D)] for _ in range(D)]
    res1 = [0 for _ in range(D*D*D)]
    res2 = [0 for _ in range(D*D*D)]
    vol1 = 0
    vol2 = 0
    adj1 = [[[0 for _ in range(D)] for _ in range(D)] for _ in range(D)]
    adj2 = [[[0 for _ in range(D)] for _ in range(D)] for _ in range(D)]

    for i in range(D):
        for j in range(D):
            f1[i][j] = 1 if f[0][i][j] == '1' else 0
            r1[i][j] = 1 if r[0][i][j] == '1' else 0
            f2[i][j] = 1 if f[1][i][j] == '1' else 0
            r2[i][j] = 1 if r[1][i][j] == '1' else 0
    
    for i, (x, y, z) in enumerate(product(range(D), range(D), range(D))):
        if f1[z][x] and r1[z][y]:
            h_f1[z][x] += 1
            h_r1[z][y] += 1
            exist1[x][y][z] = 1
            vol1 += 1

        if f2[z][x] and r2[z][y]:
            h_f2[z][x] += 1
            h_r2[z][y] += 1
            exist2[x][y][z] = 1
            vol2 += 1

    # check adj
    for i, (x,y,z) in enumerate(product(range(D), range(D), range(D))):
        if exist1[x][y][z]:
            for d in range(6):
                nx=x+di[d]
                ny=y+dj[d]
                nz=z+dk[d]
                if is_out(nx,ny,nz,D): continue
                if exist1[nx][ny][nz]:
                    adj1[x][y][z]+=1

        if exist2[x][y][z]:
            for d in range(6):
                nx=x+di[d]
                ny=y+dj[d]
                nz=z+dk[d]
                if is_out(nx,ny,nz,D): continue
                if exist2[nx][ny][nz]:
                    adj2[x][y][z] += 1

    sil1 = Silhouette(f1, r1, h_f1, h_r1, exist1, adj1, vol1)
    sil2 = Silhouette(f2, r2, h_f2, h_r2, exist2, adj2, vol2)
    print("Before:", file=sys.stderr)
    print(f"Volume1: {vol1}", file=sys.stderr)
    print(f"Volume2: {vol2}", file=sys.stderr)

    # isolate check
    for i, (x, y, z) in enumerate(product(range(D), range(D), range(D))):
        if sil1.exist[x][y][z] and sil1.adj[x][y][z] == 0:
            sil1.exist[x][y][z] = 0
            sil1.h_f[z][x] -= 1
            sil1.h_r[z][y] -= 1
            sil1.vol -= 1


        if sil2.exist[x][y][z] and sil2.adj[x][y][z] == 0:
            sil2.exist[x][y][z] = 0
            sil2.h_f[z][x] -= 1
            sil2.h_r[z][y] -= 1
            sil2.vol -= 1

    maxIter = 300
    iter = 0
    if sil1.vol > sil2.vol:
        while sil1.vol != sil2.vol:
            iter += 1
            if maxIter==0:
                break
            maxIter-=1
            pos=(-1,-1,-1)
            best_score = 0
            for i, (x,y,z) in enumerate(product(range(D), range(D), range(D))):
                if sil1.h_f[z][x] == 1 or sil1.h_r[z][y] == 1: continue
                if sil1.exist[x][y][z] == 0: continue

                ok = True
                for d in range(6):
                    nx=x+di[d]
                    ny=y+dj[d]
                    nz=z+dk[d]
                    if is_out(nx,ny,nz,D): continue
                    if sil1.adj[nx][ny][nz]==1 and (sil1.h_f[nz][nx] <= 2 or sil1.h_r[nz][ny] <= 2):
                        ok = False
                        break
                
                if not ok: continue
                
                score = sil1.h_f[z][x] + sil1.h_r[z][y]
                if sil1.adj[x][y][z] == 1 and sil1.h_f[z][x] > 1 and sil1.h_r[z][y] > 1:
                    print(f"Iteration {iter}: {x},{y},{z} will be deleted", file=sys.stderr)
                    score = 200 + abs((D-1)/2-x) + abs((D-1)/2-y) + abs((D-1)/2-z)


                if score > best_score:
                    best_score = score
                    pos = (x,y,z)

            if pos == (-1,-1,-1): continue
            x,y,z=pos
            print(f"Iteration {iter}: {x}, {y}, {z} is deleted!!.. (score: {best_score})", file=sys.stderr)
            sil1.exist[x][y][z] = 0
            sil1.h_f[z][x] -= 1
            sil1.h_r[z][y] -= 1
            sil1.vol -= 1
            for d in range(6):
                nx=x+di[d]
                ny=y+dj[d]
                nz=z+dk[d]
                if is_out(nx,ny,nz,D): continue
                if sil1.exist[nx][ny][nz]:
                    sil1.adj[nx][ny][nz] -= 1

    elif sil1.vol < sil2.vol:
        while sil1.vol != sil2.vol:
            iter += 1
            if maxIter==0:
                break
            maxIter-=1
            pos = (-1,-1,-1)
            best_score = 0
            for i, (x,y,z) in enumerate(product(range(D), range(D), range(D))):
                if sil2.h_f[z][x] == 1 or sil2.h_r[z][y] == 1: continue
                if sil2.exist[x][y][z] == 0: continue

                ok = True
                for d in range(6):
                    nx=x+di[d]
                    ny=y+dj[d]
                    nz=z+dk[d]
                    if is_out(nx,ny,nz,D): continue
                    if sil2.adj[nx][ny][nz] == 1:
                        ok = False
                        break
                
                if not ok: continue

                score = sil2.h_f[z][x] + sil2.h_r[z][y]
                if sil2.adj[x][y][z] == 1 and sil2.h_f[z][x] > 1 and sil2.h_r[z][y] > 1:
                    print(f"Iteration {iter}: {x},{y},{z} will be deleted", file=sys.stderr)
                    score = 200 + abs((D-1)/2-x) + abs((D-1)/2-y) + abs((D-1)/2-z)

                if score > best_score:
                    best_score = score
                    pos = (x,y,z)
            
            if pos == (-1,-1,-1): continue
            x,y,z=pos

            print(f"Iteration {iter}: {x}, {y}, {z} is deleted!!.. (score: {best_score})", file=sys.stderr)

            sil2.exist[x][y][z] = 0
            sil2.h_f[z][x] -= 1
            sil2.h_r[z][y] -= 1
            sil2.vol -= 1
            for d in range(6):
                nx=x+di[d]
                ny=y+dj[d]
                nz=z+dk[d]
                if is_out(nx,ny,nz,D): continue
                if sil2.exist[nx][ny][nz]:
                    sil2.adj[nx][ny][nz] -= 1

    n_ans = 0
    cid = 1
    for i, (x,y,z) in enumerate(product(range(D), range(D), range(D))):
        if sil1.exist[x][y][z]:
            res1[i] = cid
            cid += 1
    n_ans = max(n_ans, cid-1)
    cid = 1
    for i, (x,y,z) in enumerate(product(range(D), range(D), range(D))):
        if sil2.exist[x][y][z]:
            res2[i] = cid
            cid += 1
    n_ans = max(n_ans, cid-1)

    print(n_ans)
    print(*res1)
    print(*res2)

    print(f"max ID: {n_ans}", file=sys.stderr)
    print(f"Volume 1: {sil1.vol}", file=sys.stderr)
    print(f"Volume 2: {sil2.vol}", file=sys.stderr)

    

if __name__ == "__main__":
    naive_solve1()