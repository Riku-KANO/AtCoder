import numpy as np
vecs = np.array([
    [ 1, 0, 0],
    [-1, 0, 0],
    [0,  1, 0],
    [0, -1, 0],
    [0, 0,  1],
    [0, 0, -1],
  ], dtype=int
  )

Rx = np.array([
    [1,0,0],  
    [0,0,-1],
    [0,1,0]
  ], dtype=int
)

Ry = np.array([
    [0,0,1],
    [0,1,0],
    [-1,0,0]
  ], dtype=int
)

Rz = np.array([
    [0,-1,0],
    [1,0,0],
    [0,0,1],
  ], dtype=int
)

memo = {}

for i in range(40000):
    r = np.random.randint(0, 3)
    if r == 0:
        vecs = vecs @ Rx
    elif r == 1:
        vecs = vecs @ Ry
    else:
        vecs = vecs @ Rz
    
    t = tuple(map(lambda x: tuple(x), vecs))
    memo[t] = 1

print(len(memo.keys()))

print("int dirs[24][3][6] = {")
for k in sorted(memo.keys())[::-1]:
    print("  {")
    for i in range(3):
        print("    {", end="")
        for d in range(6):
            print(k[d][i], end=",")
        print("},")
    print("  },")

print("};")

