import numpy as np
from itertools import product
import argparse
import pathlib

def get_arg():
  parser = argparse.ArgumentParser()
  parser.add_argument("--size", default = 10, type=int)
  parser.add_argument("-n", default=15, type=int)
  parser.add_argument("-p", "--path", default=".")
  args = parser.parse_args()
  return args

if __name__ == "__main__":
  args = get_arg()
  n = args.n
  size = args.size
  path = pathlib.Path(args.path)
  for i in range(size):
    m = np.random.randint(n, n**2//12 + 1)
    x_range = list(range(n//4, 3*n//4+1))
    idx = list(range(len(x_range) * len(x_range)))
    target_idx = np.random.choice(idx, m, replace=False)
    points = [x for x in product(range(n//4, 3*n//4 + 1), range(n//4, 3*n//4 + 1))]

    with open(path / f"in{i:04d}.txt", "w") as f:
      f.write(f"{n} {m}\n")
      for j in range(m):
        point = points[target_idx[j]]
        f.write("{} {}\n".format(point[0], point[1]))
