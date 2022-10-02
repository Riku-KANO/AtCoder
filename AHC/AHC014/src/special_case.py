from itertools import product
n = 31
xargs = list(range(n//4, n*3//4+1))

def manhattan(x, y, a, b):
  return abs(x-a) + abs(y-b)

lst = []

def pattern1():
  for x in product(xargs, xargs):
    if manhattan(x[0], x[1], (n-1)//2, (n-1)//2) % 5 == 0:
      lst.append(x)




def pattern2():
  for x in product(xargs, xargs):
    if (x[0]//2 + x[1]//2) % 2 == 0:
      lst.append(x)

if __name__ == "__main__":
  pattern2()

  with open("../input/special/special009.txt", "w") as f:
    f.write(f"{n} {len(lst)}\n")
    for x in lst:
      f.write(f"{x[0]} {x[1]}\n")