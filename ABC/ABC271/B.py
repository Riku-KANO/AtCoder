N, Q = map(int, input().split())
lst = []
for i in range(N):
  buf = list(map(int, input().split()))
  L = buf[0]
  buf.pop(0)
  lst.append(buf)

for q in range(Q):
  s, t = map(int,input().split())
  s -= 1
  t -=1
  print(lst[s][t])