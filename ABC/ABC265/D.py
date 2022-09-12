from itertools import accumulate
N, P, Q, R = map(int, input().split())
A = list(map(int, input().split()))
cum = [0]+list(accumulate(A))
s = set(cum)
for i in range(N):
  if (cum[i]+P) in s and (cum[i]+P+Q) in s and (cum[i]+P+Q+R) in s:
    print("Yes")
    exit(0)
print("No")