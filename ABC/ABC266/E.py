N = int(input())
ans = 3.5
for i in range(N-1):
  ans = ans * int(ans) / 6 + sum(range(int(ans) + 1, 7)) / 6
print(ans)