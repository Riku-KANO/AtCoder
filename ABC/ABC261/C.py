from collections import Counter, defaultdict
N = int(input())
cnt = Counter()
for i in range(N):
    s = input()
    if cnt[s] == 0:
        print(s)
    else:
        print(f"{s}({cnt[s]})")
    cnt[s]+=1