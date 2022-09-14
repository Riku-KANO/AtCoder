# DPでは間に合わない
# combination nCr
# フェルマーの小定理
# inv MOD

H, W, A, B = map(int, input().split())
MOD = int(1e9+7)
factorial = [0] * (H+W+1)
inv = [0] * (H+W+1)

for i in range(0, H + W + 1):
  if i == 0:
    factorial[i] = 1
    inv[i] = pow(factorial[i], MOD-2, MOD)
  else:
    factorial[i] = factorial[i-1] * i % MOD
    inv[i] = pow(factorial[i], MOD-2, MOD)

def comb(n, r):
  return factorial[n] * inv[r] * inv[n-r] % MOD

ans = 0
for i in range(B, W):
  ans += comb(H-A+i-1, i) * comb(W-i+A-2, A-1)
  ans %= MOD
  
print(ans)