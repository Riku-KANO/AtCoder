n = int(input())
s = [input() for _ in range(n)]
ok = True
for i in range(n):
    for j in range(i+1, n):
        if s[i][j] == 'W' and s[j][i] != 'L': ok = False
        elif s[i][j] == 'L' and s[j][i] != 'W': ok = False
        elif s[i][j] == 'D' and s[j][i] != 'D': ok = False

        if not ok:
            print("incorrect")
            exit(0)

print("correct")