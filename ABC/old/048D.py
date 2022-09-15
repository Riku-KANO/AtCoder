s = input()
n = len(s)
print("Second" if (n % 2 == 0) ^ (s[0] == s[-1]) else "First")