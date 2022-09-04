s = input()
exist = [False] * 10
if s[0] == '1':
  print("No")
  exit(0)

for i in range(len(s)):
  num = str(s[i])
  if num:
    if i == 0 or i == 4: exist[3]=True
    if i == 1 or i == 7: exist[2]=True
    if i == 2 or i == 8: exist[4]=True
    if i == 3: exist[1]=True
    if i == 5: exist[5]=True
    if i == 6: exist[0]=True
    if i == 9: exist[6]=True

for i in range(7):
  for j in range(i + 1, 7):
    for k in range(i + 1, j):
      if exist[i] and exist[j] and not exist[k]:
        print("Yes")
        exit(0)

print("No")  