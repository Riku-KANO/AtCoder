import os
import re

EXE_MAIN = "./bin/solve"
REF = "out/sub4/"

for i in range(100):
  print(f"### CASE {i:03d} ######")
  os.system(f"{EXE_MAIN} < in/{i:04d}.txt 1> out/simple/{i:04d}.txt 2> out/simple/{i:04d}debug.txt")
  with open(f"./out/simple/{i:04d}debug.txt", "r") as f:
    s_target = f.read()
  with open(f"{REF}/{i:04d}debug.txt", "r") as f:
    s_ref = f.read()
  score_pattern = r".*Score\s*:\s*(\d+)"
  match_target = re.search(score_pattern, s_target)
  match_ref = re.search(score_pattern, s_ref)

  if match_target:
    score_target = int(match_target.group(1))
  else:
    score_target = -1
    print("SCORE not found")
  if match_ref:
    score_ref = int(match_ref.group(1))
  else:
    score_ref = -1
    print("SCORE not found")


  print(f"score(main/ref): ({score_target}, {score_ref}), ratio: {score_target / score_ref}")