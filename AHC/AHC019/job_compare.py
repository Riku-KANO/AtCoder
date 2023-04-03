import os
import re

EXE_MAIN = "./bin/solve7"
REF = "./out/sub10/"

TL = 6.

for i in range(100):
  print(f"### CASE {i:03d} ######")
  os.system(f"{EXE_MAIN} --TL {TL} < in/{i:04d}.txt 1> out/tmp/out{i:04d}.txt 2> out/tmp/debug{i:04d}.txt")
  with open(f"./out/tmp/debug{i:04d}.txt", "r") as f:
    s_target = f.read()
  with open(f"{REF}/debug{i:04d}.txt", "r") as f:
    s_ref = f.read()
  score_pattern = r"SCORE\s*:\s*(\d+)"
  time_pattern = r"TIME\s*:\s*(\d+\.\d+)\ss"
  match_target = re.search(score_pattern, s_target)
  match_ref = re.search(score_pattern, s_ref)
  match_time = re.search(time_pattern, s_target)

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

  if match_time:
    elapsed_time = float(match_time.group(1))
  else:
    elapsed_time = 0

  print(f"score(main/ref): ({score_target}, {score_ref}), ratio: {score_ref / score_target}")
  print(f"elapsed time: {elapsed_time} s")
  if elapsed_time > TL:
    print("########################")
    print("   TIME LIMIT OVER !!!!!")
    print("########################\n")
  
  if not os.path.exists(f"./out/best/seed{i:04d}.txt"):
    with open(f"./out/best/seed{i:04d}.txt", "w") as f:
      f.write(f"{score_target}\n")
    
    os.system(f"cat out/tmp/out{i:04d}.txt >> out/best/seed{i:04d}.txt")
  
  else:
    with open(f"./out/best/seed{i:04d}.txt") as f:
      s = f.readlines()
      if s != []:
        best_score = int(s[0].strip())
        if best_score == -1:
          best_score = 1e19
      else:
        best_score = 1e19

    if score_target < best_score:
      with open(f"./out/best/seed{i:04d}.txt", "w") as f:
        f.write(f"{score_target}\n")
      os.system(f"echo {score_target} > out/best/seed{i:04d}.txt")
      os.system(f"cat out/tmp/out{i:04d}.txt >> out/best/seed{i:04d}.txt")


  