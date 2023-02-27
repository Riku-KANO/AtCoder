import os
import optuna
import matplotlib.pyplot as plot
import glob
import threading
import re
import pickle
C = 1
DEBUG = False


def objective(trial):
  total_score = 0
  thred_id = threading.get_ident()
  file_names = glob.glob(f"./in/default/*")
  
  start = trial.suggest_int("start", 10, 500)
  coef = trial.suggest_uniform("coef", 1., 4.)
  init_stiff_exp = trial.suggest_int("stiff_init", 10, 1000)

  if DEBUG:
    file_names = file_names[:5]
  

  for file_name in file_names:
    os.system(f"./bin/solve3_opt --start {start} --coef {coef} --stiff {init_stiff_exp} --C {C} < {file_name} 1> ./out/opt/solve2/out{thred_id}.txt 2> ./out/opt/solve2/debug{thred_id}.txt")
    with open(f"./out/opt/solve3/debug{thred_id}.txt", "r") as f:
      s = f.read()

    score_pattern = r"SCORE:\s*(\d+)"
    match = re.search(score_pattern, s)
    if match:
      score = int(match.group(1))
    else:
      print("SCORE not found")
      score = 10 ** 8
    total_score += score

  return total_score

def objective2(trial):
  total_score = 0
  thred_id = threading.get_ident()
  file_name = "./in/default/0000.txt"
  
  start = trial.suggest_int("start", 10, 800)
  coef = trial.suggest_uniform("coef", 1., 2.5)

  num_iter = 20 if not DEBUG else 5  

  for i in range(num_iter):
    os.system(f"./bin/solve3_opt --seed {i} --start {start} --coef {coef} --C {C} < {file_name} 1> ./out/opt/solve3/out{thred_id}.txt 2> ./out/opt/solve3/debug{thred_id}.txt")
    with open(f"./out/opt/solve3/debug{thred_id}.txt", "r") as f:
      s = f.read()

    score_pattern = r"SCORE:\s*(\d+)"
    match = re.search(score_pattern, s)
    if match:
      score = int(match.group(1))
    else:
      print("SCORE not found")
      score = 10 ** 8
    total_score += score

  return total_score

if __name__ == "__main__":
  os.system("g++ ./src/solve3.cpp -O2 -o ./bin/solve3_opt -D_LOCAL -D_OPTUNA -std=c++17")
  for i in range(1):
    
    i = 4
    C = 2 ** i
    study = optuna.create_study()
    study.optimize(objective2, n_trials=100, n_jobs=2)
    with open(f"study{2**i}.pkl", "wb") as f:
      pickle.dump(study, f)
    