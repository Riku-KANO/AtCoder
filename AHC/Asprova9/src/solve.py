import optuna

X, M, C, E = map(int, input().split())
costs = [[list(map(int, input().split())) for i in range(9)] for m in range(M)]

def objective(trial):
  patterns = []
  for i in range(M):
    exec(f"ma_{i} = trial.suggest_int('ma_{i}', 1, 9)")
    exec(f"mb_{i} = trial.suggest_int('mb_{i}', 1, 9)")
    exec(f"patterns.append((ma_{i}, mb_{i}))")
  for m in range(M):
    print("".join(map(str, patterns[m])) * X)
    
  score, V, D = map(int, input().split())
  for _ in range(X*M): input()
  return - score

study = optuna.create_study()
study.optimize(objective, n_trials=E)
