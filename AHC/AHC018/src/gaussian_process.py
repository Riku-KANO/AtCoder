import numpy as np

class GaussianProcessRegressor:
  def __init__(self, kernel, noise=0.0):
    self.kernel = kernel
    self.noise = noise
    self.X_train = None
    self.y_train = None
    self.alpha = None

  def fit(self, X_train, y_train):
    self.X_train = X_train
    self.y_train = y_train
    K = self.kernel(X_train, X_train) + np.eye(len(X_train)) * self.noise
    self.alpha = np.linalg.solve(K, y_train)

  def predict(self, X_test):
    K_s = self.kernel(self.X_train, X_test)
    K_ss = self.kernel(X_test, X_test)
    y_mean = K_s.T.dot(self.alpha)
    y_cov = K_ss - K_s.T.dot(np.linalg.solve(self.kernel(self.X_train, self.X_train) + np.eye(len(self.X_train)) * self.noise, K_s))
    return y_mean, y_cov
