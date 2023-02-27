// debug code for gaussian process
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>
#include <random>
#include <queue>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <utility>
#include <iomanip>
#include <iterator>
#include <stack>
#include <cmath>
#include <numeric>
#include <chrono>
#include <sstream>
#include <optional>

// atcoder
#include <atcoder/dsu>

// boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#define pii std::pair<int,int>
#define Vec std::vector<int>
#define rep(i, s, n) for(int i = s; i < (n); i++)
#define Grid std::vector<std::vector<bool>> 
using namespace std;
using namespace atcoder;
namespace ublas = boost::numeric::ublas;
using ll = long long;

const int INF = 1e9 + 10;
const ll LINF = 1LL<<60;
const double EPS = 1e-6;
constexpr double PI = std::acos(-1);

// clockwise
static const int di[4] = {0, 1, 0, -1};
static const int dj[4] = {1, 0, -1, 0};
static const int di8[8] = {0, 1, 1,  1,  0, -1, -1, -1};
static const int dj8[8] = {1, 1, 0, -1, -1, -1,  0,  1};

std::random_device seed_gen;
std::mt19937 mt(seed_gen());

clock_t start_time, cur_time;

double get_time() {
  cur_time = clock();
  return (double)(cur_time - start_time) / CLOCKS_PER_SEC;
}

template<typename T> T pow2(T x) {return x * x;}

bool is_out(int nx, int ny, int H, int W) {
  return (nx < 0 or ny < 0 or nx >= H or ny >= W);
}

namespace gl {
  int N; // 
  int W; // 水源
  int K;
  int C;
  int wi[5];
  int wj[5];
  int hi[10];
  int hj[10];
  int S[404][404];
  const int MAXP = 5'000;
  const int MINP = 10;
}

namespace param {
  const double TL = 5.0;
  const double TL_90 = TL * 0.90;
  const double TL_95 = TL * 0.95;
  int seed;
  double COEF = 1.0;
  double LENGTH = 7.0;
}

// calc invert matrix by LU decomposition
template<typename T> ublas::matrix<T> calc_invmat(ublas::matrix<T> A){
  ublas::matrix<T> B = ublas::identity_matrix<T>(A.size1());
  ublas::permutation_matrix<std::size_t> pm(A.size1());
  int res = ublas::lu_factorize(A, pm);
  if(res != 0) {
    std::cerr << "calc inverse Failed" << std::endl;
    return B;
  }
  ublas::lu_substitute(A, pm, B);
  return B;
}

class GaussianProcess {
public:
  GaussianProcess(double coef = 1.0, double length = 5.0) {
    this->coef = coef;
    this->length = length;
    this->is_trained = false;
  }

  void fit(const std::vector<std::pair<int,int>>& pos, const std::vector<double>& y, const std::vector<double>& noises) {
    std::cerr << "start fit" << std::endl;
    size_t ker_size = pos.size();
    this->X_train = pos;
    this->y_train = y;
    this->noises = noises;
    ublas::matrix<double> y_vec(ker_size, 1);
    ublas::matrix<double> k_train(ker_size, ker_size);
    double y_train_mean_ = 0.0;
    double y_train_std_ = 0.0;

    rep(i, 0, ker_size) {
      y_vec(i, 0) = y[i];
      y_train_mean_ += y[i];
    }
    y_train_mean_ /= (int)ker_size;
    this->y_train_mean = y_train_mean_;

    
    rep(i, 0, ker_size) y_train_std_ += pow2(y[i]-y_train_mean_);
    
    y_train_std_ = std::sqrt(y_train_std_ / ker_size);
    this->y_train_std = y_train_std_;
    rep(i, 0, ker_size) {
      rep(j, 0, ker_size) {
        k_train(i, j) = this->rbf_kernel(pos[i], pos[j]);
      }
    }

    this->K = k_train;
    rep(i, 0, ker_size) k_train(i, i) += pow2(noises[i]);
    ublas::matrix<double> k_train_inv = calc_invmat(k_train);
    this->K_inv = k_train_inv;
    ublas::matrix<double> alp = ublas::prod(k_train_inv, y_vec);
    this->alpha = alp;
    this->is_trained = true;
    // std::cerr << "alpha: " << this->alpha << std::endl;
    // std::cerr << "K: " << this->K << std::endl;
    // std::cerr << "K inv: " << this->K_inv << std::endl;
    std::cerr << "y train mean: " << this->y_train_mean << std::endl;
    std::cerr << "y train std: " << this->y_train_std << std::endl;
  }

  std::pair<double, double> predict(const pii& x){
    if(! this->is_trained) {
      std::cerr << "This instance is not fitted. please fit before predict" << std::endl;
      return {-1, -1};
    }
    size_t ker_size = this->X_train.size();
    std::pair<double, double> ret;
    ublas::matrix<double> k_s(ker_size, 1);
    double k_ss;
    rep(i, 0, ker_size) k_s(i, 0) = this->rbf_kernel(this->X_train[i], x);
    k_ss = this->rbf_kernel(x, x);
    
    ublas::matrix<double> y_tmp = ublas::prod(ublas::trans(k_s), this->alpha);
    
    double y_mean = y_tmp(0,0);
    // std::cerr << y_mean << std::endl;
    ublas::matrix<double> k_t = ublas::trans(k_s);
    ublas::matrix<double> k_2 = ublas::prod(this->K_inv, k_s);
    // std::cerr << k_t << std::endl;
    // std::cerr << k_2 << std::endl;
    double y_var = k_ss - ublas::prod(k_t , k_2)(0, 0);
    ret = {y_mean, y_var};
    return ret;
  }

  std::vector<std::pair<double, double>> predict(const std::vector<pii>& x){
    if(! this->is_trained) {
      std::cerr << "This instance is not fitted. please fit before predict" << std::endl;
      std::vector<pair<double,double>> ret;
      return ret;
    }
    typedef std::pair<double, double> Result;
    size_t ker_size = this->X_train.size();
    size_t test_size = x.size();
    std::vector<Result> ret;
    ublas::matrix<double> k_s(ker_size, test_size);
    ublas::matrix<double> k_ss(test_size, test_size);
    rep(i, 0, ker_size) {
      rep(j, 0, test_size) {
        k_s(i, j) = this->rbf_kernel(this->X_train[i], x[j]);
      }
    }
    rep(i, 0, test_size) {
      rep(j, 0, test_size) {
        k_ss(i, j) = this->rbf_kernel(x[i], x[j]);
      }
    }
    ublas::matrix<double> y_mean = ublas::prod(this->alpha, k_s);
    ublas::matrix<double> k_t = ublas::trans(k_s);
    ublas::matrix<double> k = ublas::prod(this->K_inv, k_s);
    ublas::matrix<double> y_cov = k_ss - ublas::prod(k_t, k);
    for(int i = 0; i < test_size; i++) ret.push_back(Result{y_mean(i, 0), y_cov(i, i)});
    return ret;
  }

private:
  std::vector<pii> X_train;
  std::vector<double> y_train;
  std::vector<double> noises;
  ublas::matrix<double> K; // kernel(x_train, x_train)
  ublas::matrix<double> K_inv;
  ublas::matrix<double> alpha; // pre-calculated matrix. used for prediction
  bool is_trained;
  double coef;
  double length;
  double y_train_mean;
  double y_train_std;

  // radius basis function kernel
  double rbf_kernel(const pii& x1, const pii& x2){
    return this->coef * std::exp(-(pow2(x1.first-x2.first) + pow2(x1.second-x2.second)) / (2 * pow2(this->length)));
  }
};

Grid make_init_lines() {
  int N = gl::N;
  int W = gl::W;
  int K = gl::K;
  int pi[W+K];
  int pj[W+K];
  for(int i = 0; i < W+K; i++) {
    pi[i] = i < W ? gl::wi[i] : gl::hi[i];
    pj[i] = i < W ? gl::wj[i] : gl::hj[i];
  }

  // calc distances
  std::vector<std::pair<ll, pii>> dists;
  for(int i = 0; i < W + K; i++) {
    for(int j = i + 1; j < W + K; j++) {
      long long dist = pow2(pi[i] - pi[j]) + pow2(pj[i] - pj[j]);
      dists.push_back({dist, pii{i, j}});
    }
  }
  std::sort(dists.begin(), dists.end());
  std::vector<std::vector<int>> G(W+K);
  std::vector<pii> edges;
  std::vector<int> boss(W+K, -1);

  // construct basic minimum spanning tree
  dsu uf(W+K);
  int cnt = 0;
  for(auto p: dists) {
    long long dist = p.first;
    int i = p.second.first;
    int j = p.second.second;
    if(uf.same(i, j) or (i < W and j < W)) continue;
    if(i < W and boss[j] != -1) continue;
    if(boss[i] != -1 and boss[j] != -1 and boss[i] != boss[j]) continue;
    G[i].push_back(j);
    G[j].push_back(i);
    uf.merge(i, j);
    cnt++;
    edges.push_back({i, j});
    for(int w = 0; w < W; w++) {
      for(int k = 0; k < K; k++) if(uf.same(w, k+W)) boss[W+k] = w;
    }
  }

  // basic mst summary
  int M = edges.size();
  cerr << M << endl;
  for(auto e: edges) std::cerr << e.first << " " << e.second << std::endl;
  
  // find precise mst
  long long L = 1e18;
  std::vector<int> ans(M);
  std::vector<std::vector<bool>> ans_grid;
  auto fill = [](int si, int sj, int ti, int tj, std::vector<std::vector<bool>>& v) {
    if(si == ti) for(int j = std::min(sj, tj); j <= std::max(sj, tj); j++) v[si][j] = true;
    else for(int i = std::min(si, ti); i <= std::max(si, ti); i++) v[i][sj] = true;
  };
  std::vector<std::vector<bool>> grid(N, std::vector<bool>(N, false));
  for(int i = 0; i < (1<<M); i++) {
    grid.clear();
    grid.resize(N, std::vector<bool>(N, false));
    long long l = 0;
    for(int j = 0; j < M; j++) {
      int s = edges[j].first;
      int t = edges[j].second;
      
      int psi, psj;
      int pti, ptj;
      if(s < W) psi = pi[s], psj = pj[s];
      else psi = pi[s], psj = pj[s];
      if(t < W) pti = pi[t], ptj = pj[t];
      else pti = pi[t], ptj = pj[t];
      // std::cerr << "("<<psx << ","<<psy<<"), ("<<ptx<<","<<pty<<")" << std::endl;
      if((1<<j)&i) {
        int pui = pti;
        int puj = psj;
        fill(psi, psj, pui, puj, grid);
        fill(pui, puj, pti, ptj, grid);
      } else {
        int pui = psi;
        int puj = ptj;
        fill(psi, psj, pui, puj, grid);
        fill(pui, puj, pti, ptj, grid);
      }
    }
    for(int r = 0; r < N; r++) for(int c = 0; c < N; c++) if(grid[r][c]==1) l++;
    if(l < L) {
      L = l;
      ans_grid = grid;
      for(int j = 0; j < M; j++) ans[j] = i&(1<<j) ? 1 : 0;
    }
  }
  std::cerr << "Shortest Line: " << L << std::endl;

  return ans_grid;
}

int main() {
  std::cin >> gl::N >> gl::W >> gl::K >> gl::C;
  for(int i = 0; i < gl::N; i++) {
    for(int j = 0; j < gl::N; j++) {
      std::cin >> gl::S[i][j];
    }
  }
  for(int i = 0; i < gl::W; i++) std::cin >> gl::wi[i] >> gl::wj[i];
  for(int j = 0; j < gl::K; j++) std::cin >> gl::hi[j] >> gl::hj[j];
  
  Grid init_lines = make_init_lines();

  int path_length = 0;
  for(int i = 0; i < gl::N; i++) {
    for(int j = 0; j < gl::N; j++) {
      if(init_lines[i][j]) {
        for(int w = 0; w < gl::W; w++) if(gl::wi[w] == i and gl::wj[w] == j) continue;
        path_length++;
      }
    }
  }
  Vec range, targets;
  for(int i = 0; i < path_length; i++) range.push_back(i);
  std::sample(
    range.begin(), 
    range.end(), 
    std::back_inserter(targets), 
    path_length / 3, 
    mt
  );
  
  std::vector<pii> X;
  std::vector<double> y;
  std::vector<double> noises;

  auto it = targets.begin();
  int cnt = 0;
  std::cerr << "FIND CANDIDATE" << std::endl;
  for(int i = 0; i < gl::N; i++) {
    for(int j = 0; j < gl::N; j++) {
      if(init_lines[i][j]) {
        if(*it == cnt) {
          int power = std::max(gl::C * 5, 10);
          int res = 0;
          int total_power = (gl::S[i][j] + power - 1) / power * power;

          X.push_back({i, j});
          y.push_back((double)total_power - (double)power / 2.);
          noises.push_back(std::sqrt(pow2((double)power) / 12.));
          std::cerr << i << " " << j << " " << y.back() << " " << noises.back() << std::endl;
          it++;
        }
        cnt++;
      }
      if(it == targets.end()) break;
    }
    if(it == targets.end()) break;
  }

  GaussianProcess gp(param::COEF, param::LENGTH);
  std::vector<std::vector<std::pair<double, double>>> preds(gl::N, std::vector<std::pair<double, double>>(gl::N));
  gp.fit(X, y, noises);
  for(int i = 0; i < gl::N; i++) {
    for(int j = 0; j < gl::N; j++) {
      std::pair<double, double> res = gp.predict({i, j});
      std::cerr << res.first << " ";
      preds[i][j] = res;
    }
    std::cerr << std::endl;
  }

  for(int i = 0; i < gl::N; i++) {
    for(int j = 0; j < gl::N; j++) {
      std::cout << preds[i][j].first - gl::S[i][j] << " ";
    }
    std::cout << std::endl;
  }



  std::cerr << get_time() << std::endl;
  return 0;
}