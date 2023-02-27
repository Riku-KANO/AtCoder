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
#include <memory>
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

bool operator==(const pii& lhs, const pii& rhs) {
  return (lhs.first == rhs.first && lhs.second == rhs.second);
}

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
  int starts[8] = {13, 34, 41, 58, 65, 116, 165, 205};
  int coefs[8] = {
    1.0292408601924954,
    1.0220879491309989,
    1.0404330624789304,
    1.0525996160890472,
    1.1033116911939742,
    1.1139859957590141,
    1.1678180810131298,
    1.2544750301323389
  };
}

struct IOServer {
  int S[204][204];
  int total_cost = 0;
  void init() {
    for(int i = 0; i < gl::N; i++) {
      for(int j = 0; j < gl::N; j++) {
        std::cin >> S[i][j];
      }
    }
  }

  int read_output(int i, int j, int power) {
    int res;
    if(power > 5000 or power < 1) {
      std::cerr << "Invalid power" << std::endl;
      res = -1;
    } else if(i < 0 or i >= gl::N or j < 0 or j >= gl::N) {
      std::cerr << "Invalid position" << std::endl;
      res = -1;
    } else if(this->S[i][j] <= 0) {
      std::cerr << "The bedrock is already broken" << std::endl;
      res = -1;
    } else if(S[i][j] - power <= 0) {
      S[i][j] -= power;
      if(is_connected()) res = 2;
      else res = 1;
    } else {
      S[i][j] -= power;
      res = 0;
    }
    total_cost += power + gl::C;
    return res;
  }

  bool is_connected() {
    std::vector<std::vector<bool>> visited(gl::N, std::vector<bool>(gl::N, false));

    // bfs
    for(int i = 0; i < gl::W; i++) {
      if(S[gl::wi[i]][gl::wj[i]] > 0) continue;
      std::queue<pii> q;
      q.push({gl::wi[i], gl::wj[i]});
      while(!q.empty()) {
        auto p = q.front();
        int r = p.first;
        int c = p.second;
        if(visited[r][c]) continue;
        visited[r][c] = true;
        for(int d = 0; d < 4; d++) {
          int ni = r + di[d];
          int nj = c + dj[d];
          if(is_out(ni, nj, gl::N, gl::N)) continue;
          if(S[ni][nj] <= 0) q.push({ni, nj});
        }
      }
    }

    for(int i = 0; i < gl::W; i++) if(!visited[gl::wi[i]][gl::wj[i]]) return false;
    for(int i = 0; i < gl::K; i++) if(!visited[gl::hi[i]][gl::hj[i]]) return false;
    return true;
  }
};


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

// 
double acquisition_function() {
  double ret;

  return ret;
}

// # memo
// 水源と家の最短経路が必ずしも良いとは限らない(その間に固い岩盤がある可能性がある)。
// 逆行列計算(200*200): C++(603ms), python(110ms)...C++なのに遅い？？
//
// C=1: a=38, coef=1.071071836416865
// C=2: a=42, coef=1.1103374315829098
// C=4: a=48, coef=1.1581850602369728
// C=8: a=49, coef=1.2295262864675027
// C=16: a=87, coef=1.2309734139372803
// C=32: a=109, coef=1.4340920873990537
// C=64: a=152, coef=1.4505429558219907
// C=128: a=234, coef=1.7479564275184625
// 最小全域木のノード数は多くても1000個ほどになりそう。

// # idea
// 最小全域木っぽくしたほうがいい（家を経由する。無駄な掘削を減らす）
// 
// Cが大きいときはなるべく少ない回数で掘削を行った方がいい。
// 逆にCが小さいときは何回でも地盤調査をした方がいい。
// 観測値と実測値に差がある(観測値は必ず正の誤差をもつ)。
// 観測値のエラーは掘削の回数に依存する(細かい回数掘削したら岩盤の正確な硬さが分かるが一発で掘削しきると正確な硬さは不明)。
// 正確な硬さを知る必要性はあるのか？
// ガウス過程回帰は計算コストが高いので測定できるポイントは限られる。
// 掘削しながらガウス過程回帰もコストが高すぎる。
// -> 観測点を近いものだけ選ぶようにしたら計算コストがへらせそう！
// -> オリジナルの獲得関数を定義したらよさそう。重点的に考える。
// 不確実性探索もよさそう。
// 実行時間が短すぎるので探索する余地があまりなさそう。

// ## idea1
// 1. MSTで水源の道の決定（これをベースとして変更を加えたりする。）
// 2. 水源と家付近の掘削。岩盤の固さとその誤差を測定
// 3. 
// 結果：pythonで試作品を動かしてみたがガウス過程回帰の精度もあまり高くはない。
// 岩盤が固いところを避けていく方がハイスコアに結び付くかもしれない。

// ## idea2
// 1. 200*200の区画を5*5の小さい区画に分けてそこでMSTをつくる（マンハッタン距離が400/(W+K)以上であることが保証されている）。
// つまり40*40のグリッド上でMSTを構成する(この時辺の総数3120)。
// 2. MST上のノードのareaを一回ずつ叩く->岩盤の固さの更新
// 3. MSTを再構成
// 4. MSTが確定するまで2,3を繰り返し行う。
// 5. 

struct Area {
  double stiff_exp;
  int total_power;
  int cnt;
  bool has_point;
  bool done;
  Area() {
    stiff_exp = 2505;
    total_power = 0;
    cnt = 0;
    has_point = false;
    done = false;
  }
  Area(double stiff_exp_, int total_power_, int cnt_, bool has_point_, bool done_) {
    stiff_exp = stiff_exp_;
    total_power = total_power_;
    cnt = cnt_;
    has_point = has_point_;
    done = done_;
  }
};


class Solver {
public:
  typedef std::vector<std::vector<bool>> Grid;
  #ifdef _LOCAL
  IOServer judge;
  #endif 
  std::vector<std::vector<double>> mu;
  std::vector<std::vector<double>> mu_pred;
  std::vector<std::vector<double>> noises;
  std::vector<std::vector<double>> noises_pred;
  std::vector<double> powers;
  bool done;

  Solver() {
    this->done = false;
  } // constructor

  void init() {
    cerr << "INITIALIZATION STARTED" << endl;
    std::cin >> gl::N >> gl::W >> gl::K >> gl::C;
    #ifdef _LOCAL
    judge.init();
    #endif
    for(int i = 0; i < gl::W; i++) std::cin >> gl::wi[i] >> gl::wj[i];
    for(int i = 0; i < gl::K; i++) std::cin >> gl::hi[i] >> gl::hj[i];
    this->done = false;
    this->mu.resize(gl::N, std::vector<double>(gl::N, (-1.)));
    this->mu_pred.resize(gl::N, std::vector<double>(gl::N, ((double)gl::MAXP+gl::MINP)/2));
    this->noises.resize(gl::N, std::vector<double>(gl::N, -1.));
    this->noises_pred.resize(gl::N, std::vector<double>(gl::N, std::sqrt(pow2((double)gl::MAXP-gl::MINP)/12)));
    for(int i = 0; i < 8; i++) {
      if((gl::C >> i) & 1) {
        double power = param::starts[i];
        double coef = param::coefs[i];
        double total_power = 0.0;
        while(total_power < 5000) {
          this->powers.push_back(power);
          total_power += power;
          power *= coef;
        }
      }
    }
    
    cerr << "INITIALIZATION ENDED" << endl;
  }


  // main solve
  // idea2
  void solve() {
    // 5*5の小さい区画に分割
    int small_grid_size = 5;
    std::vector<std::vector<Area>> meta_grid(gl::N / small_grid_size, std::vector<Area>(gl::N / small_grid_size));
    // meta area of sources
    for(int i = 0; i < gl::W; i++) {
      int mi = gl::wi[i] / small_grid_size;
      int mj = gl::wj[i] / small_grid_size;
      auto t = this->excavation_till_break(gl::wi[i], gl::wj[i], this->powers); 
      int stiff = std::get<0>(t);
      int total_power = std::get<1>(t);
      int cnt = std::get<2>(t);
      meta_grid[mi][mj] = Area(stiff, total_power, cnt, true, true);
    }
    // meta area of houses
    for(int i = 0; i < gl::K; i++) {
      int mi = gl::hi[i] / small_grid_size;
      int mj = gl::hj[i] / small_grid_size;
      auto t = this->excavation_till_break(gl::hi[i], gl::hj[i], this->powers);
      double stiff = std::get<0>(t);
      int total_power = std::get<1>(t);
      int cnt = std::get<1>(2);
      meta_grid[mi][mj] = Area(stiff, total_power, cnt, true, true);
    }

    // find meta grid
    Grid cur_mst = this->construct_mst(meta_grid);
    Grid next_mst = cur_mst;
    while(true) {
      this->hit_on_mst(meta_grid);
      next_mst = this->construct_mst(meta_grid);
      if(this->same_mst(cur_mst, next_mst)) break;
      cur_mst = next_mst;
    }

    
    
  } // solve

  // return summaries of solve
  void summary() {

  }


  // debug gaussian process
  void check_gp() {
    Grid init_lines = this->make_init_lines();
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
            int total_power = 0;
            while(res == 0) {
              res = this->excavation(i, j, power);
              total_power += power;
            }

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
    gp.fit(X, y, noises);
    for(int i = 0; i < gl::N; i++) {
      for(int j = 0; j < gl::N; j++) {
        std::pair<double, double> res = gp.predict({i, j});
        std::cerr << res.first << " ";
      }
      std::cerr << std::endl;
    }
    std::cerr << get_time() << std::endl;
  }

private:

  Grid construct_mst(const std::vector<std::vector<Area>> &grid) {
    int W = gl::W;
    int K = gl::K;
    int sz = grid.size();
    Grid ret(sz, std::vector<bool>(sz, false));
    std::vector<std::vector<double>> dists(gl::W+gl::K, std::vector<double>(gl::W+gl::K));
    for(int i = 0; i < gl::W; i++) {
      pii start = {gl::wi[i] / 5, gl::wj[i] / 5};
      std::vector<std::vector<double>> dist = this->dijkstra(start, grid);
      for(int j = 0; j < gl::W+gl::K; j++) {
        int ti = j < gl::W ? gl::wi[j] / 5 : gl::hi[j-gl::W] / 5;
        int tj = j < gl::W ? gl::wj[j] / 5 : gl::hj[j-gl::W] / 5;
        dists[i][j] = dists[j][i] = dist[ti][tj];
      }
    }
    for(int i = 0; i < gl::K; i++) {
      pii start = {gl::hi[i] / 5, gl::hj[i] / 5};
      std::vector<std::vector<double>> dist = this->dijkstra(start, grid);
      for(int j = 0; j < W+K; j++) {
        int ti = j < W ? gl::wi[j] / 5 : gl::hi[j-gl::W] / 5;
        int tj = j < W ? gl::wj[j] / 5 : gl::hj[j-gl::W] / 5;
        dists[i][j] = dists[j][i] = dist[ti][tj];
      }
    }

    std::vector<pair<double, pii>> dist_list;
    for(int i = 0; i < W+K; i++) {
      for(int j = i + 1; j < W+K; j++) {
        dist_list.push_back(dist[i][j], {i, j});
      }
    }

    std::sort(dist_list.begin(), dist_list.end());
    std::vector<pii> edges;
    std::vector<int> boss(W+K, -1);

    // construct basic minimum spanning tree
    dsu uf(W+K);
    int cnt = 0;
    for(auto p: dist_list) {
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

    for(auto edge: edges) {
      pii start = edge.first < W ? gl::wi[edge.first] : gl::hi[edge.first];
      pii goal = edge.second < W ? gl::wj[edge.second] : gl::hj[edge.second];
      Grid path = dijkstra(start, goal, grid);
      for(int i = 0; i < sz; i++) {
        for(int j = 0; j < sz; j++) {
          if(path[i][j])continue;
          ret[i][j] = true;
        }
      }
    }

    return ret;
  }

  void hit_on_mst(const Grid &mst, std::vector<std::vector<Area>> &grid) {
    int sz = mst.size();
    for(int i = 0; i < sz; i++) {
      for(int j = 0; j < sz; j++) {
        if(!mst[i][j] or grid[i][j].done) continue;
        int power = (int)this->powers[grid[i][j].cnt];
        power = std::min(power, 5000-grid[i][j].total_power);
        int res = this->excavation(i * 5 + 2, j * 5 + 2, power);
        grid[i][j].cnt++;
        grid[i][j].total_power += power;
        if(res == 1) {
          grid[i][j].done = true;
          grid[i][j].stiff_exp = (double)grid[i][j].total_power - (double)power/2;
        } else if(res == 0) {
          grid[i][j].stiff_exp = (double)(5000 + grid[i][j].total_power) / 2;
        } else if(res == -1) {
          std::cerr << "Invalid output from hit_on_mst function" << std::endl;
          std::exit(1);
        }
      }
    }
  }

  bool same_mst(const Grid &a, const Grid &b) {
    int sz = a.size();
    for(int i = 0; i < sz; i++) for(int j = 0; j < sz; j++) if(a[i][j] != b[i][j]) return false;
    return true;
  }

  std::vector<std::vector<double>> dijkstra(const pii& start, const std::vector<std::vector<Area>> &grid) {
    typedef std::pair<double, pair<pii, pii>> Edge;
    int sz = grid.size();
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, 1e18));
    std::priority_queue<Edge> pq;
    dist[start.first][start.second] = 0.0;
    pq.push({0., {{-1, -1}, start}});
    while(!pq.empty()) {
      auto edge = pq.top(); pq.pop();
      double cost = edge.first;
      int ci = edge.second.second.first;
      int cj = edge.second.second.second;
      if(cost > dist[ci][cj] + EPS) continue;
      for(int d = 0; d < 4; d++) {
        int ni = ci + di[d];
        int nj = cj + dj[d];
        double w = (grid[ci][cj].stiff_exp + grid[ni][nj].stiff_exp) / 2;
        if(is_out(ni, j, sz, sz)) continue;
        if(dist[ni][nj] > dist[ci][cj] + w + EPS) {
          dist[ni][nj] = dist[ci][cj] + w;
          pq.push(dist[ni][nj], {{ci, cj}, {ni, nj}});
        }
      }
    }
    return dist;
  }

  // calc Grid path
  Grid dijkstra(const pii& s, const pii& t, const std::vector<std::vector<Area>& gtid) {
    typedef std::pair<double, pair<pii, pii>> Edge;
    int sz = grid.size();
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, std::vector<double>(sz, 1e18)));
    std::vector<std::vector<pii>> pre(sz, std::vector<pii>(sz, {-1, -1}));
    Grid ret(sz, std::vector<bool>(sz, false));
    std::priority_queue<Edge> pq;
    dist[s.first][s.second] = 0.0;
    pre[s.first][s.second] = s;
    pq.push({0, {{-1, -1}, s}});
    while(!pq.empty()) {
      auto edge = pq.top(); pq.pop();
      double cost = edge.first;
      int ci = edge.second.second.first;
      int cj = edge.second.second.second;
      if(cost > dist[ci][cj] + EPS or cost > dist[t.first][t.second] + EPS) continue;
      for(int d = 0; d < 4; d++) {
        int ni = ci + di[d];
        int nj = cj + dj[d];
        double w = (grid[ci][cj].stiff_exp + grid[ni][nj].stiff_exp) / 2;
        if(is_out(ni, j, sz, sz)) continue;
        if(dist[ni][nj] > dist[ci][cj] + w + EPS) {
          dist[ni][nj] = dist[ci][cj] + w;
          pq.push(dist[ni][nj], {{ci, cj}, {ni, nj}});
          pre[ni][nj] = {ci, cj};
        }
      }
    }

    // calc Grid
    pii cur = t;
    while(true) {
      ret[cur.first][cur.second] = true;
      if(cur == s) break;
      cur = pre[cur.first][cur.second];
    }

    return ret;
  }



  // return tiredness and total power
  std::tuple<int, int> excavation_till_break(int x, int y, int power) {
    assert(0 <= power and power <= 5000);
    int total_power = 0;
    int tiredness = 0;
    bool broken = false;
    while(!broken) {
      std::cout << x << " " << y << " " << power << std::endl;
      total_power += power;
      tiredness += (power + gl::C);
      int res;
      std::cin >> res;
      if(res == 0) continue;
      else if(res == 1) break;
      else if(res == 2) {
        std::cerr << "All done!" << std::endl;
        this->done = true;
        break;
      } else if(res == -1) {
        std::cerr << "Invarid Output!" << std::endl;
        std::exit(1);
      }
    }
    return {tiredness, total_power};
  }

  // powers from parameters
  std::tuole<double, int, int> excavation_till_break(int x, int y, cosnt std::vector<double>&powers) {
    int total_power = 0;
    int cnt = 0;
    bool broken = false;
    auto it = powers.begin();
    int power = (int)*it
    while(!broken) {
      power = std::min(power, total_power - power);
      int res;
      #ifdef _LOCAL
      res = this->judge.read_output(x, y, power);
      #else
      std::cout << x << " " << y << " " << power << std::endl;
      std::cin >> res;
      #endif
      total_power += power;
      cnt++;
      if(res == 0) {
        it++;
        power = (int)*it
        continue;
      }
      else if(res == 1) break;
      else if(res == 2) {
        std::cerr << "All done!" << std::endl;
        this->done = true;
        break;
      } else if(res == -1) {
        std::cerr << "Invarid Output!" << std::endl;
        std::exit(1);
      }
    }
    return {(double)total_power - (double)power / 2, total_power, cnt};
  }

  // normal excavation
  int excavation(int x, int y, int power) {
    assert(0 <= power and power <= 5000);
    int res;
    #ifdef _LOCAL
      res = this->judge.read_output(x, y, power)
    #else
      std::cout << x << " " << y << " " << power << std::endl;
      std::cin >> res;
    #endif
    return res;
  }

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

};


int main(int argc, char* argv[]) {
  start_time = clock();
  Solver solver;

  #ifdef _OPTUNA
  for(int i = 1; i < argc; i++) {
    if(argv[i] == "--seed") {
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return 1;
      } else {
        int seed = std::stoi(argv[i+1]);
        mt = std::mt19937(seed);
      }
      i += 2;
    } else {
      std::cerr << "Invalid argument" << std::endl;
      return 1;
    }
  }
  #endif

  #ifdef _DEBUG
  solver.init();
  solver.check_gp();
  return 0;
  #endif
  solver.init();
  solver.solve();
  solver.summary();
  return 0;
}