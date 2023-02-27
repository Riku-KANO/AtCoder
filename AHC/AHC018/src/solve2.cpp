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
#include <cstdlib>

// atcoder
#include <atcoder/dsu>

// // boost
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/lu.hpp>
// #include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/io.hpp>

#define pii std::pair<int,int>
#define Vec std::vector<int>
#define rep(i, s, n) for(int i = s; i < (n); i++)

using namespace std;
using namespace atcoder;
// namespace ublas = boost::numeric::ublas;
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

bool is_out(int ni, int nj, int N) {
  return (nj < 0 or ni < 0 or ni >= N or nj >= N);
}

bool is_out(int ni, int nj, int H, int W) {
  return (nj < 0 or ni < 0 or ni >= H or nj >= W);
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
  const int MAXP = 5000;
  const int MINP = 10;
}

namespace param {
  const double TL = 5.0;
  const double TL_90 = TL * 0.90;
  const double TL_95 = TL * 0.95;
  int seed;
  double COEF = 1.0;
  double LENGTH = 7.0;
  int SMALL_GRID = 10;
  int INIT_STIF_EXP = 200;
  int starts[8] = {13, 34, 41, 58, 65, 116, 165, 205};
  double coefs[8] = {
    1.0292408601924954,
    1.0220879491309989,
    1.0404330624789304,
    1.0525996160890472,
    1.1033116911939742,
    1.1139859957590141,
    1.1678180810131298,
    1.2544750301323389
  };

  int pre_defined_C = -1;
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
      if(this->is_connected()) res = 2;
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
        auto p = q.front(); q.pop();
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


// // calc invert matrix by LU decomposition
// template<typename T> ublas::matrix<T> calc_invmat(ublas::matrix<T> A){
//   ublas::matrix<T> B = ublas::identity_matrix<T>(A.size1());
//   ublas::permutation_matrix<std::size_t> pm(A.size1());
//   int res = ublas::lu_factorize(A, pm);
//   if(res != 0) {
//     std::cerr << "calc inverse Failed" << std::endl;
//     return B;
//   }
//   ublas::lu_substitute(A, pm, B);
//   return B;
// }

// class GaussianProcess {
// public:
//   GaussianProcess(double coef = 1.0, double length = 5.0) {
//     this->coef = coef;
//     this->length = length;
//     this->is_trained = false;
//   }

//   void fit(const std::vector<std::pair<int,int>>& pos, const std::vector<double>& y, const std::vector<double>& noises) {
//     std::cerr << "start fit" << std::endl;
//     size_t ker_size = pos.size();
//     this->X_train = pos;
//     this->y_train = y;
//     this->noises = noises;
//     ublas::matrix<double> y_vec(ker_size, 1);
//     ublas::matrix<double> k_train(ker_size, ker_size);
//     double y_train_mean_ = 0.0;
//     double y_train_std_ = 0.0;

//     rep(i, 0, ker_size) {
//       y_vec(i, 0) = y[i];
//       y_train_mean_ += y[i];
//     }
//     y_train_mean_ /= (int)ker_size;
//     this->y_train_mean = y_train_mean_;

    
//     rep(i, 0, ker_size) y_train_std_ += pow2(y[i]-y_train_mean_);
    
//     y_train_std_ = std::sqrt(y_train_std_ / ker_size);
//     this->y_train_std = y_train_std_;
//     rep(i, 0, ker_size) {
//       rep(j, 0, ker_size) {
//         k_train(i, j) = this->rbf_kernel(pos[i], pos[j]);
//       }
//     }

//     this->K = k_train;
//     rep(i, 0, ker_size) k_train(i, i) += pow2(noises[i]);
//     ublas::matrix<double> k_train_inv = calc_invmat(k_train);
//     this->K_inv = k_train_inv;
//     ublas::matrix<double> alp = ublas::prod(k_train_inv, y_vec);
//     this->alpha = alp;
//     this->is_trained = true;
//     // std::cerr << "alpha: " << this->alpha << std::endl;
//     // std::cerr << "K: " << this->K << std::endl;
//     // std::cerr << "K inv: " << this->K_inv << std::endl;
//     std::cerr << "y train mean: " << this->y_train_mean << std::endl;
//     std::cerr << "y train std: " << this->y_train_std << std::endl;
//   }

//   std::pair<double, double> predict(const pii& x){
//     if(! this->is_trained) {
//       std::cerr << "This instance is not fitted. please fit before predict" << std::endl;
//       return {-1, -1};
//     }
//     size_t ker_size = this->X_train.size();
//     std::pair<double, double> ret;
//     ublas::matrix<double> k_s(ker_size, 1);
//     double k_ss;
//     rep(i, 0, ker_size) k_s(i, 0) = this->rbf_kernel(this->X_train[i], x);
//     k_ss = this->rbf_kernel(x, x);
    
//     ublas::matrix<double> y_tmp = ublas::prod(ublas::trans(k_s), this->alpha);
    
//     double y_mean = y_tmp(0,0);
//     // std::cerr << y_mean << std::endl;
//     ublas::matrix<double> k_t = ublas::trans(k_s);
//     ublas::matrix<double> k_2 = ublas::prod(this->K_inv, k_s);
//     // std::cerr << k_t << std::endl;
//     // std::cerr << k_2 << std::endl;
//     double y_var = k_ss - ublas::prod(k_t , k_2)(0, 0);
//     ret = {y_mean, y_var};
//     return ret;
//   }

//   std::vector<std::pair<double, double>> predict(const std::vector<pii>& x){
//     if(! this->is_trained) {
//       std::cerr << "This instance is not fitted. please fit before predict" << std::endl;
//       std::vector<pair<double,double>> ret;
//       return ret;
//     }
//     typedef std::pair<double, double> Result;
//     size_t ker_size = this->X_train.size();
//     size_t test_size = x.size();
//     std::vector<Result> ret;
//     ublas::matrix<double> k_s(ker_size, test_size);
//     ublas::matrix<double> k_ss(test_size, test_size);
//     rep(i, 0, ker_size) {
//       rep(j, 0, test_size) {
//         k_s(i, j) = this->rbf_kernel(this->X_train[i], x[j]);
//       }
//     }
//     rep(i, 0, test_size) {
//       rep(j, 0, test_size) {
//         k_ss(i, j) = this->rbf_kernel(x[i], x[j]);
//       }
//     }
//     ublas::matrix<double> y_mean = ublas::prod(this->alpha, k_s);
//     ublas::matrix<double> k_t = ublas::trans(k_s);
//     ublas::matrix<double> k = ublas::prod(this->K_inv, k_s);
//     ublas::matrix<double> y_cov = k_ss - ublas::prod(k_t, k);
//     for(int i = 0; i < test_size; i++) ret.push_back(Result{y_mean(i, 0), y_cov(i, i)});
//     return ret;
//   }

// private:
//   std::vector<pii> X_train;
//   std::vector<double> y_train;
//   std::vector<double> noises;
//   ublas::matrix<double> K; // kernel(x_train, x_train)
//   ublas::matrix<double> K_inv;
//   ublas::matrix<double> alpha; // pre-calculated matrix. used for prediction
//   bool is_trained;
//   double coef;
//   double length;
//   double y_train_mean;
//   double y_train_std;

//   // radius basis function kernel
//   double rbf_kernel(const pii& x1, const pii& x2){
//     return this->coef * std::exp(-(pow2(x1.first-x2.first) + pow2(x1.second-x2.second)) / (2 * pow2(this->length)));
//   }
// };

// // 
// double acquisition_function() {
//   double ret;

//   return ret;
// }

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
// 探索と活用のうまい塩梅を考えるべきなのだろうが、難しい。
// シミュレーションを行って、木のコストの期待改善率など計算出来たらなぁと思うが5秒では無理。

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
// 一回探索した後の硬さの推定値をtotal_powerの数倍にしたらseed0でスコアが半分に。良い探索ができたみたい。


// ## idea1
// 1. MSTで水源の道の決定（これをベースとして変更を加えたりする。）
// 2. 水源と家付近の掘削。岩盤の固さとその誤差を測定
// 3. 
// 結果：pythonで試作品を動かしてみたがガウス過程回帰の精度もあまり高くはない。
// 岩盤が固いところを避けていく方がハイスコアに結び付くかもしれない。
// 実行時間制限が長ければこの方法がベストソリューションな気もする。

// ## idea2
// 1. 200*200の区画を5*5の小さい区画に分けてそこでMSTをつくる（マンハッタン距離が400/(W+K)以上であることが保証されている）。
// つまり40*40のグリッド上でMSTを構成する(この時辺の総数3120)。 
// (5 * 5の時のグリッドを可視化したら意外に小さかったので8*8でもよさそう。このとき、25*25)
// 2. MST上のノードのareaを一回ずつ叩く->岩盤の固さの更新
// 3. MSTを再構成
// 4. MSTが確定するまで2,3を繰り返し行う。
// 5. 大雑把なMSTをもとに辺を張る。
// xxxxxxxxxooooxxxxxxxxxxxx
// xxooxxxxxoxxxxxxxxxxxxxxx
// xxxoooooooxxxxxxxxxxxxxxx
// xxxxxxoxxxxxxxxxxxxxxxxxx
// xxxxxxoxxxxxxoxxxxxxxxxxx
// xxxxxxoxxxxxxoxxxxxxxxxxx
// xxxxooooooooooxxxxxxxxxxx
// xxxxxxxxxxxxxxxxxxxxxxxxx

// 2/22 掘った点の間を補間するのもありかもしれない。ポテンシャルダイクストラもあり。

// 課題: 
// グリッドの木をよりMSTに近づけるべき。
// おそらくC=1, 2, 4のときでは大差をつけられていない。()
// C=128の時は乱数次第で10~20%ぶれる。
// 探索の回数を最小化

// ここまでのまとめ: Cが小さいときは地盤の探索をするようにした方がハイスコアに結び付く。
// (seed=0で硬いところを迂回したら半分以上スコアが下がった)
// Cが大きいときはできるだけ少ない回数で探索をし、経路を見つける。
// パラメータの最適化
// 岩盤の硬さを推定して効率的に破壊


struct Area {
  pii pos;
  double stiff_exp;
  int total_power;
  int cnt;
  bool has_point;
  bool done;
  Area() {
    stiff_exp = param::INIT_STIF_EXP;
    total_power = 0;
    cnt = 0;
    has_point = false;
    done = false;
  }
  Area(pii pos_, double stiff_exp_, int total_power_, int cnt_, bool has_point_, bool done_) {
    pos = pos_;
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
  std::vector<double> powers;
  int num_search = 0;
  int search_cost = 0;
  bool done;

  Solver() {
    this->done = false;
  } // constructor

  void init() {
    std::cerr << "INITIALIZATION STARTED" << std::endl;
    std::cin >> gl::N >> gl::W >> gl::K >> gl::C;
    if(param::pre_defined_C != -1) gl::C = param::pre_defined_C; 
    #ifdef _LOCAL
    judge.init();
    #endif
    for(int i = 0; i < gl::W; i++) std::cin >> gl::wi[i] >> gl::wj[i];
    for(int i = 0; i < gl::K; i++) std::cin >> gl::hi[i] >> gl::hj[i];
    this->done = false;
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
    std::cerr << "POWERS: " << std::endl;
    for(double power: this->powers) std::cerr << power << " ";
    std::cerr << std::endl;

    std::cerr << "INITIALIZATION ENDED" << std::endl;
  }


  // main solve
  // idea2
  void solve() {
    // S_GRID_SIZE * S_GRID_SIZEの小さい区画に分割
    const int S_GRID_SIZE = param::SMALL_GRID;
    assert(gl::N % S_GRID_SIZE == 0);

    std::vector<std::vector<Area>> meta_grid(gl::N / S_GRID_SIZE, std::vector<Area>(gl::N / S_GRID_SIZE));
    // meta area of sources
    for(int i = 0; i < gl::W; i++) {
      int mi = gl::wi[i] / S_GRID_SIZE;
      int mj = gl::wj[i] / S_GRID_SIZE;
      auto t = this->excavation_till_break(gl::wi[i], gl::wj[i]); 
      double stiff = std::get<0>(t);
      int total_power = std::get<1>(t);
      this->search_cost += total_power;
      int cnt = std::get<2>(t);
      meta_grid[mi][mj] = Area(pii{gl::wi[i], gl::wj[i]}, stiff, total_power, cnt, true, true);
    }
    // meta area of houses
    for(int i = 0; i < gl::K; i++) {
      int mi = gl::hi[i] / S_GRID_SIZE;
      int mj = gl::hj[i] / S_GRID_SIZE;
      auto t = this->excavation_till_break(gl::hi[i], gl::hj[i]);
      double stiff = std::get<0>(t);
      int total_power = std::get<1>(t);
      this->search_cost += total_power;
      int cnt = std::get<2>(t);
      meta_grid[mi][mj] = Area(pii{gl::hi[i], gl::hj[i]}, stiff, total_power, cnt, true, true);
    }

    std::cerr << "FINISHED FIRST EXCAVATION" << std::endl;


    // find meta grid
    Grid cur_meta_mst = std::get<0>(this->construct_mst(meta_grid));
    Grid next_meta_mst = cur_meta_mst;
    std::vector<pii> mst_edges;
    std::cerr << "START FINDING meta MST" << std::endl;
    int maxIter = this->powers.size();
    while(maxIter) {
      std::cerr << "OK, " << maxIter << std::endl;
      this->hit_on_mst(cur_meta_mst, meta_grid, false, true);
      auto t = this->construct_mst(meta_grid);
      next_meta_mst = std::get<0>(t);
      mst_edges = std::get<1>(t);
      if(this->same_mst(cur_meta_mst, next_meta_mst)) {
        maxIter--;
      } else {
        maxIter = this->powers.size();
        cur_meta_mst = next_meta_mst;
      }
    }

    std::cerr << "FOUND meta MST" << std::endl;


    std::cerr << "START FINDING ORIGINAL GRID TREE" << std::endl;
    // construct original grid info
    std::vector<std::vector<Area>> original_grid(gl::N, std::vector<Area>(gl::N));
    for(int i = 0; i < gl::N / S_GRID_SIZE; i++) {
      for(int j = 0; j < gl::N / S_GRID_SIZE; j++) {
        if(meta_grid[i][j].has_point) {
          pii pos = meta_grid[i][j].pos;
          original_grid[pos.first][pos.second] = meta_grid[i][j];
        } else {
          int ii = i * S_GRID_SIZE + S_GRID_SIZE / 2;
          int jj = j * S_GRID_SIZE + S_GRID_SIZE / 2;
          original_grid[ii][jj] = meta_grid[i][j];
        }

        for(int ii = S_GRID_SIZE * i; ii < (i+1)*S_GRID_SIZE; ii++) {
          for(int jj = S_GRID_SIZE * j; jj < (j + 1) * S_GRID_SIZE; jj++) {
            if(meta_grid[i][j].stiff_exp - EPS< param::INIT_STIF_EXP and meta_grid[i][j].stiff_exp + EPS > param::INIT_STIF_EXP) {
              original_grid[ii][jj].stiff_exp = 4000.;
            } else {
              if(meta_grid[i][j].done) original_grid[ii][jj].stiff_exp = meta_grid[i][j].stiff_exp * 1.1;    
              else original_grid[ii][jj].stiff_exp = std::min(meta_grid[i][j].stiff_exp * 2, 5000.);
            }
          }
        }
      }
    }

    // define valid area from meta-grid mst result
    std::vector<std::vector<bool>> valid_area(gl::N, std::vector<bool>(gl::N, false));
    for(int i = 0; i < gl::N / S_GRID_SIZE; i++) {
      for(int j = 0; j < gl::N / S_GRID_SIZE; j++) {
        if(cur_meta_mst[i][j]) {
          for(int y = i * S_GRID_SIZE; y < (i+1)*S_GRID_SIZE; y++) {
            for(int x = j * S_GRID_SIZE; x < (j + 1) * S_GRID_SIZE; x++) {
              valid_area[y][x] = true;
            }
          }
        }

      }
    }
  
    // 
    Grid cur_mst = this->construct_mst(original_grid, mst_edges, valid_area, true);
    Grid next_mst = cur_mst;
    // while(true) {
    //   bool all_done = this->hit_on_mst(cur_mst, original_grid, true);
    //   next_mst = this->construct_mst(original_grid, mvalid_area, true);
    //   if(this->same_mst(cur_mst, next_mst)) break;
    //   cur_mst = next_mst;
    // }

    while(!this->hit_on_mst(cur_mst, original_grid, true, false));
    
    
  } // solve

  // return summaries of solve
  void summary() {
    std::cerr << "\n### SUMMARY ####################" << std::endl;
    #ifdef _LOCAL
    double ratio = (double) (this->search_cost + this->num_search * gl::C) / this->judge.total_cost;
    std::cerr << "SCORE: " << this->judge.total_cost << std::endl;
    std::cerr << "NUMBER OF SEARCH: " << this->num_search << std::endl;
    std::cerr << "SEARCH COST: " << this->search_cost << std::endl;

    std::cerr << "SEARCH RATIO: " << std::fixed << std::setprecision(4) << ratio << std::endl;
    #endif
    std::cerr << "ELAPSED TIME: " << get_time() << std::endl;
    std::cerr << "################################\n" << std::endl;
  }

private:

  std::tuple<Grid, std::vector<pii>> construct_mst(const std::vector<std::vector<Area>> &grid) {
    int W = gl::W;
    int K = gl::K;
    int S_GRID_SIZE = param::SMALL_GRID;
    int sz = grid.size();
    Grid ret(sz, std::vector<bool>(sz, false));
    std::vector<std::vector<double>> dists(W+K, std::vector<double>(W+K));
    for(int i = 0; i < W; i++) {
      pii start = {gl::wi[i] / S_GRID_SIZE, gl::wj[i] / S_GRID_SIZE};
      std::vector<std::vector<double>> dist = this->dijkstra(start, grid, S_GRID_SIZE);
      for(int j = 0; j < W+K; j++) {
        int ti = j < W ? gl::wi[j] / S_GRID_SIZE : gl::hi[j-W] / S_GRID_SIZE;
        int tj = j < W ? gl::wj[j] / S_GRID_SIZE : gl::hj[j-W] / S_GRID_SIZE;
        dists[i][j] = dists[j][i] = dist[ti][tj];
      }
    }
    std::cerr << "FOUND SOURCES DISTANCE (FROM construct_mst)" << std::endl;
    for(int i = 0; i < K; i++) {
      pii start = {gl::hi[i] / S_GRID_SIZE, gl::hj[i] / S_GRID_SIZE};
      std::vector<std::vector<double>> dist = this->dijkstra(start, grid, S_GRID_SIZE);
      for(int j = 0; j < W+K; j++) {
        int ti = j < W ? gl::wi[j] / S_GRID_SIZE : gl::hi[j-W] / S_GRID_SIZE;
        int tj = j < W ? gl::wj[j] / S_GRID_SIZE : gl::hj[j-W] / S_GRID_SIZE;
        dists[i+W][j] = dists[j][i+W] = dist[ti][tj];
      }
    }
    std::cerr << "FOUND HOUSES DISTANCE (FROM construct_mst)" << std::endl;
    std::vector<pair<double, pii>> dist_list;
    for(int i = 0; i < W+K; i++) {
      for(int j = i + 1; j < W+K; j++) {
        dist_list.push_back({dists[i][j], {i, j}});
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
      uf.merge(i, j);
      cnt++;
      edges.push_back({i, j});
      for(int w = 0; w < W; w++) {
        for(int k = 0; k < K; k++) if(uf.same(w, k+W)) boss[W+k] = w;
      }
    }
    
    std::cerr << "FOUND EDGES OF MST (from construct_mst)" << std::endl;
    std::cerr << "EDGE LIST: " << std::endl;
    for(pii edge: edges) {
      std::cerr << edge.first << " " << edge.second << endl;
    }
    for(pii edge: edges) {
      pii start = edge.first < W ? pii{gl::wi[edge.first], gl::wj[edge.first]} : pii{gl::hi[edge.first-W], gl::hj[edge.first-W]};
      pii goal = edge.second < W ? pii{gl::wi[edge.second], gl::wj[edge.second]} : pii{gl::hi[edge.second-W], gl::hj[edge.second-W]};
      start = {start.first / S_GRID_SIZE, start.second / S_GRID_SIZE};
      goal = {goal.first / S_GRID_SIZE, goal.second / S_GRID_SIZE};
      Grid path = dijkstra(start, goal, grid, S_GRID_SIZE);
      for(int i = 0; i < sz; i++) {
        for(int j = 0; j < sz; j++) {
          if(!path[i][j]) continue;
          ret[i][j] = true;
        }
      }
    }

    return {ret, edges};
  }


  Grid construct_mst(const std::vector<std::vector<Area>> &grid, const std::vector<pii> &edges, const std::vector<std::vector<bool>> &valid_area, bool is_original) {
    int W = gl::W;
    int K = gl::K;
    int S_GRID_SIZE = is_original ? 1 : param::SMALL_GRID;
    int sz = grid.size();
    Grid ret(sz, std::vector<bool>(sz, false));
    
    
    for(pii edge: edges) {
      pii start = edge.first < W ? pii{gl::wi[edge.first], gl::wj[edge.first]} : pii{gl::hi[edge.first-W], gl::hj[edge.first-W]};
      pii goal = edge.second < W ? pii{gl::wi[edge.second], gl::wj[edge.second]} : pii{gl::hi[edge.second-W], gl::hj[edge.second-W]};
      start = {start.first / S_GRID_SIZE, start.second / S_GRID_SIZE};
      goal = {goal.first / S_GRID_SIZE, goal.second / S_GRID_SIZE};
      Grid path = dijkstra(start, goal, grid, S_GRID_SIZE);
      for(int i = 0; i < sz; i++) {
        for(int j = 0; j < sz; j++) {
          if(!path[i][j])continue;
          ret[i][j] = true;
        }
      }
    }

    return ret;
  }

  bool hit_on_mst(const Grid &mst, std::vector<std::vector<Area>> &grid, bool is_original, bool random) {
    int sz = mst.size();
    int S_GRID_SIZE = is_original ? 1 : param::SMALL_GRID;
    for(int i = 0; i < sz; i++) {
      for(int j = 0; j < sz; j++) {
        if(!mst[i][j] or grid[i][j].done) continue;
        if(random and mt()%10 < 2) continue;
        int power = (int)this->powers[grid[i][j].cnt];
        power = std::min(power, 5000-grid[i][j].total_power);
        int res = this->excavation(i * S_GRID_SIZE + S_GRID_SIZE / 2, j * S_GRID_SIZE + S_GRID_SIZE/2, power);
        grid[i][j].cnt++;
        grid[i][j].total_power += power;
        if(not is_original) {
          this->num_search++;
          this->search_cost += power;
        }
        if(res == 1) {
          grid[i][j].done = true;
          grid[i][j].stiff_exp = (double)grid[i][j].total_power - (double)power/2;
          if(grid[i][j].cnt == 1) {
            // 一発で終わった場合
            for(int d = 0; d < 4; d++) {
              int ni = i + di[d];
              int nj = j + dj[d];
              if(is_out(ni, nj, gl::N/S_GRID_SIZE, gl::N/S_GRID_SIZE)) continue;
              if(not grid[ni][nj].done) grid[ni][nj].stiff_exp = power;
            }
          }
        } else if(res == 0) {
          grid[i][j].stiff_exp = std::min(grid[i][j].total_power * 3.0, 5000.);
        } else if(res == -1) {
          std::cerr << "Invalid output from hit_on_mst function" << std::endl;
          std::exit(1);
        } else if(res == 2) {
          return true;
        }
      }
    }
    return false;
  }

  bool same_mst(const Grid &a, const Grid &b) {
    int sz = a.size();
    for(int i = 0; i < sz; i++) for(int j = 0; j < sz; j++) if(a[i][j] != b[i][j]) return false;
    return true;
  }

  std::vector<std::vector<double>> dijkstra(const pii& start, const std::vector<std::vector<Area>> &grid, int c) {
    typedef std::pair<double, pair<pii, pii>> Edge;
    int sz = grid.size();
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, 1e18));
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    dist[start.first][start.second] = 0.0;
    pq.push({0., {pii{-1, -1}, start}});
    while(!pq.empty()) {
      Edge edge = pq.top(); pq.pop();
      double cost = edge.first;
      int ci = edge.second.second.first;
      int cj = edge.second.second.second;
      if(cost > dist[ci][cj] + EPS) continue;
      for(int d = 0; d < 4; d++) {
        int ni = ci + di[d];
        int nj = cj + dj[d];
        if(is_out(ni, nj, sz, sz)) continue;
        //int bonus = gl::C*(grid[ci][cj].cnt + grid[ni][nj].cnt)/2;
        double w = (grid[ci][cj].stiff_exp + grid[ni][nj].stiff_exp) / 2 * c + gl::C * c;
        if(dist[ni][nj] > dist[ci][cj] + w + EPS) {
          dist[ni][nj] = dist[ci][cj] + w;
          pq.push({dist[ni][nj], {pii{ci, cj}, pii{ni, nj}}});
        }
      }
    }
    return dist;
  }

  // include valid area arguments
  std::vector<std::vector<double>> dijkstra(const pii& start, const std::vector<std::vector<Area>> &grid, const std::vector<std::vector<bool>> &valid_area) {
    typedef std::pair<double, pair<pii, pii>> Edge;
    int sz = grid.size();
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, 1e18));
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    dist[start.first][start.second] = 0.0;
    pq.push({0., {pii{-1, -1}, start}});
    while(!pq.empty()) {
      Edge edge = pq.top(); pq.pop();
      double cost = edge.first;
      int ci = edge.second.second.first;
      int cj = edge.second.second.second;
      if(cost > dist[ci][cj] + EPS) continue;
      for(int d = 0; d < 4; d++) {
        int ni = ci + di[d];
        int nj = cj + dj[d];
        if(is_out(ni, nj, sz, sz)) continue;
        int bonus = gl::C*(grid[ci][cj].cnt + grid[ni][nj].cnt)/2;
        double w = std::max((grid[ci][cj].stiff_exp + grid[ni][nj].stiff_exp) / 2-bonus, 0.) + gl::C;
        if(!valid_area[ni][nj]) continue;
        if(dist[ni][nj] > dist[ci][cj] + w + EPS) {
          dist[ni][nj] = dist[ci][cj] + w;
          pq.push({dist[ni][nj], {pii{ci, cj}, pii{ni, nj}}});
        }
      }
    }
    return dist;
  }


  // calc Grid path
  Grid dijkstra(const pii& s, const pii& t, const std::vector<std::vector<Area>>& grid, int c) {
    typedef std::pair<double, pair<pii, pii>> Edge;
    int sz = grid.size();
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, 1e18));
    std::vector<std::vector<pii>> pre(sz, std::vector<pii>(sz, {-1, -1}));
    Grid ret(sz, std::vector<bool>(sz, false));
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    dist[s.first][s.second] = 0.0;
    pre[s.first][s.second] = s;
    pq.push({0., {pii{-1, -1}, s}});
    while(!pq.empty()) {
      Edge edge = pq.top(); pq.pop();
      double cost = edge.first;
      int ci = edge.second.second.first;
      int cj = edge.second.second.second;
      if(cost > dist[ci][cj] + EPS or cost > dist[t.first][t.second] + EPS) continue;
      for(int d = 0; d < 4; d++) {
        int ni = ci + di[d];
        int nj = cj + dj[d];
        if(is_out(ni, nj, sz, sz)) continue;
        // int bonus = gl::C*(grid[ci][cj].cnt + grid[ni][nj].cnt)/2;
        double w = (grid[ci][cj].stiff_exp + grid[ni][nj].stiff_exp) / 2 * c + gl::C * c;
        if(dist[ni][nj] > dist[ci][cj] + w + EPS) {
          dist[ni][nj] = dist[ci][cj] + w;
          pq.push({dist[ni][nj], {pii{ci, cj}, pii{ni, nj}}});
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

  // based on pre-defined powers
  std::tuple<double, int, int> excavation_till_break(int i, int j) {
    int total_power = 0;
    int cnt = 0;
    bool broken = false;
    auto it = this->powers.begin();
    int power = (int)(*it);
    while(!broken) {
      power = std::min(power, gl::MAXP - total_power);
      int res;
      std::cout << i << " " << j << " " << power << std::endl;
      #ifdef _LOCAL
      res = this->judge.read_output(i, j, power);
      // std::cerr << "POWER: " << power << ", REMAIN S: " << std::max(judge.S[i][j], 0) << std::endl;
      #else
      std::cin >> res;
      #endif
      total_power += power;
      cnt++;
      if(res == 0) {
        it++;
        power = (int)*it;
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
  int excavation(int i, int j, int power) {
    assert(1 <= power and power <= 5000);
    int res;
    std::cout << i << " " << j << " " << power << std::endl;
    #ifdef _LOCAL
    res = this->judge.read_output(i, j, power);
    #else
    std::cin >> res;
    #endif
    return res;
  }

};

bool arg_parse(int argc, char* argv[]) {
  for(int i = 1; i < argc; i++) {
    if(argv[i] == "--seed") {
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return false;
      } else {
        int seed = std::stoi(argv[i+1]);
        mt = std::mt19937(seed);
      }
      i++;
    } else if(argv[i] == "--start") {
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return false;
      } else {
        for(int j = 0; j < 8; j++) {
          param::starts[j] = std::stoi(argv[i+1]);
        }
      }
      i++;
    } else if(argv[i] == "--coef") {
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return false;
      } else {
        for(int j = 0; j < 8; j++) {
          param::coefs[j] = std::stod(argv[i+1]);
        }
      }
      i++;
    } else if(argv[i] == "--stiff"){
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return false;
      } else {
        for(int j = 0; j < 8; j++) {
          param::INIT_STIF_EXP = std::stoi(argv[i+1]);
        }
      }
      i++;
    } else if(argv[i] == "--C") {
      if(i + 1  >= argc) {
        std::cerr << "Invalid argument" << std::endl;
        return false;
      } else {
        param::pre_defined_C = std::stoi(argv[i+1]);
      }
      i++;

    }
  }
  return true;
}

int main(int argc, char* argv[]) {
  start_time = clock();
  Solver solver;

  #ifdef _OPTUNA
  if(!arg_parse(argc, argv)) return -1;
  #endif

  solver.init();
  solver.solve();
  solver.summary();
  return 0;
}