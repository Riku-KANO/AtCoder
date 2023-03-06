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
#include <memory>
#include <cstdio>
#include <sstream>
#include <tuple>
#include <optional>

// atcoder
#include <atcoder/fenwicktree>
#include <atcoder/segtree>

// boost
#include <boost/polygon/polygon.hpp>
#include <boost/geometry.hpp>

using namespace std;
namespace bg = boost::geometry;
using ll = long long;

#define pii std::pair<int,int>
#define Vec std::vector<int>
#define rep(i, s, n) for(int i = s; i < (n); i++)
#define BBox bg::model::box
#define Cartesian bg::cs::cartesian
#define BPoint bg::model::point

// # memo -------
// Vmaxの範囲はVを全体積としたとき、０．４V～０．９V程度。0.9Vのときが難しそうではある。
// ri =0 の時、横(wpi)の辺が横軸、縦(hpi)の辺が縦軸に平行になるように置く。
// ri =1 の時、縦hpi)の辺が横軸、横(wpi)の辺が縦軸に平行になるように置く。
// ri =2 の時、高さ(dpi)の辺が横軸、縦(hpi)の辺が縦軸に平行になるように置く。
// ri =3 の時、縦(hpi)の辺が横軸、高さ(dpi)の辺が縦軸に平行になるように置く。
// ri =4 の時、高さ(dpi)の辺が横軸、横(wpi)の辺が縦軸に平行になるように置く。
// ri =5 の時、横(wpi)の辺が横軸、高さ(dpi)の辺が縦軸に平行になるように置く。

// # idea
// 強化学習が出来そうだと思ったが実行時間短すぎて使えない（１つの推論で0.2sぐらい）。
// C++で解を作ってその教師データをDQNに食わせて状態価値を推論させることはできるかも。
// それができたら無駄な探索が減る可能性があり、速くベスト解が見つかる。
// 
//---------------

// ----------------------- constant variables ------------------------
static const int di[4] = {0, 1, 0, -1};
static const int dj[4] = {1, 0, -1, 0};
static const int di8[8] = {0, 1, 1,  1,  0, -1, -1, -1};
static const int di8[8] = {1, 1, 0, -1, -1, -1,  0,  1};
const int INF = 1e9 + 10;
const ll LINF = 1LL<<60;
const double EPS = 1e-8;
const double PI = std::acos(-1);
const double MIN_OVERLAP_RATIO = 0.60;
const int NUM_ACTION = 6;
const int MAX_WIDTH = 1120;
const int MAX_HEIGHT = 680;
int MAX_DEPTH;
#ifdef ONLINE_JUDGE
const double TL = 2.0;
#elif _TRAINING
const double TL = 100.0;
#else
const double TL = 2.0;
#endif
const double TL80=TL*0.80;
const double TL85=TL*0.85;
const double TL90=TL*0.90;
const double TL95=TL*0.95;
const double TL98=TL*0.98;
// ------------------------ utility -------------------------------------
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
std::uniform_real_distribution<> dist01(0,1);
clock_t start_time, cur_time;

// ------------------------ functions -----------------------------------

double get_time() {
  return (double)(cur_time - start_time) / CLOCKS_PER_SEC;
}

bool is_out(int i, int j, int H, int W) {
  return (i < 0 or j < 0 or i >= H or j >= W);
}

// ----------------------- global variables -----------------------------

int M, W, H, B, D;

// ----------------------------------------------------------------------

struct Vec3D {
  int x, y, z;
};

Vec3D operator+(const Vec3D& lhs, const Vec3D& rhs) {
  return Vec3D{lhs.x + rhs.y, lhs.y + rhs.y, lhs.z + rhs.z};
}

Vec3D operator-(const Vec3D& lhs, const Vec3D& rhs) {
  return Vec3D{lhs.x + rhs.y, lhs.y + rhs.y, lhs.z + rhs.z};
}


struct Baggage {
  int id;
  int h, w, d;
  bool f; // can rotate around z-axis
  bool g; // can load on other baggage;
  vector<bool> valid_act;

  Baggage() {
    this->valid_act.resize(NUM_ACTION, true);
  }

  std::tuple<int,int,int> rotate(int r) {
    if(not valid_act[r]) {
      std::cerr << "The rotation is not allowed for baggage(" << bid << ")" << std::endl;
      std::exit(-1);
      return {h, w, d};
    }

    switch(r) {
      case 0: // hwd -> hwd
        return {h, w, d};

      case 1: // hwd -> whd
        return {w, h, d};

      case 2: // hwd -> hdw
        return {h, d, w};

      case 3: // hwd -> dhw
        return {d, h, w};

      case 4: // hwd -> wdh
        return {w, d, h};

      case 5: // hwd -> dwh
        return {d, w, h};

      default:
      {
        std::cerr << "Invalid rotation" << std::endl;
        std::exit(-1);
        return {h, w, d};
      }
    }
  } // rotate function

};

struct Action {
  int bid;
  int r;
  int x, y, z;
};

struct Result {
  long long score;
  std::vector<Action> actions;
};

class State {
public:
  int score;
  std::vector<BBox> loaded_box;
  std::vector<int> remain_box; // number of remaining of each box id;
  std::vector<Action> actions; // done actions
  std::vector<Vec3D> pos_cand; // next candidate position
  atcoder::fenwick_tree<int> bit;
  int max_height;
  int min_height;

  State() {

  } // constructor

private:

};


class Solver {
public:
  int N;
  std::vector<Baggage> baggages;
  Result best_result;

  Solver() {
    
  } // constructor

  void init() {
    cerr << "INITIALIZATION STARTED" << endl;
    std::cin >> M >> W >> H >> B >> D;
    MAX_DEPTH = D;
    N = 0;

    for(int bid = 0; bid < M; bid++) {
      char f, g;
      Baggage b;
      int h, w, d;
      int num;

      b.id = bid;
      std::cin >> b.h >> b.w >> b.d >> num >> f >> g;
      b.f = f == 'Y' ? true : false;
      b.g = g == 'Y' ? true : false;
      if(b.f) {
        b.valid_act[2] = b.valid_act[3] = b.valid_act[4] = b.valid_act[5] = false;
      }
      for(int i = 0; i < num; i++) {
        baggages.push_back(b);
      }

      N += b.num;


    }
    cerr << "INITIALIZATION ENDED" << endl;
  }

  // main solve
  void solve() {


  } // solve

  // trainig solver for DQN
  void solve_for_training() {

  } // solve for training

  // return summaries of solve
  void summary() {
    std::cerr << "### SUMMARY ##################" << std::endl;
    // std::cerr << "BEST SCORE       : " << this->best_score << std::endl;
    std::cerr << "ELAPSED TIME     : " << get_time() << " s" << std::endl;
    std::cerr << "##############################" << std::endl;
  }

  // show output
  void show() {
    assert(best_result.actions == this->N);
    for(Action a: this->best_result.actions) {
      std::cout << a.bid << " " << a.r << " " << a.x << " " << a.y << " " << a.z << std::endl;
    }
  }


  ll calc_score() {

  }


  double eval_score(){

  }
  
private:

};


bool arg_parse(int argc, char* argv[]) {

}

int main(int argc, char* argv[]) {
  start_time = clock();
  Solver solver();

  #ifdef _OPTUNA
  for(int i = 1; i < argc; i++) {
    if(std::strcmp(argv[i], "--seed")) {
      if(i + 1 < argc) {
        int seed = std::stoi(argv[i+1]);
        mt = mt19937(seed);
      }
    } 
  }
  #endif

  solver.init();

  #ifdef _TRAINING
  solver.solve_for_training();
  #else
  solver.solve();
  solver.show();
  solver.summary();
  #endif
  return 0;
}