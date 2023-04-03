// baseline
// 普通の2部マッチング
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
#include <cstring>
#include <memory>
#include <tuple>

// atcoder
#include <atcoder/dsu>
#include <atcoder/maxflow>

// >---------------  macro ---------------------
#define pii std::pair<int,int>
#define Vec std::vector<int>
#define rep(i, s, n) for(int i = s; i < (n); i++)


using namespace std;
using ll = long long;

// >---------------  constants ---------------------
const int INF = 1e9 + 10;
const ll LINF = 1LL<<60;
const double EPS = 1e-6;
const double PI = std::acos(-1);
static const int di[6] = {0, 1, 0, -1, 0, 0};
static const int dj[6] = {1, 0, -1, 0, 0, 0};
static const int dk[6] = {0, 0, 0, 0, 1, -1};
static const int di8[8] = {0, 1, 1,  1,  0, -1, -1, -1};
static const int dj8[8] = {1, 1, 0, -1, -1, -1,  0,  1};
const double TL = 6.0;
const double TL98 = TL * 0.98;
const double TL95 = TL * 0.95;
const double TL90 = TL * 0.90;
const double TL85 = TL * 0.85;



// >---------------  global variables ---------------------
int D;
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
clock_t start_time, cur_time;



// >---------------  functionss  ---------------------
double get_time() {
  cur_time = clock();
  return (double)(cur_time - start_time) / CLOCKS_PER_SEC;
}

bool is_out(int i, int j, int k, int d) {
  return (i < 0 or j < 0 or k < 0 or i >= d or j >= d or k >= d);
}

int pos2ID(int x, int y, int z, int d) {
  return x * d * d + y * d + z;
}

std::tuple<int,int,int> id2Pos(int id, int d) {
  int x, y, z;
  z = id % d;
  id /= d;
  y = id % d;
  id /= d;
  x = id;
  return {x, y, z};
}



// -------------------------------------------------------


// # memo
// 両方のシルエットに適用できるようななるべく大きいブロックを作った方がいい。
// ブロックが余らないようにしたい。
// 2部マッチング


struct Silhouette {
  bool f[15][15]; // (z,x)
  bool r[15][15]; // (z,y)
  Silhouette() {
    for(int i = 0; i < D; i++) {
      for(int j = 0; j < D; j++) {
        f[i][j] = r[i][j] = false;
      }
    }
  }
};


struct Result {
  long long score;
  int n;
  std::vector<int> res1, res2;
};

class Solver {
public:
  Silhouette s1, s2;
  Result best_result;

  Solver() {
  } // constructor

  void init() {
    cerr << "INITIALIZATION STARTED" << endl;
    std::cin >> D;
    std::string buf;
    // silhouette 1
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) s1.f[i][j] = buf[j] == '1' ? true : false; 
    }
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) s1.r[i][j] = buf[j] == '1' ? true : false; 
    }
    

    // silhouette 2
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) s2.f[i][j] = buf[j] == '1' ? true : false; 
    }
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) s2.r[i][j] = buf[j] == '1' ? true : false; 
    }

    cerr << "INITIALIZATION ENDED" << endl;
  }

  // main solve
  void solve() {
    int v1[D][D][D] = {};
    int v2[D][D][D] = {};
    int S = D * D * D;
    int T = D * D * D + 1;
    int maxV = D * D * D;
    atcoder::mf_graph<int> g1(D * D * D + 10), g2(D * D * D + 10);

    for(int x = 0; x < D; x++) {
      for(int y = 0; y < D; y++) {
        for(int z = 0; z < D; z++) {
          if(s1.f[z][x] and s1.r[z][y]) v1[x][y][z] = 1;
          if(s2.f[z][x] and s2.r[z][y]) v2[x][y][z] = 1;
        }
      }
    }

    // construct max flow graph
    for(int x = 0; x < D; x++) {
      for(int y = 0; y < D; y++) {
        for(int z = 0; z < D; z++) {
          int u = pos2ID(x, y, z, D);
          if(v1[x][y][z] == 1) {
            if((x+y+z)%2==0) {
              g1.add_edge(S, u, 1);
              for(int d = 0; d < 6; d++) {
                int nx = x + di[d];
                int ny = y + dj[d];
                int nz = z + dk[d];
                int v = pos2ID(nx, ny, nz, D);
                if(is_out(nx, ny, nz, D)) continue;
                if(v1[nx][ny][nz] == 1) g1.add_edge(u, v, 1);
              }
            } else g1.add_edge(u, T, 1);
          }

          if(v2[x][y][z] == 1) {
            if((x+y+z)%2==0) {
              g2.add_edge(S, u, 1);
              for(int d = 0; d < 6; d++) {
                int nx = x + di[d];
                int ny = y + dj[d];
                int nz = z + dk[d];
                int v = pos2ID(nx, ny, nz, D);
                if(is_out(nx, ny, nz, D)) continue;
                if(v2[nx][ny][nz] == 1) g2.add_edge(u, v, 1);
              }
            } else g2.add_edge(u, T, 1);

          }
        }
      }
    } // construct max flow graph

    // flow
    int f1 = g1.flow(S, T);
    int f2 = g2.flow(S, T);
    std::cerr << "Flow 1: " << f1 << std::endl;
    std::cerr << "Flow 2: " << f1 << std::endl;

    // 
    int min_pair = std::min(f1, f2);
    int maxID1 = min_pair, maxID2 = min_pair;
    std::vector<int> res1(maxV, 0), res2(maxV, 0);
    int cnt = 1;

    for(auto e: g1.edges()) {
      if(cnt > min_pair) break;
      if(e.from == S or e.to == T or e.flow == 0) continue;
      res1[e.from] = cnt;
      res1[e.to] = cnt;
      cnt++;
    }

    cnt = 0;
    for(auto e: g2.edges()) {
      if(cnt > min_pair) break;
      if(e.from == S or e.to == T or e.flow == 0) continue;
      res2[e.from] = cnt;
      res2[e.to] = cnt;
      cnt++;
    }

    // fill remain block with  1*1*1 block
    for(int x = 0; x < D; x++) {
      for(int y = 0; y < D; y++) {
        for(int z = 0; z < D; z++) {
          int u = pos2ID(x, y, z, D);
          if(v1[x][y][z] and res1[u] == 0) {
            res1[u] = maxID1 + 1;
            maxID1++;
          }
          if(v2[x][y][z] and res2[u] == 0) {
            res2[u] = maxID2 + 1;
            maxID2++;
          }
        }
      }
    }

    int n_ans = std::max(maxID1, maxID2);
    long long score = this->calc_score(n_ans, res1, res2);
    best_result = {score, n_ans, res1, res2};

  } // solve

  // return summaries of solve
  void summary() {
    std::cerr << "\n##### SUMMARY #######################\n";
    std::cerr << "ELAPSED TIME : " << get_time() << " s\n";
    std::cerr << "BEST SCORE   : " << this->best_result.score << "\n";
    std::cerr << "##### SUMMARY #######################\n";
  }

  // show output
  void show() {
    std::cout << best_result.n << "\n";
    for(int v: best_result.res1) std::cout << v << " ";
    std::cout << "\n"; 
    for(int v: best_result.res2) std::cout << v << " ";
    std::cout << "\n"; 
  }


  ll calc_score(int n, std::vector<int> &v1, std::vector<int> &v2) {
    long long ret = 0;
    int r = 0;
    double v_tot = 0.0;
    std::vector<int> block1(D * D * D + 1, 0), block2(D * D * D + 1, 0);
    for(int v: v1) block1[v]++;
    for(int v: v2) block2[v]++;
    for(int i = 1; i <= D * D * D; i++) {
      if(block1[i] == 0 or block2[i] == 0) r += block1[i] + block2[i];
      else if(block1[i] == block2[i]) v_tot += 1.0 / block1[i];
    }
    return std::round(1e9 * ((double)r + v_tot));
  }


  double eval_score(){

  }
  
private:

};


bool arg_parse(int argv, char* argc[]) {
  for(int i = 0; i < argv; i++) {
    if(std::strcmp(argc[i], "--seed")) {
      if(i + 1 > argv) {
        std::cerr << "no arguments." << std::endl;
        return false;
      }
      int _seed = std::stoi(argc[i+1]);
      mt = std::mt19937(_seed);
    }
  }
  return true;
}

int main(int argv, char* argc[]) {
  start_time = clock();
  Solver solver;

  #ifdef _OPTUNA
  bool ok = arg_parse(argv, argc)
  #endif

  solver.init();
  solver.solve();
  solver.show();
  solver.summary();
  return 0;
}