// includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <stack>
#include <bitset>

#include <random>
#include <algorithm>
#include <limits>
#include <utility>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <chrono>
#include <memory>
#include <optional>

#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>


// --------------------- macros ----------------------
using namespace std;
#define pii std::pair<int, int>
#define Vec std::vector
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define FOR(i, s, n) for(int i = (int)s; i < (int)(n); i++)

using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int16_t;
using i64 = int16_t;

// --------------------- constants ----------------------
constexpr int INF = 1<<30;
constexpr long long LINF = 1LL << 60;
constexpr double EPS = 1e-6;
const double PI = std::acos(-1);

double TL = 1.98;

// --------------------- global variables ----------------------
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
std::uniform_int_distribution<> rand01(0, 1);
std::uniform_real_distribution<> randReal(0, 1);
clock_t start_time, cur_time;

// >--- ↓ input variables ----- <
int N; // 30<=N<=100
int D; // 2 <=D<=N/4
int Q; // 2N<=Q<=32N
// --------------------- functionss ----------------------

inline double get_time(clock_t start_time)
{
  return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

bool arg_parse(int argv, char *argc[])
{
  for (int i = 0; i < argv; i++)
  {
    if (std::strcmp(argc[i], "--seed") == 0)
    {
      if (i + 1 > argv)
      {
        std::cerr << "no arguments." << std::endl;
        return false;
      }
      int _seed = std::stoi(argc[i + 1]);
      mt = std::mt19937(_seed);
    }
    if (std::strcmp(argc[i], "--TL") == 0)
    {
      if(i + 1 > argv) 
      {
        std::cerr << "No arguments." << std::endl;
      }
      double _TL = std::stod(argc[i + 1]);
      TL = _TL;

    }
  }
  return true;
}

// --------------------- classes ----------------------
enum class Comp {
  LARGER,
  LESS,
  EQUAL
};

struct IOServer {
  int n, d, q;
  vector<int> w;
  IOServer(){}
  IOServer(int n_, int d_, int q_, const std::vector<int> &w_) {
    n = n_;
    d = d_;
    q = q_;
    w = w_;
  }

  Comp query(int nl, int nr, const std::vector<int> &vl, const std::vector<int> &vr) {
    assert(vl.size() != 0 && vr.size() != 0);
    ll lsum = 0LL;
    ll rsum = 0LL;
    for(int idx: vl) lsum += w[idx];
    for(int idx: vr) rsum += w[idx];
    if(lsum > rsum) {
      return Comp::LARGER;
    } else if(lsum < rsum) {
      return Comp::LESS;
    } else {
      return Comp::EQUAL;
    }
  }

  long long calc_score(const std::vector<int> &ans) {
    long long t_sum = 0LL;
    double var = 0.0;
    std::vector<long long> t(D, 0);
    FOR(i, 0, N) t[ans[i]] += w[i];
    FOR(d, 0, D) t_sum += t[d];
    double t_mean = (double)t_sum / D;
    FOR(d, 0, D) var += pow(t[d] - t_mean, 2);
    var /= D;
    return 1 + std::round(100.0 * std::sqrt(var));
  }
};
// ----------------------------------------------------


class Solver
{
public:
#ifdef _LOCAL
  IOServer server;
#endif
  int num_search = 0;

  Solver()
  {
  } // constructor

  void init()
  {
    std::cin >> N >> D >> Q; 

#ifdef _LOCAL
    std::vector<int> w(N);
    FOR(i, 0, N) std::cin >> w[i];
    server = IOServer(N, D, Q, w);
#endif

    fprintf(stderr, "[case] N: %d, D: %d, Q: %d \n",N, D, Q);
    std::cerr << "Initialization done\n";
  }

  //
  void solve()
  {
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::vector<std::pair<int,int>> history;
    FOR(i, 0, N) cluster[i%D].push_back(i);
    std::uniform_int_distribution<> randD(0, D-1);
    std::set<int> fixed_cluster_indices;
    int q = 0;
    while(q < Q) {
      cerr << "query: " << q << endl;
      int nl, nr;
      std::vector<int> vl, vr;
      int largest_cluster_idx;
      int smallest_cluster_idx;
      bool suspend = false;
      vector<vector<int>> graph(D);
      std::vector<int> cluster_indices;
      FOR(d, 0, D) if(fixed_cluster_indices.find(d) == fixed_cluster_indices.end()) cluster_indices.emplace_back(d);
      largest_cluster_idx = cluster_indices.front();
      FOR(d, 0, cluster_indices.size() - 1) {
        int u = largest_cluster_idx;
        int v = cluster_indices[d + 1];
        
        nl = cluster[u].size();
        nr = cluster[v].size();
        vl = cluster[u];
        vr = cluster[v];

        Comp res = query(nl, nr, vl, vr, q);
        if(res == Comp::LARGER) {
          graph[u].push_back(v);
        } else if(res == Comp::LESS) {
          graph[v].push_back(u);
          largest_cluster_idx = v;
        } else {
          graph[v].push_back(u);
          largest_cluster_idx = v;
        }

        if(q == Q && d < D - 1) {
          suspend = true;
        }
        if(q == Q) break;
      }

      if(q == Q) break;
      std::vector<int> leaves;
      FOR(d, 0, D) if(graph[d].empty() && fixed_cluster_indices.find(d) == fixed_cluster_indices.end()) leaves.push_back(d);
      smallest_cluster_idx = leaves.front();
      FOR(d, 0, leaves.size() - 1) {
        int u = smallest_cluster_idx;
        int v = leaves[d + 1];
        
        nl = cluster[u].size();
        nr = cluster[v].size();
        vl = cluster[u];
        vr = cluster[v];

        Comp res = query(nl, nr, vl, vr, q);
        if(res == Comp::LARGER) {
          smallest_cluster_idx = v;
        }

        if(q == Q && d < leaves.size() - 1) {
          suspend = true;
        }
        if(q == Q) break;
      }

      if(!suspend) {
        int from = largest_cluster_idx;
        int to = smallest_cluster_idx;
        if(!history.empty()) {
          std::pair<int,int> pre = history.back();
          if(pre.first == from && pre.second == to) {
            best_cluster = cluster;
          } else if(!(pre.first == to && pre.second == from)) {
            best_cluster = cluster;
          } else {
            cluster = pre_cluster;
            continue;
          }
        }
        int fsize = cluster[from].size();
        int tsize = cluster[to].size();
        if(fsize == 1) {
          fixed_cluster_indices.insert(from);
          continue;
        }
        int idx = (u32)mt() % cluster[from].size();
        int target_item_idx = cluster[from][idx];
        pre_cluster = cluster;
        cluster[from].erase(cluster[from].begin() + idx);
        cluster[to].push_back(target_item_idx);
        history.emplace_back(from, to);
      }
    }

    std::vector<int> ans(N, 0);
    for(int d = 0; d < D; d++) {
      for(int idx: best_cluster[d]) ans[idx] = d;
    }

    for(int i = 0; i < N; i++) std::cout << ans[i] << " \n"[i + 1 == N];
#ifdef _LOCAL
    std::cerr << "Score: " << server.calc_score(ans) << "\n";
#endif

  } // solve

  Comp query(int nl, int nr, const std::vector<int>&vl, const std::vector<int>&vr, int& num_query) {
    num_query++;
    std::cout << nl << " " << nr << " ";
    for(int idx: vl) std::cout << idx << " ";
    for(int idx: vr) std::cout << idx << " ";
    std::cout << std::endl;
#ifdef _LOCAL
    return this->server.query(nl, nr, vl, vr);
#else
    char res;
    std::cin >> res;
    switch (res) {
      case '>':
        return Comp::LARGER;
      case '<':
        return Comp::LESS;
      case '=':
        return Comp::EQUAL;
      default:
        std::cerr << "Invalid result" << std::endl;
        std::exit(-1);
    }
#endif
  } 

  // return summaries of solve
  void summary()
  {
    std::cerr << "\n##### SUMMARY #######################\n";
    fprintf(stderr, "ELAPSED TIME          : %6.4f s\n", get_time(start_time));
    fprintf(stderr, "NUMBER OF SEARCH      : %d\n", this->num_search);
    std::cerr << "######################################\n";
  }


private:
  
};

int main(int argv, char *argc[])
{
  start_time = clock();
  Solver solver;

  solver.init();
  solver.solve();
  solver.summary();
  return 0;
}