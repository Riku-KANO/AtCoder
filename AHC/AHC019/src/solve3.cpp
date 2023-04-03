// sub10
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
#include <optional>
#include <tuple>

// atcoder
#include <atcoder/dsu>
#include <atcoder/maxflow>

// >---------------  macro ---------------------
#define pii std::pair<int,int>
#define Vec std::vector<int>
#define rep(i, s, n) for(int i = s; i < (n); i++)
#define all_space(x, y, z, d) for(int x = 0; x < d; x++) for(int y = 0; y < d; y++) for(int z = 0; z < d; z++)

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
double TL = 6.0;
double TL98 = TL * 0.98;
double TL95 = TL * 0.95;
double TL90 = TL * 0.90;
double TL85 = TL * 0.85;
double TL10 = TL * 0.10;

constexpr int MAX_ROT = 24;
constexpr int MAX_DIM = 3;
constexpr int MAX_DIR = 6;
const int dirs[MAX_ROT][MAX_DIM][MAX_DIR] = {
  {
    {1,-1,0,0,0,0,},
    {0,0,1,-1,0,0,},
    {0,0,0,0,1,-1,},
  },
  {
    {1,-1,0,0,0,0,},
    {0,0,0,0,-1,1,},
    {0,0,1,-1,0,0,},
  },
  {
    {1,-1,0,0,0,0,},
    {0,0,0,0,1,-1,},
    {0,0,-1,1,0,0,},
  },
  {
    {1,-1,0,0,0,0,},
    {0,0,-1,1,0,0,},
    {0,0,0,0,-1,1,},
  },
  {
    {0,0,1,-1,0,0,},
    {1,-1,0,0,0,0,},
    {0,0,0,0,-1,1,},
  },
  {
    {0,0,0,0,1,-1,},
    {1,-1,0,0,0,0,},
    {0,0,1,-1,0,0,},
  },
  {
    {0,0,0,0,-1,1,},
    {1,-1,0,0,0,0,},
    {0,0,-1,1,0,0,},
  },
  {
    {0,0,-1,1,0,0,},
    {1,-1,0,0,0,0,},
    {0,0,0,0,1,-1,},
  },
  {
    {0,0,1,-1,0,0,},
    {0,0,0,0,1,-1,},
    {1,-1,0,0,0,0,},
  },
  {
    {0,0,0,0,-1,1,},
    {0,0,1,-1,0,0,},
    {1,-1,0,0,0,0,},
  },
  {
    {0,0,0,0,1,-1,},
    {0,0,-1,1,0,0,},
    {1,-1,0,0,0,0,},
  },
  {
    {0,0,-1,1,0,0,},
    {0,0,0,0,-1,1,},
    {1,-1,0,0,0,0,},
  },
  {
    {0,0,1,-1,0,0,},
    {0,0,0,0,-1,1,},
    {-1,1,0,0,0,0,},
  },
  {
    {0,0,0,0,1,-1,},
    {0,0,1,-1,0,0,},
    {-1,1,0,0,0,0,},
  },
  {
    {0,0,0,0,-1,1,},
    {0,0,-1,1,0,0,},
    {-1,1,0,0,0,0,},
  },
  {
    {0,0,-1,1,0,0,},
    {0,0,0,0,1,-1,},
    {-1,1,0,0,0,0,},
  },
  {
    {0,0,1,-1,0,0,},
    {-1,1,0,0,0,0,},
    {0,0,0,0,1,-1,},
  },
  {
    {0,0,0,0,-1,1,},
    {-1,1,0,0,0,0,},
    {0,0,1,-1,0,0,},
  },
  {
    {0,0,0,0,1,-1,},
    {-1,1,0,0,0,0,},
    {0,0,-1,1,0,0,},
  },
  {
    {0,0,-1,1,0,0,},
    {-1,1,0,0,0,0,},
    {0,0,0,0,-1,1,},
  },
  {
    {-1,1,0,0,0,0,},
    {0,0,1,-1,0,0,},
    {0,0,0,0,-1,1,},
  },
  {
    {-1,1,0,0,0,0,},
    {0,0,0,0,1,-1,},
    {0,0,1,-1,0,0,},
  },
  {
    {-1,1,0,0,0,0,},
    {0,0,0,0,-1,1,},
    {0,0,-1,1,0,0,},
  },
  {
    {-1,1,0,0,0,0,},
    {0,0,-1,1,0,0,},
    {0,0,0,0,1,-1,},
  },
};

// >---------------  global variables ---------------------
int D;
int MAX_V;
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
clock_t start_time, cur_time;
uniform_real_distribution<> randReal(0,1);

// >---------------  functionss  ---------------------
double get_time() {
  cur_time = clock();
  return (double)(cur_time - start_time) / CLOCKS_PER_SEC;
}

bool is_out(int i, int j, int k, int d) {
  return (i < 0 or j < 0 or k < 0 or i >= d or j >= d or k >= d);
}

template<typename T> T pow2(T x) {
  return x * x;
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

// axis=0: x, axis=1: y, axis=2: z
std::tuple<int,int,int> rotate(int x, int y, int z, int axis=2, int degree=90) {
  if(axis==0) {
    if(degree==90) {
      return {x, -z, y};
    } else if(degree==180) {
      return {x, -y, -z};
    } else if(degree == -90 || degree == 270) {
      return {x, z, -y};
    }
  } else if(axis == 1) {
    if(degree==90) {
      return {-z, y, x};
    } else if(degree==180) {
      return {-x, y, -z};
    } else if(degree == -90 || degree == 270) {
      return {z, y, -x};
    }
  } else if(axis == 2) {
    if(degree==90) {
      return {-y, x, z};
    } else if(degree==180) {
      return {-x, -y, z};
    } else if(degree == -90 || degree == 270) {
      return {y, -x, z};
    }
  } else {
    std::cerr << "Invalid arguments!!\n";
    return {-1,-1,-1};
  }
}

bool arg_parse(int argv, char* argc[]) {
  for(int i = 0; i < argv; i++) {
    if(std::strcmp(argc[i], "--seed") == 0) {
      if(i + 1 > argv) {
        std::cerr << "no arguments." << std::endl;
        return false;
      }
      int _seed = std::stoi(argc[i+1]);
      mt = std::mt19937(_seed);
    }

    if(std::strcmp(argc[i], "--TL") == 0) {
      if(i + 1 > argv) {
        std::cerr << "no arguments." << std::endl;
        return false;
      }
      double _TL = std::stod(argc[i+1]);
      TL = _TL;
      TL98 = TL * 0.98;
      TL95 = TL * 0.95;
      TL90 = TL * 0.90;
      TL85 = TL * 0.85;
      TL10 = TL * 0.10;
    }
  }
  return true;
}

// Rolling hash2d
const int BASE1 = 13;
const int BASE2 = 19;
const int MOD = 1e9 + 7;
std::vector<long long> p_table, q_table;

void prepare(int L) {
    p_table.resize(L+1); q_table.resize(L+1);
    p_table[0] = q_table[0] = 1;
    for (int i = 0; i < L; i++) {
        p_table[i+1] = (1LL * p_table[i] * BASE1) % MOD;
        q_table[i+1] = (1LL * q_table[i] * BASE2) % MOD;
    }
}

std::vector<std::vector<int>> rolling_hash(const std::vector<std::vector<int>>& S) {
    int W = S[0].size();
    int H = S.size();
    std::vector<std::vector<int>> ret(H+1, std::vector<int>(W+1, 0));
    for (int i = 0; i < H; i++) {
        int su = 0;
        std::vector<int>& dp = ret[i];
        std::vector<int>& di = ret[i+1];
        const std::vector<int>& si = S[i];
        for (int j = 0; j < W; j++) {
            int v = si[j];
            su = (1LL * su * BASE1 % MOD + v) % MOD;
            di[j+1] = (su + 1LL * dp[j+1] * BASE2 % MOD) % MOD;
        }
    }
    return ret;
}

/*
  0-indexed
  [i0, i1)×[j1, j2)の範囲のハッシュ値を計算
  S: 2次元データ
  i0: 第1軸の最小値, i1: 第1軸の上限
  j0: 第2軸の最小値, j1: 第2軸の上限
*/
int get_hash(const vector<vector<int>> &S, int i0, int i1, int j0, int j1) {
    int P = p_table[j1 - j0];
    int Q = q_table[i1 - i0];
    return (S[i1][j1] - 1LL * S[i1][j0] * P % MOD - 1LL * S[i0][j1] * Q % MOD + 1LL * S[i0][j0] * (P * Q % MOD) % MOD + MOD) % MOD;
}


std::vector<int> gen_shuffle_arr(int n) {
  std::vector<int> ret(n);
  std::iota(ret.begin(), ret.end(), 0);
  std::shuffle(ret.begin(), ret.end(), mt);
  return ret;
}


// >---------------  classes  ---------------------

/* Lowlink: グラフの関節点・橋を列挙する構造体
    作成: O(E+V)
    関節点の集合: vector<int> aps
    橋の集合: vector<P> bridges
    https://algo-logic.info/articulation-points/
    ※一部改変
*/
struct LowLink {
    const std::vector<std::vector<int>> &G;
    std::vector<int> used, ord, low;
    std::set<int> aps;  // articulation points
    std::vector<pii> bridges;
    LowLink(const std::vector<std::vector<int>> &G_) : G(G_) {
        used.assign(G.size(), 0);
        ord.assign(G.size(), 0);
        low.assign(G.size(), 0);
        int k = 0;
        for (int i = 0; i < (int)G.size(); i++) {
            if (!used[i]) k = dfs(i, k, -1);
        }
        sort(bridges.begin(), bridges.end()); // 必要ならソートする
    }

    int dfs(int id, int k, int par) { // id:探索中の頂点, k:dfsで何番目に探索するか, par:idの親
        used[id] = true;
        ord[id] = k++;
        low[id] = ord[id];
        bool is_aps = false;
        int count = 0; // 子の数
        for (auto &e : G[id]) {
            if (!used[e]) {
                count++;
                k = dfs(e, k, id);
                low[id] = min(low[id], low[e]);
                if (par != -1 && ord[id] <= low[e]) is_aps = true; // 条件2を満たすので関節点
                if (ord[id] < low[e]) bridges.emplace_back(min(id, e), max(id, e));  
            } else if (e != par) { // eが後退辺の時
                low[id] = min(low[id], ord[e]);
            }
        }
        if (par == -1 && count >= 2) is_aps = true; // 条件1を満たすので関節点
        if (is_aps) aps.insert(id);
        return k;
    }
};

// -------------------------------------------------------


// # memo
// 両方のシルエットに適用できるようななるべく大きいブロックを作った方がいい。
// ブロックが余らないようにしたい。
// 2部マッチング（引き算的な考え）
// BFS(足し算的な考え)
// seed=1の2つのオブジェクトの体積をそろえ且つ連結にするのが結構難しい。今後の課題の一つ(連結にしなくてもハイスコアになることが分かった。)
// オブジェクトを連結にするとスコアが悪くなってしまうケースも存在する。（sub5）
// Dが小さいとき、最適値を見つけやすい（探索空間が小さいから）。
// 勝負になるのはDが大きいときかもしれない。
// seed=1, 11のスコアが悪い。
// z軸方向に延びる棒もシルエットを覆う効率が高い。
// 無駄な探索が多いので効率を上げる必要がある。(change_state2の３重ループとか)
// sub8は１回の探索の範囲がひろいのでsub7よりもDが大きいときにスコアが落ちがち。
// 一方でsub8はDが小さいときの解でハイスコアが出せている傾向にある。

// # 考察
// どちらかのオブジェクトのみを考える。任意の３次元上の２点で同じピースの領域に入り得ないペアがある。
// 可能なピースの集合をPとするとき、任意のピースPi(∈P)が占める領域に対して(x1,y1,z1)∉Pi, (x2,y2,z2)∉Piとなるペア(x1,y1,z1), (x2,y2,z2)が存在する。
// つまり、ある地点からBFSをした場合、たどり着くことができない点が存在する。それらの点には別のピースが必ず割り当てられる。
// BFSでたどり着けない点を列挙し続け、それらの点からブロックを構築すると考えやすいかもしれない。（3/23）
// まぐれでハイスコアになった結果を見てみるとピースの収まり方が整っている傾向がある。
// →ある一つのピースを特定の位置に置いたとき、シルエットの頂点数（？）が減っている。角張りが無くなっている傾向がある。

struct Point {
  int x, y, z;
};

struct Block {
  int size;
  std::vector<int> ids;
};

struct Silhouette {
  std::vector<std::vector<int>> f; // (z,x)
  std::vector<std::vector<int>> r; // (z,y)
  std::vector<std::vector<int>> rh_f; // rolling hash
  std::vector<std::vector<int>> rh_r; // rolling hash
  std::vector<bool> unnecessary;
  std::vector<int> pid_list;
  std::vector<int> pm;
};

struct State {
  std::vector<std::vector<std::vector<int>>> exist;
  std::vector<bool> valid; // flag for position
  std::vector<std::vector<int>> h_f;
  std::vector<std::vector<int>> h_r;
  std::vector<std::vector<bool>> f_filled;
  std::vector<std::vector<bool>> r_filled;
  std::vector<int> res;
  int max_vol;
  State() {
    exist.resize(D, std::vector<std::vector<int>>(D, std::vector<int>(D, 0)));
    valid.resize(MAX_V, false);
    h_f.resize(D, std::vector<int>(D, 0));
    h_r.resize(D, std::vector<int>(D, 0));
    f_filled.resize(D, std::vector<bool>(D, false));
    r_filled.resize(D, std::vector<bool>(D, false));

    res.resize(D * D * D, 0);
    max_vol=0;
  }
};


struct Result {
  int n;
  long long score;
  std::vector<int> res1, res2;

  Result() {
    score = 1e18;
  }

  Result(int _n, long long _score, const Vec& _res1, const Vec& _res2) {
    this->n = _n;
    this->score = _score;
    this->res1 = _res1;
    this->res2 = _res2;
  }
};


class Solver {
public:
  Silhouette sil1, sil2;
  int total_vol_max1 = 0;
  int total_vol_max2 = 0;
  int total_vol_min1 = 0;
  int total_vol_min2 = 0;
  Result best_result;
  int num_search = 0;
  Solver() {

  } // constructor

  void init() {
    std::cerr << "INITIALIZATION STARTED\n";

    std::cin >> D;
    sil1.f.resize(D, std::vector<int>(D));
    sil1.r.resize(D, std::vector<int>(D));
    sil2.f.resize(D, std::vector<int>(D));
    sil2.r.resize(D, std::vector<int>(D));
    MAX_V = D * D * D;

    std::string buf;

    // silhouette 1
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) sil1.f[i][j] = buf[j] == '1' ? 1 : 0; 
    }
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) sil1.r[i][j] = buf[j] == '1' ? 1 : 0; 
    }
    

    // silhouette 2
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) sil2.f[i][j] = buf[j] == '1' ? 1 : 0; 
    }
    for(int i = 0; i < D; i++) {
      std::cin >> buf;
      for(int j = 0; j < D; j++) sil2.r[i][j] = buf[j] == '1' ? 1 : 0; 
    }

    // construct rolling hash table
    prepare(D);
    sil1.rh_f = rolling_hash(sil1.f);
    sil1.rh_r = rolling_hash(sil1.r);
    sil2.rh_f = rolling_hash(sil2.f);
    sil2.rh_r = rolling_hash(sil2.r);

    // unnecessary
    sil1.unnecessary.resize(MAX_V, false);
    sil2.unnecessary.resize(MAX_V, false);

    this->init_unnecessary(sil1);
    this->init_unnecessary(sil2);

    // 最大容量と最小容量の確認

    for(int z = 0; z < D; z++) {
      int numf1 = 0;
      int numf2 = 0;
      int numr1 = 0;
      int numr2 = 0;
      for(int x = 0; x < D; x++) {
        if(sil1.f[z][x] > 0) numf1++;
        if(sil2.f[z][x] > 0) numf2++;
      }
      for(int y = 0; y < D; y++) {
        if(sil1.r[z][y] > 0) numr1++;
        if(sil2.r[z][y] > 0) numr2++;
      }
      this->total_vol_max1 += numf1 * numr1;
      this->total_vol_max2 += numf2 * numr2;
      this->total_vol_min1 += std::max(numf1, numr1);
      this->total_vol_min2 += std::max(numf2, numr2);
    }

    // initialization of pid list and permutation array
    all_space(x,y,z,D) {
      int pid = pos2ID(x, y, z, D);
      if(!sil1.unnecessary[pid]) sil1.pid_list.push_back(pid);
      if(!sil2.unnecessary[pid]) sil2.pid_list.push_back(pid);
    }

    sil1.pm.resize(sil1.pid_list.size());
    sil2.pm.resize(sil2.pid_list.size());

    std::iota(sil1.pm.begin(), sil1.pm.end(), 0);
    std::iota(sil2.pm.begin(), sil2.pm.end(), 0);

    std::cerr << "INITIALIZATION ENDED\n";
  }

  /*
    idea: construct and destruct
    形成したブロックを壊してまた作り上げる。
    去年のAHCのestieコンのwriter解に類似した解法をやってみる
    消し方を工夫しないといけない。
  */
  void solve4() {
    std::pair<State, State> best_state;
    State cur_state1, cur_state2;
    State pre_state1, pre_state2;
    std::vector<int> pid_list1, pid_list2;
    long long best_score = 1e18;
    long long pre_score = 1e18;
    int id = 1;
    int target_state = 1;
    double iter_time = 0.0;
    double iter_start_time = get_time();
    double iter_end_time;

    all_space(x,y,z,D) {
      int uid = pos2ID(x, y, z, D);
      if(!sil1.unnecessary[uid]) {
        cur_state1.valid[uid] = true;
        cur_state1.max_vol++;
      }
      if(!sil2.unnecessary[uid]) {
        cur_state2.valid[uid] = true;
        cur_state2.max_vol++;
      }
    }

    std::shuffle(sil1.pm.begin(), sil1.pm.end(), mt);
    std::shuffle(sil2.pm.begin(), sil2.pm.end(), mt);
    std::vector<std::pair<double, int>> tmp1, tmp2;
    for(int &i: sil1.pm) {
      int s1 = sil1.pid_list[i];
      int x, y, z;
      std::tie(x,y,z) = id2Pos(s1, D);
      double score = pow2((double)x-(double)D/2) + pow2((double)y-(double)D/2) + pow2((double)z-(double)D/2);
      tmp1.push_back(std::make_pair(score, s1));
    }

    for(int &i: sil2.pm) {
      int s1 = sil2.pid_list[i];
      int x, y, z;
      std::tie(x,y,z) = id2Pos(s1, D);
      double score = pow2((double)x-(double)D/2) + pow2((double)y-(double)D/2) + pow2((double)z-(double)D/2);
      tmp2.push_back(std::make_pair(score, s1));
    }
    std::sort(tmp1.rbegin(), tmp1.rend());
    std::sort(tmp2.rbegin(), tmp2.rend());

    for(auto& p: tmp1) pid_list1.push_back(p.second);
    for(auto& p: tmp2) pid_list2.push_back(p.second);

    // make initial solution
    // this->make_init_sol(cur_state1, cur_state2);
    // long long init_score = this->calc_score(-1, cur_state1.res, cur_state2.res);
    // best_score = init_score;
    // best_state = {cur_state1, cur_state2};
    // int n_init = 0;
    // for(int pid = 0; pid < MAX_V; pid++) n_init = std::max({n_init, cur_state1.res[pid], cur_state2.res[pid]});
    // std::cerr << "initial score: " << init_score << "\n";
    // std::cerr << n_init << std::endl;
    // for(int v: cur_state1.res) std::cerr << v << " ";
    // std::cerr << "\n";
    // for(int v: cur_state2.res) std::cerr << v << " ";
    // std::cerr << "\n\n";

    while(get_time() + iter_time < TL95) {
      this->num_search++;
      if(target_state == 1) this->change_state4(cur_state1, cur_state2, sil1, sil2, id, pid_list1, pid_list2);
      else if(target_state == 2) this->change_state4(cur_state2, cur_state1, sil2, sil1, id, pid_list2, pid_list1);


      // 座標圧縮
      std::vector<int> id_list;
      int n_id;

      for(int pid = 0; pid < MAX_V; pid++) {
        id_list.push_back(cur_state1.res[pid]);
        id_list.push_back(cur_state2.res[pid]);
      }
      std::sort(id_list.begin(), id_list.end());
      id_list.erase(std::unique(id_list.begin(), id_list.end()), id_list.end());
      n_id = id_list.size() - 1;

      // relabel
      for(int pid = 0; pid < MAX_V; pid++) {
        int id1 = std::lower_bound(id_list.begin(), id_list.end(), cur_state1.res[pid]) - id_list.begin();
        int id2 = std::lower_bound(id_list.begin(), id_list.end(), cur_state2.res[pid]) - id_list.begin();
        cur_state1.res[pid] = id1;
        cur_state2.res[pid] = id2;
      }

      long long score = this->calc_score(0, cur_state1.res, cur_state2.res);
      if(score < best_score) {
        best_score = score;
        best_state = {cur_state1, cur_state2};
        int n_tmp = n_id;

      } 
      if(score < pre_score) {
        pre_state1 = cur_state1;
        pre_state2 = cur_state2;
        pre_score = score;
      } else {
        if(mt() % 20 != 0) {
          cur_state1 = pre_state1;
          cur_state2 = pre_state2;
          
          id_list.clear();
          for(int pid = 0; pid < MAX_V; pid++) {
            id_list.push_back(cur_state1.res[pid]);
            id_list.push_back(cur_state2.res[pid]);
          }
          std::sort(id_list.begin(), id_list.end());
          id_list.erase(std::unique(id_list.begin(), id_list.end()), id_list.end());
          n_id = id_list.size() - 1;
          
        } else {
          pre_score = score;
          pre_state1 = cur_state1;
          pre_state2 = cur_state2;
        }

      }


      int num_single1 = 0;
      int num_single2 = 0;

      std::vector<int> cnt1(MAX_V, 0), cnt2(MAX_V, 0);
      for(int pid = 0; pid < MAX_V; pid++) {
        cnt1[cur_state1.res[pid]]++;
        cnt2[cur_state2.res[pid]]++;
        id = std::max({id, cur_state1.res[pid], cur_state2.res[pid]});
      }

      id++;

      // search bad ids
      std::set<int> bad_ids;
      int largest_id = -1;
      int largest_size = 0;
      int smallest_id = -1;
      int smallest_size = INF;
      std::vector<pii> num_vert_list = this->calc_num_vert(n_id, cur_state1.res, cur_state2.res);
      std::sort(num_vert_list.rbegin(), num_vert_list.rend());
      int diff = n_id - n_id / 2;
      int erase_size = std::min((int)mt() % diff + n_id / 2, n_id);
      for(int i = 0; i <= erase_size; i++) {
        bad_ids.insert(num_vert_list[i].second);
      }
      for(int v = 1; v < MAX_V; v++) {
        if(cnt1[v] > largest_size) {
          largest_id = v;
          largest_size = cnt1[v];
        }
        if(cnt1[v] > 0 && cnt1[v] < smallest_size) {
          smallest_id = v;
          smallest_size = cnt1[v];
        }
      }
      if(mt()%2==0)bad_ids.insert(largest_id);
      bad_ids.insert(smallest_id);


      bad_ids.insert(mt() % id_list.back() + 1); // ランダムに1個選んで消すようにする。
      if(mt()%2==0) bad_ids.insert(mt() % id_list.back() + 1);
      if(mt()%2==0) bad_ids.insert(mt() % id_list.back() + 1);
      if(mt()%2==0) bad_ids.insert(mt() % id_list.back() + 1);
      if(mt()%2==0) bad_ids.insert(mt() % id_list.back() + 1);


      // delete blocks
      all_space(x, y, z, D) {
        int bid1, bid2;
        int pid = pos2ID(x, y, z, D);

        if(bad_ids.count(cur_state1.res[pid])) {
          cur_state1.res[pid] = 0;
          cur_state1.exist[x][y][z] = 0;
        }

        if(bad_ids.count(cur_state2.res[pid])) {
          cur_state2.res[pid] = 0;
          cur_state2.exist[x][y][z] = 0;
        }
      }

      if(num_single1 > num_single2) target_state = 1;
      else if(num_single2 > num_single1) target_state = 2;
      else target_state = mt() % 2 + 1;
      // refine state
      for(int z = 0; z < D; z++) {
        for(int x = 0; x < D; x++) {
          cur_state1.f_filled[z][x] = false;
          cur_state2.f_filled[z][x] = false;
          cur_state1.h_f[z][x] = 0;
          cur_state2.h_f[z][x] = 0;
        }
        
        for(int y = 0; y < D; y++) {
          cur_state1.r_filled[z][y] = false;
          cur_state2.r_filled[z][y] = false;
          cur_state1.h_r[z][y] = 0;
          cur_state2.h_r[z][y] = 0;
        }
      }

      all_space(x, y, z, D) {
        if(cur_state1.exist[x][y][z] > 0) {
          cur_state1.f_filled[z][x] = true;
          cur_state1.r_filled[z][y] = true;
          cur_state1.h_f[z][x]++;
          cur_state1.h_r[z][y]++;
        }
        if(cur_state2.exist[x][y][z] > 0) {
          cur_state2.f_filled[z][x] = true;
          cur_state2.r_filled[z][y] = true;
          cur_state2.h_f[z][x]++;
          cur_state2.h_r[z][y]++;
        }
      }


      std::shuffle(sil1.pm.begin(), sil1.pm.end(), mt);
      std::shuffle(sil2.pm.begin(), sil2.pm.end(), mt);
      tmp1.clear();
      tmp2.clear();
      pid_list1.clear();
      pid_list2.clear();

      for(int &i: sil1.pm) {
        int s1 = sil1.pid_list[i];
        int x, y, z;
        std::tie(x,y,z) = id2Pos(s1, D);
        if(cur_state1.exist[x][y][z] > 0) continue;
        double score = pow2((double)x-(double)D/2) + pow2((double)y-(double)D/2) + pow2((double)z-(double)D/2);
        tmp1.push_back(std::make_pair(score, s1));
      }

      for(int &i: sil2.pm) {
        int s1 = sil2.pid_list[i];
        int x, y, z;
        std::tie(x,y,z) = id2Pos(s1, D);
        if(cur_state2.exist[x][y][z] > 0) continue;
        double score = pow2((double)x-(double)D/2) + pow2((double)y-(double)D/2) + pow2((double)z-(double)D/2);
        tmp2.push_back(std::make_pair(score, s1));
      }
      std::sort(tmp1.rbegin(), tmp1.rend());
      std::sort(tmp2.rbegin(), tmp2.rend());

      for(auto& p: tmp1) pid_list1.push_back(p.second);
      for(auto& p: tmp2) pid_list2.push_back(p.second);
      iter_end_time = get_time();
      iter_time = iter_end_time - iter_start_time;
      iter_start_time = iter_end_time;

    } // while(get_time() < TL90)


    // 座標圧縮
    std::vector<int> id_list;
    auto[best1, best2] = best_state;

    for(int pid = 0; pid < MAX_V; pid++) {
      id_list.push_back(best1.res[pid]);
      id_list.push_back(best2.res[pid]);
    }
    std::sort(id_list.begin(), id_list.end());
    id_list.erase(std::unique(id_list.begin(), id_list.end()), id_list.end());

    // relabel
    for(int pid = 0; pid < MAX_V; pid++) {
      int id1 = std::lower_bound(id_list.begin(), id_list.end(), best1.res[pid]) - id_list.begin();
      int id2 = std::lower_bound(id_list.begin(), id_list.end(), best2.res[pid]) - id_list.begin();
      best1.res[pid] = id1;
      best2.res[pid] = id2;
    }

    int n_ans = id_list.back();
    long long score = this->calc_score(n_ans, best1.res, best2.res);
    best_result.n = n_ans;
    best_result.score = score;
    best_result.res1 = best1.res;
    best_result.res2 = best2.res;
  }


  void reduce_volume(State& st1, State& st2) {
    auto opt_conn = [&](State& state1, const State& state2, const Silhouette& sil1, const Silhouette& sil2) -> void {
      int maxIter = std::abs(state1.max_vol - state2.max_vol);
      while(maxIter--) {
        std::vector<std::vector<int>> graph(MAX_V);
        all_space(x, y, z, D) {
          int uid = pos2ID(x, y, z, D);
          if(state1.exist[x][y][z]) {
            int nx, ny, nz;
            int vid;
            for(int d = 0; d < 6; d++) {
              nx = x + di[d];
              ny = y + dj[d];
              nz = z + dk[d];
              if(is_out(nx, ny, nz, D)) continue;
              vid = pos2ID(nx, ny, nz, D);
              if(state1.exist[nx][ny][nz]) {
                graph[uid].push_back(vid);
              }
            }
          }
        }

        LowLink lowlink(graph);
        for(int& pid: gen_shuffle_arr(MAX_V)) {
          int x, y, z;
          std::tie(x, y, z) = id2Pos(pid, D);
          if(state1.exist[x][y][z] == 0 || lowlink.aps.count(pid)) continue;
          if(state1.h_f[z][x] == 1 || state1.h_r[z][y] == 1) continue;
          state1.exist[x][y][z] = 0;
          state1.valid[pid] = false;
          state1.h_f[z][x]--;
          state1.h_r[z][y]--;
          state1.max_vol--;
          break;
        }
      }

    };

    if(st1.max_vol > st2.max_vol) opt_conn(st1, st2, sil1, sil2);
    else if (st2.max_vol > st1.max_vol) opt_conn(st2, st1, sil2, sil1);
  }


  // return summaries of solve
  void summary() {
    int v1 = 0, v2 = 0;
    int n_single1 = 0, n_single2 = 0;
    int min_vol1 = INF;
    int min_vol2 = INF;
    int max_vol1 = 0;
    int max_vol2 = 0;
    std::vector<int> block1(MAX_V + 1, 0);
    std::vector<int> block2(MAX_V + 1, 0);

    for(int i = 0; i < MAX_V; i++) {
      block1[best_result.res1[i]]++;
      block2[best_result.res2[i]]++;
      if(best_result.res1[i] > 0) v1++;
      if(best_result.res2[i] > 0) v2++;
    }

    for(int bid = 1; bid <= MAX_V; bid++) {
      if(block1[bid] == 1) n_single1++;
      if(block2[bid] == 1) n_single2++;
      if(block1[bid] >= 1) {
        min_vol1 = std::min(min_vol1, block1[bid]);
        max_vol1 = std::max(max_vol1, block1[bid]);
      }
      if(block2[bid] >= 1) {
        min_vol2 = std::min(min_vol2, block1[bid]);
        max_vol2 = std::max(max_vol2, block1[bid]);
      }
    }


    std::cerr << "\n##### SUMMARY #######################\n";
    std::cerr << "ELAPSED TIME     : " << get_time() << " s\n";
    std::cerr << "BEST SCORE       : " << this->best_result.score << "\n";
    std::cerr << "NUMBER OF SEARCH : " << this->num_search << "\n";
    std::cerr << "--- Object 1-------\n";
    std::cerr << "TOTAL VOLUME     : " << v1 << "\n";
    std::cerr << "TOTAL VOLUME MAX : " << this->total_vol_max1 << "\n";
    std::cerr << "TOTAL VOLUME MIN : " << this->total_vol_min1 << "\n";
    std::cerr << "NUM 1*1*1 CUBE   : " << n_single1 << "\n";
    std::cerr << "MAX SIZE CUBE    : " << max_vol1 << "\n";
    std::cerr << "MIN SIZE CUBE    : " << min_vol1 << "\n";
    std::cerr << "--- Object 2-------\n";
    std::cerr << "TOTAL VOLUME  : " << v2 << "\n";
    std::cerr << "TOTAL VOLUME MAX : " << this->total_vol_max2 << "\n";
    std::cerr << "TOTAL VOLUME MIN : " << this->total_vol_min2 << "\n";
    std::cerr << "NUM 1*1*1 CUBE   : " << n_single2 << "\n";
    std::cerr << "MAX SIZE CUBE    : " << max_vol2 << "\n";
    std::cerr << "MIN SIZE CUBE    : " << min_vol2 << "\n";
    std::cerr << "#####################################\n";
  }

  // show output
  void show() {
    std::cout << best_result.n << std::endl;
    for(int v: best_result.res1) std::cout << v << " ";
    std::cout << std::endl; 
    for(int v: best_result.res2) std::cout << v << " ";
    std::cout << std::endl; 
  }


  ll calc_score(int n, std::vector<int> &v1, std::vector<int> &v2) {
    long long ret = 0;
    int r = 0;
    double v_tot = 0.0;
    std::vector<int> block1(D * D * D + 1, 0), block2(D * D * D + 1, 0);
    for(int v: v1) block1[v]++;
    for(int v: v2) block2[v]++;
    for(int i = 1; i <= D * D * D; i++) {
      if(block1[i] == 0 || block2[i] == 0) r += block1[i] + block2[i];
      else if(block1[i] == block2[i]) v_tot += 1.0 / block1[i];
    }
    return std::round(1e9 * ((double)r + v_tot));
  }


  bool judge_WA() {
    int x,y,z;
    for(int i = 0; i < D * D * D; i++) {
      std::tie(x,y,z) = id2Pos(i, D);
      if(this->sil1.f[z][x] == 0 && this->best_result.res1[i] > 0) {
        std::cerr << "Object1: FRONT SILHOUETTE: --OUTLIAR--\n";
        std::cerr << "x,y,z:("<<x<<","<<y<<","<<z<<")\n";
        return false;
      }
      if(this->sil1.r[z][y] == 0 && this->best_result.res1[i] > 0) {
        std::cerr << "Object1: RIGHT SILHOUETTE: --OUTLIAR--\n";
        std::cerr << "x,y,z:("<<x<<","<<y<<","<<z<<")\n";
        return false;
      }

      if(this->sil2.f[z][x] == 0 && this->best_result.res2[i] > 0) {
        std::cerr << "Object2: FRONT SILHOUETTE: --OUTLIAR--\n";
        std::cerr << "x,y,z:("<<x<<","<<y<<","<<z<<")\n";
        return false;
      }
      if(this->sil2.r[z][y] == 0 && this->best_result.res2[i] > 0) {
        std::cerr << "Object2: RIGHT SILHOUETTE: --OUTLIAR--\n";
        std::cerr << "x,y,z:("<<x<<","<<y<<","<<z<<")\n";
        return false;
      }
      
    }
    return true;
  }

  /*
    ローカル調査用
  */
  void find_largest() {
    struct Info {
      int uid, vid;
      int size;
      int rot;
    };

    State state1, state2;
    std::vector<Info> blocks;

    all_space(x,y,z,D) {
      int uid = pos2ID(x, y, z, D);
      if(!sil1.unnecessary[uid]) {
        state1.valid[uid] = true;
        state1.max_vol++;
      }
      if(!sil2.unnecessary[uid]) {
        state2.valid[uid] = true;
        state2.max_vol++;
      }
    }

    std::pair<int,int> best_pos_pair;
    std::pair<int,int> best_fit_pos_pair;
    int best_rot;
    int best_fit_rot;
    int best_block_size = -1;
    int best_fit_num_vert = INF;
    
    for(int s1: this->sil1.pid_list) {
      if(!state1.valid[s1]) continue;
    
      int x1,y1,z1;
      std::tie(x1,y1,z1) = id2Pos(s1,D);

      
      for(int s2: this->sil2.pid_list) {
        if(!state2.valid[s2]) continue;

        int x2, y2, z2;
        std::tie(x2, y2, z2) = id2Pos(s2, D);
        std::vector<bool> used1, used2;

        // search all rotation
        for(int r = 0; r < MAX_ROT; r++) {
          used1.clear();
          used2.clear();
          used1.resize(MAX_V, false);
          used2.resize(MAX_V, false);

          std::queue<pii> q;
          int block_size=1;

          q.push({s1,s2});
          used1[s1] = true;
          used2[s2] = true;

          while(!q.empty()) {
            int cx1,cx2,cy1,cy2,cz1,cz2;
            auto&&[uid1,uid2] = q.front(); q.pop();
            std::tie(cx1,cy1,cz1) = id2Pos(uid1, D);
            std::tie(cx2,cy2,cz2) = id2Pos(uid2, D);

            for(int d = 0; d < 6; d++) {
              int nx1,nx2,ny1,ny2,nz1,nz2;
              int vid1,vid2;
              nx1=cx1+dirs[0][0][d]; ny1=cy1+dirs[0][1][d]; nz1=cz1+dirs[0][2][d];
              nx2=cx2+dirs[r][0][d]; ny2=cy2+dirs[r][1][d]; nz2=cz2+dirs[r][2][d];
              vid1=pos2ID(nx1,ny1,nz1,D);
              vid2=pos2ID(nx2,ny2,nz2,D);
              if(is_out(nx1,ny1,nz1,D) || is_out(nx2,ny2,nz2,D)) continue;
              if(used1[vid1] || used2[vid2]) continue;
              if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
              q.push(pii{vid1,vid2});
              used1[vid1] = true;
              used2[vid2] = true;
              block_size++;
            }
          } // while(!q.empty())



          if(block_size > best_block_size) {
            best_block_size = block_size;
            best_rot = r;
            best_pos_pair = std::make_pair(s1, s2);
          }

          if(block_size > 1) {
            blocks.push_back({s1, s2, block_size, r});
          }


        } // search_rotation

      } // for(int s2: this->sil2.pid_list)

    } // for(int s1: this->sil1.pid_list)

    std::cerr << "MAX BLOCK SIZE: " << best_block_size << "\n";
    std::cerr << "(s1,s2): " << best_pos_pair.first << "," << best_pos_pair.second << "\n\n";
    std::cerr << "FOUND BLOCK PAIR: " << blocks.size() << "\n";

    // 復元
    std::queue<pii> q;
    auto&&[s1, s2] = best_pos_pair;
    std::vector<bool> used1(MAX_V, false), used2(MAX_V, false);

    q.push({s1,s2});
    used1[s1] = true;
    used2[s2] = true;

    while(!q.empty()) {
      auto&&[uid1,uid2] = q.front(); q.pop();
      int cx1,cx2,cy1,cy2,cz1,cz2;
      std::tie(cx1,cy1,cz1) = id2Pos(uid1,D);
      std::tie(cx2,cy2,cz2) = id2Pos(uid2,D);

      for(int d = 0; d < 6; d++) {
        int nx1,nx2,ny1,ny2,nz1,nz2;
        nx1=cx1+dirs[0][0][d]; ny1=cy1+dirs[0][1][d]; nz1=cz1+dirs[0][2][d];
        nx2=cx2+dirs[best_rot][0][d]; ny2=cy2+dirs[best_rot][1][d]; nz2=cz2+dirs[best_rot][2][d];
        int vid1,vid2;
        vid1=pos2ID(nx1,ny1,nz1,D);
        vid2=pos2ID(nx2,ny2,nz2,D);
        if(is_out(nx1,ny1,nz1,D) || is_out(nx2,ny2,nz2,D)) continue;
        if(used1[vid1] || used2[vid2]) continue;
        if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
        q.push(pii{vid1,vid2});
        used1[vid1] = true;
        used2[vid2] = true;
        
      }
    } // while(!q.empty())

    for(int i = 0; i < MAX_V; i++) {
      int x,y,z;
      std::tie(x,y,z) = id2Pos(i,D);
      if(used1[i]) {
        state1.res[i] = 1;
        state1.exist[x][y][z] = 1;
        state1.f_filled[z][x] = true;
        state1.r_filled[z][y] = true;
      }
      if(used2[i]) {
        state2.res[i] = 1;
        state2.exist[x][y][z] = 1;
        state2.f_filled[z][x] = true;
        state2.r_filled[z][y] = true;
      }
    }

    long long score = this->calc_score(1, state1.res, state2.res);
    best_result.n = 1;
    best_result.res1 = state1.res;
    best_result.res2 = state2.res;
    best_result.score = score;


  }


  void find_best_fit() {
    struct Info {
      int uid, vid;
      int size;
      int rot;
    };
    typedef std::pair<int, pii> Cube;
    int MAX_BLOCK_SIZE = 6;
    int maxIter = 10;
    int bid = 1;
    State state1, state2;
    std::vector<Info> blocks;

    all_space(x,y,z,D) {
      int uid = pos2ID(x, y, z, D);
      if(!sil1.unnecessary[uid]) {
        state1.valid[uid] = true;
        state1.max_vol++;
      }
      if(!sil2.unnecessary[uid]) {
        state2.valid[uid] = true;
        state2.max_vol++;
      }
    }

    while(maxIter--) {

      std::pair<int,int> best_pos_pair;
      int best_rot;
      int best_fit_score = INF;
      
      for(int s1: this->sil1.pid_list) {
        if(!state1.valid[s1]) continue;
        

        int x1,y1,z1;
        std::tie(x1,y1,z1) = id2Pos(s1,D);

        if(state1.exist[x1][y1][z1] > 0) continue;

        int g1 = 0, h1 = 0;
        if(!is_out(x1-1, 0, z1, D)) g1 += state1.f_filled[z1][x1-1];
        if(!is_out(x1+1, 0, z1, D)) g1 += state1.f_filled[z1][x1+1];
        if(!is_out(x1, 0, z1-1, D)) g1 += state1.f_filled[z1-1][x1];
        if(!is_out(x1, 0, z1+1, D)) g1 += state1.f_filled[z1+1][x1];
        if(!is_out(0, y1-1, z1, D)) h1 += state1.r_filled[z1][y1-1];
        if(!is_out(0, y1+1, z1, D)) h1 += state1.r_filled[z1][y1+1];
        if(!is_out(0, y1, z1-1, D)) h1 += state1.r_filled[z1-1][y1];
        if(!is_out(0, y1, z1+1, D)) h1 += state1.r_filled[z1+1][y1];
        if(g1 >= 3 || h1 >= 3) continue;
        
        for(int s2: this->sil2.pid_list) {
          if(!state2.valid[s2]) continue;

          int x2, y2, z2;
          std::tie(x2, y2, z2) = id2Pos(s2, D);

          if(state2.exist[x2][y2][z2] > 0) continue;

          int g2 = 0, h2 = 0;
          if(!is_out(x2-1, 0, z2, D)) g2 += state2.f_filled[z2][x2-1];
          if(!is_out(x2+1, 0, z2, D)) g2 += state2.f_filled[z2][x2+1];
          if(!is_out(x2, 0, z2-1, D)) g2 += state2.f_filled[z2-1][x2];
          if(!is_out(x2, 0, z2+1, D)) g2 += state2.f_filled[z2+1][x2];
          if(!is_out(0, y2-1, z2, D)) h2 += state2.r_filled[z2][y2-1];
          if(!is_out(0, y2+1, z2, D)) h2 += state2.r_filled[z2][y2+1];
          if(!is_out(0, y2, z2-1, D)) h2 += state2.r_filled[z2-1][y2];
          if(!is_out(0, y2, z2+1, D)) h2 += state2.r_filled[z2+1][y2];
          if(g2 >= 3 || h2 >= 3) continue;

          std::vector<bool> used1, used2;

          // search all rotation
          for(int rot = 0; rot < MAX_ROT; rot++) {
            used1.clear();
            used2.clear();
            used1.resize(MAX_V, false);
            used2.resize(MAX_V, false);
            std::vector<std::vector<bool>> f_filled1 = state1.f_filled;
            std::vector<std::vector<bool>> r_filled1 = state1.r_filled;
            std::vector<std::vector<bool>> f_filled2 = state2.f_filled;
            std::vector<std::vector<bool>> r_filled2 = state2.r_filled;
            
            std::queue<pii> q;
            int block_size=1;

            q.push({s1,s2});
            used1[s1] = true;
            used2[s2] = true;
            f_filled1[z1][x1] = r_filled1[z1][y1] = true;
            f_filled2[z2][x2] = r_filled2[z2][y2] = true;

            while(!q.empty()) {
              int cx1,cx2,cy1,cy2,cz1,cz2;
              auto&&[uid1,uid2] = q.front(); q.pop();
              std::tie(cx1,cy1,cz1) = id2Pos(uid1, D);
              std::tie(cx2,cy2,cz2) = id2Pos(uid2, D);

              for(int d = 0; d < 6; d++) {
                int nx1,nx2,ny1,ny2,nz1,nz2;
                int vid1,vid2;
                nx1=cx1+dirs[0][0][d]; ny1=cy1+dirs[0][1][d]; nz1=cz1+dirs[0][2][d];
                nx2=cx2+dirs[rot][0][d]; ny2=cy2+dirs[rot][1][d]; nz2=cz2+dirs[rot][2][d];
                vid1=pos2ID(nx1,ny1,nz1,D);
                vid2=pos2ID(nx2,ny2,nz2,D);
                if(is_out(nx1,ny1,nz1,D) || is_out(nx2,ny2,nz2,D)) continue;
                if(used1[vid1] || used2[vid2]) continue;
                if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
                if(state1.exist[nx1][ny1][nz1] > 0 || state2.exist[nx2][ny2][nz2] > 0) continue;
                if(f_filled1[nz1][nx1] && f_filled2[nz2][nx2] && r_filled1[nz1][ny1] && r_filled2[nz2][ny2]) continue;
                q.push(std::make_pair(vid1, vid2));
                used1[vid1] = used2[vid2] = true;
                f_filled1[nz1][nx1] = f_filled2[nz2][nx2] = r_filled1[nz1][ny1] = r_filled2[nz2][ny2] = true;
                block_size++;
                if(block_size >= MAX_BLOCK_SIZE) break;
              }

              if(block_size >= MAX_BLOCK_SIZE) break;
            } // while(!q.empty())

            // find best fit
            std::vector<std::vector<bool>> r(D, std::vector<bool>(D, false));
            std::vector<std::vector<bool>> f(D, std::vector<bool>(D, false));

            int n_vert = 0;
            int xi, yi, zi;

            // object 1
            for(int pid = 0; pid < MAX_V; pid++) {
              std::tie(xi, yi, zi) = id2Pos(pid, D);
              if(state1.valid[pid]) {
                f[zi][xi] = true;
                r[zi][yi] = true;
              }
            }

            for(int pid = 0; pid < MAX_V; pid++) {
              std::tie(xi, yi, zi) = id2Pos(pid, D);
              if(used1[pid] || state1.res[pid] > 0) {
                f[zi][xi] = false;
                r[zi][yi] = false;
              }
            }
              
            int f1, f2, f3, f4;
            for(int z = 0; z <= D; z++) {
              for(int x = 0; x <= D; x++) {
                f1 = is_out(x-1, 0, z-1, D) ? 0 : f[z-1][x-1];
                f2 = is_out(x, 0, z-1, D) ? 0 : f[z-1][x];
                f3 = is_out(x-1, 0, z, D) ? 0 : f[z][x-1];
                f4 = is_out(x, 0, z, D) ? 0 : f[z][x];
                if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
                else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
              }

              for(int y = 0; y <= D; y++) {
                f1 = is_out(0, y-1, z-1, D) ? 0 : r[z-1][y-1];
                f2 = is_out(0, y, z-1, D) ? 0 : r[z-1][y];
                f3 = is_out(0, y-1, z, D) ? 0 : r[z][y-1];
                f4 = is_out(0, y, z, D) ? 0 : r[z][y];
                if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
                else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 

              }
            }

            // object2
            f.clear(); r.clear();
            f.resize(D, std::vector<bool>(D, false));
            r.resize(D, std::vector<bool>(D, false));
            for(int pid = 0; pid < MAX_V; pid++) {
              std::tie(xi, yi, zi) = id2Pos(pid, D);
              if(state2.valid[pid]) {
                f[zi][xi] = true;
                r[zi][yi] = true;
              }
            }

            for(int pid = 0; pid < MAX_V; pid++) {
              std::tie(xi, yi, zi) = id2Pos(pid, D);
              if(used2[pid] || state2.res[pid] > 0) {
                f[zi][xi] = false;
                r[zi][yi] = false;
              }
            }
              
            for(int z = 0; z <= D; z++) {
              for(int x = 0; x <= D; x++) {
                f1 = is_out(x-1, 0, z-1, D) ? 0 : f[z-1][x-1];
                f2 = is_out(x, 0, z-1, D) ? 0 : f[z-1][x];
                f3 = is_out(x-1, 0, z, D) ? 0 : f[z][x-1];
                f4 = is_out(x, 0, z, D) ? 0 : f[z][x];
                if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
                else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
              }

              for(int y = 0; y <= D; y++) {
                f1 = is_out(0, y-1, z-1, D) ? 0 : r[z-1][y-1];
                f2 = is_out(0, y, z-1, D) ? 0 : r[z-1][y];
                f3 = is_out(0, y-1, z, D) ? 0 : r[z][y-1];
                f4 = is_out(0, y, z, D) ? 0 : r[z][y];
                if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
                else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
              }
            }
            

            if(n_vert < best_fit_score) {
              best_fit_score = n_vert;
              best_rot = rot;
              best_pos_pair = std::make_pair(s1, s2);
            }

              

          } // search_rotation

        } // for(int s2: this->sil2.pid_list)

      } // for(int s1: this->sil1.pid_list)

      std::cerr << "MAX fit score: " << best_fit_score << "\n";
      std::cerr << "(s1,s2): " << best_pos_pair.first << "," << best_pos_pair.second << "\n\n";

      if(best_fit_score == INF) break;

      // 復元
      std::queue<pii> q;
      auto&&[s1, s2] = best_pos_pair;
      std::vector<bool> used1(MAX_V, false), used2(MAX_V, false);
      std::vector<std::vector<bool>> f_filled1 = state1.f_filled;
      std::vector<std::vector<bool>> r_filled1 = state1.r_filled;
      std::vector<std::vector<bool>> f_filled2 = state2.f_filled;
      std::vector<std::vector<bool>> r_filled2 = state2.r_filled;

      q.push({s1,s2});
      used1[s1] = true;
      used2[s2] = true;

      while(!q.empty()) {
        auto&&[uid1,uid2] = q.front(); q.pop();
        int cx1,cx2,cy1,cy2,cz1,cz2;
        std::tie(cx1,cy1,cz1) = id2Pos(uid1,D);
        std::tie(cx2,cy2,cz2) = id2Pos(uid2,D);

        for(int d = 0; d < 6; d++) {
          int nx1,nx2,ny1,ny2,nz1,nz2;
          nx1=cx1+dirs[0][0][d]; ny1=cy1+dirs[0][1][d]; nz1=cz1+dirs[0][2][d];
          nx2=cx2+dirs[best_rot][0][d]; ny2=cy2+dirs[best_rot][1][d]; nz2=cz2+dirs[best_rot][2][d];
          int vid1,vid2;
          vid1=pos2ID(nx1,ny1,nz1,D);
          vid2=pos2ID(nx2,ny2,nz2,D);
          if(is_out(nx1,ny1,nz1,D) || is_out(nx2,ny2,nz2,D)) continue;
          if(used1[vid1] || used2[vid2]) continue;
          if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
          if(state1.exist[nx1][ny1][nz1] > 0 || state2.exist[nx2][ny2][nz2] > 0) continue;
          if(f_filled1[nz1][nx1] && f_filled2[nz2][nx2] && r_filled1[nz1][ny1] && r_filled2[nz2][ny2]) continue;
          q.push(std::make_pair(vid1,vid2));
          used1[vid1] = used2[vid2] = true;
          f_filled1[nz1][nx1] = f_filled2[nz2][nx2] = r_filled1[nz1][ny1] = r_filled2[nz2][ny2] = true;
          
        }
      } // while(!q.empty())

      for(int i = 0; i < MAX_V; i++) {
        int x,y,z;
        std::tie(x,y,z) = id2Pos(i, D);
        if(used1[i]) {
          state1.res[i] = bid;
          state1.exist[x][y][z] = bid;
          state1.f_filled[z][x] = true;
          state1.r_filled[z][y] = true;
        }
        if(used2[i]) {
          state2.res[i] = bid;
          state2.exist[x][y][z] = bid;
          state2.f_filled[z][x] = true;
          state2.r_filled[z][y] = true;
        }
      }

      bid++;
    } // while(maxIter--)




    long long score = this->calc_score(1, state1.res, state2.res);
    best_result.n = bid-1;
    best_result.res1 = state1.res;
    best_result.res2 = state2.res;
    best_result.score = score;


  }
  

  void test() {
    State state1, state2;
    all_space(x,y,z,D) {
      int uid = pos2ID(x, y, z, D);
      if(!sil1.unnecessary[uid]) {
        state1.valid[uid] = true;
        state1.max_vol++;
      }
      if(!sil2.unnecessary[uid]) {
        state2.valid[uid] = true;
        state2.max_vol++;
      }
    }
    std::vector<int> starts1 = find_start_points(state1, state2, this->sil1, this->sil2);
    std::vector<int> starts2 = find_start_points(state2, state1, this->sil2, this->sil1);
    std::cerr << "STARTS 1\n";
    for(auto v1: starts1) {
      std::cerr << v1 << " ";
    }
    std::cerr << "\nSTART 2\n";
    for(auto v2: starts2) {
      std::cerr << v2 << " ";
    }
    std::cerr << "\n";

    int cid1 = 0;
    int cid2 = 0;
    int n_ans = 0;
    for(int v: starts1) state1.res[v] = ++cid1;
    for(int v: starts2) state2.res[v] = ++cid2;
    n_ans = std::max(cid1, cid2);
    best_result.n = n_ans;
    best_result.res1 = state1.res;
    best_result.res2 = state2.res;
  }

  std::vector<int> find_start_points(const State& state1, const State& state2, const Silhouette& sil1, const Silhouette& sil2) {
    std::vector<int> ret;
    const int LIMIT = 2;
    int n_pos = sil1.pid_list.size();
    int idx = mt() % n_pos;
    int s1 = sil1.pid_list[idx];
    std::vector<std::vector<int>> f_filled(D, std::vector<int>(D, 0));
    std::vector<std::vector<int>> r_filled(D, std::vector<int>(D, 0));

    ret.push_back(s1);

    auto&&[x_ini, y_ini, z_ini] = id2Pos(s1, D);
    f_filled[z_ini][x_ini] = r_filled[z_ini][y_ini] = true;

    std::vector<bool> used1, used2;

    while(true) {
      for(int s2: sil2.pid_list) {
        
        for(int r = 0; r < MAX_ROT; r++) {
          used1.clear(); used2.clear();
          used1.resize(MAX_V, false); used2.resize(MAX_V, false);

          std::queue<pii> q;
          q.push(std::make_pair(s1, s2));
          used1[s1] = used2[s2] = true;

          while(!q.empty()) {
            auto&&[u1, u2] = q.front(); q.pop();
            int x1,x2,y1,y2,z1,z2;
            int nx1,nx2,ny1,ny2,nz1,nz2;
            int v1, v2;

            std::tie(x1, y1, z1) = id2Pos(u1, D);
            std::tie(x2, y2, z2) = id2Pos(u2, D);

            for(int d = 0; d < 6; d++) {
              nx1 = x1 + dirs[0][0][d]; ny1 = y1 + dirs[0][1][d]; nz1 = z1 + dirs[0][2][d];
              nx2 = x2 + dirs[r][0][d]; ny2 = y2 + dirs[r][1][d]; nz2 = z2 + dirs[r][2][d];
              v1 = pos2ID(nx1, ny1, nz1, D);
              v2 = pos2ID(nx2, ny2, nz2, D);

              if(is_out(nx1, ny1, nz1, D) || is_out(nx2, ny2, nz2, D)) continue;
              if(used1[v1] || used2[v2] || !state1.valid[v1] || !state2.valid[v2]) continue;

              used1[v1] = used2[v2] = true;
              f_filled[nz1][nx1]++;
              r_filled[nz1][ny1]++;
              q.push(std::make_pair(v1, v2));
            }
          }
        }
      }

      bool ok = true;
      for(int z = 0; z < D && ok; z++) {
        for(int x = 0; x < D; x++) if(sil1.f[z][x] == 1 && f_filled[z][x] < LIMIT) ok = false;
        for(int y = 0; y < D; y++) if(sil1.r[z][y] == 1 && r_filled[z][y] < LIMIT) ok = false;
      }

      if(ok) break;

      for(int i: gen_shuffle_arr(n_pos)) {
        int nx_s1 = sil1.pid_list[i];
        int x, y, z;
        std::tie(x, y, z) = id2Pos(nx_s1, D);
        if(f_filled[z][x] >= LIMIT && r_filled[z][y] >= LIMIT) continue;
        s1 = nx_s1;
        f_filled[z][x]++;
        r_filled[z][y]++;
        ret.push_back(s1);
        break;
      }
    }

    return ret;
  }

  

  std::vector<pii> calc_num_vert(int n_vert, const std::vector<int> res1, const std::vector<int> res2) {
    std::vector<pii> ret;
    std::vector<std::vector<int>> r, f;

    for(int bid = 1; bid <= n_vert; bid++) {
      int n_vert = 0;
      int xi, yi, zi;

      // object 1
      f = sil1.f;
      r = sil1.r;
      for(int pid = 0; pid < MAX_V; pid++) {
        std::tie(xi, yi, zi) = id2Pos(pid, D);
        if(res1[pid] == bid) {
          f[zi][xi] = 0;
          r[zi][yi] = 0;
        }
      }
        
      int f1, f2, f3, f4;
      for(int z = 0; z <= D; z++) {
        for(int x = 0; x <= D; x++) {
          f1 = is_out(x-1, 0, z-1, D) ? 0 : f[z-1][x-1];
          f2 = is_out(x, 0, z-1, D) ? 0 : f[z-1][x];
          f3 = is_out(x-1, 0, z, D) ? 0 : f[z][x-1];
          f4 = is_out(x, 0, z, D) ? 0 : f[z][x];
          if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
          else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
        }

        for(int y = 0; y <= D; y++) {
          f1 = is_out(0, y-1, z-1, D) ? 0 : r[z-1][y-1];
          f2 = is_out(0, y, z-1, D) ? 0 : r[z-1][y];
          f3 = is_out(0, y-1, z, D) ? 0 : r[z][y-1];
          f4 = is_out(0, y, z, D) ? 0 : r[z][y];
          if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
          else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 

        }
      }

      // object2
      f = sil2.f;
      r = sil2.r;
      for(int pid = 0; pid < MAX_V; pid++) {
        std::tie(xi, yi, zi) = id2Pos(pid, D);
        if(res2[pid] == bid) {
          f[zi][xi] = 0;
          r[zi][yi] = 0;
        }
      }
        
      for(int z = 0; z <= D; z++) {
        for(int x = 0; x <= D; x++) {
          f1 = is_out(x-1, 0, z-1, D) ? 0 : f[z-1][x-1];
          f2 = is_out(x, 0, z-1, D) ? 0 : f[z-1][x];
          f3 = is_out(x-1, 0, z, D) ? 0 : f[z][x-1];
          f4 = is_out(x, 0, z, D) ? 0 : f[z][x];
          if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
          else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
        }

        for(int y = 0; y <= D; y++) {
          f1 = is_out(0, y-1, z-1, D) ? 0 : r[z-1][y-1];
          f2 = is_out(0, y, z-1, D) ? 0 : r[z-1][y];
          f3 = is_out(0, y-1, z, D) ? 0 : r[z][y-1];
          f4 = is_out(0, y, z, D) ? 0 : r[z][y];
          if(f1 + f2 + f3 + f4 == 1 || f1 + f2 + f3 + f4 == 3) n_vert++;
          else if(f1 + f2 + f3 + f4 == 2 && f1 == f4) n_vert += 2; 
        }
      }

      ret.push_back(std::make_pair(n_vert, bid));
      
    }


    return ret;
  }
  
private:
  void init_unnecessary(Silhouette& sil) {
    atcoder::dsu uf(MAX_V);
    all_space(x, y, z, D) {
      if(sil.f[z][x] == 0 || sil.r[z][y] == 0) continue;
      int nx, ny, nz;
      int uid = pos2ID(x, y, z, D);
      nx = x + 1; ny = y + 1, nz = z + 1;
      if(!is_out(nx, y, z, D)) if(sil.f[z][nx] && sil.r[z][y])uf.merge(uid, uid + D * D);
      if(!is_out(x, ny, z, D)) if(sil.f[z][x] && sil.r[z][ny])uf.merge(uid, uid + D);
      if(!is_out(x, y, nz, D)) if(sil.f[nz][x] && sil.r[nz][y])uf.merge(uid, uid + 1);
    }


    for(auto &g: uf.groups()) {
      if(g.size() <= 5) {
        for(int v: g) sil.unnecessary[v] = true;
      }
    }
  }

  /*
    最大流で初期解を作る。
  */
  void make_init_sol(State& state1, State& state2) {
    int S = MAX_V + 10;
    int T = MAX_V + 11;
    int MAX_N = MAX_V + 20;
    atcoder::mf_graph<int> graph1(MAX_N), graph2(MAX_N);

    all_space(x, y, z, D) {
      int uid = pos2ID(x, y, z, D);
      int vid;
      if((x+y+z)%2==0) {
        if(state1.valid[uid]) {
          graph1.add_edge(S, uid, 1);
          for(int d = 0; d < 6; d++) {
            int nx = x + di[d];
            int ny = y + dj[d];
            int nz = z + dk[d];
            vid = pos2ID(nx, ny, nz, D);
            if(is_out(nx, ny, nz, D)) continue;
            if(state1.valid[vid]) graph1.add_edge(uid, vid, 1);
          }
        }
        if(state2.valid[uid]) {
          graph2.add_edge(S, uid, 1);
          for(int d = 0; d < 6; d++) {
            int nx = x + di[d];
            int ny = y + dj[d];
            int nz = z + dk[d];
            vid = pos2ID(nx, ny, nz, D);
            if(is_out(nx, ny, nz, D)) continue;
            if(state2.valid[vid]) graph2.add_edge(uid, vid, 1);
          }
        }
      } else {
        if(state1.valid[uid]) graph1.add_edge(uid, T, 1);
        if(state2.valid[uid]) graph2.add_edge(uid, T, 1);
      }

    } // all_space

    // flow
    int flow1 = graph1.flow(S, T);
    int flow2 = graph2.flow(S, T);


    int cid1 = 0;
    int cid2 = 0;

    for(auto& e: graph1.edges()) {
      if(e.from == S || e.to == T || e.flow == 0) continue;
      cid1++;
      state1.res[e.from] = cid1;
      state1.res[e.to] = cid1;
      int x, y, z;
      std::tie(x, y, z) = id2Pos(e.from, D);
      state1.exist[x][y][z] = cid1;
      std::tie(x, y, z) = id2Pos(e.to, D);
      state1.exist[x][y][z] = cid1;
    }

    for(auto& e: graph2.edges()) {
      if(e.from == S || e.to == T || e.flow == 0) continue;
      cid2++;
      state2.res[e.from] = cid2;
      state2.res[e.to] = cid2;
      int x, y, z;
      std::tie(x, y, z) = id2Pos(e.from, D);
      state2.exist[x][y][z] = cid2;
      std::tie(x, y, z) = id2Pos(e.to, D);
      state2.exist[x][y][z] = cid2;
    }

    all_space(x, y, z, D) {
      if(state1.exist[x][y][z] > 0) {
        state1.f_filled[z][x] = true;
        state1.r_filled[z][y] = true;
      }
      if(state2.exist[x][y][z] > 0) {
        state2.f_filled[z][x] = true;
        state2.r_filled[z][y] = true;
      }
    }
    
    int maxid = std::max(cid1, cid2);
    cid1 = maxid;
    cid2 = maxid;
    // check whether all space are filled correctly
    all_space(x, y, z, D) {
      int pid = pos2ID(x, y, z, D);
      if(state1.valid[pid]) {
        if(!state1.f_filled[z][x]) {
          state1.res[pid] = ++cid1;
          state1.exist[x][y][z] = cid1;
          state1.f_filled[z][x] = true;
          state1.r_filled[z][y] = true;
        }

        if(!state1.r_filled[z][y]) {
          state1.res[pid] = ++cid1;
          state1.exist[x][y][z] = cid1;
          state1.f_filled[z][x] = true;
          state1.r_filled[z][y] = true;
        }
      }

      if(state2.valid[pid]) {
        if(!state2.f_filled[z][x]) {
          state2.res[pid] = ++cid2;
          state2.exist[x][y][z] = cid2;
          state2.f_filled[z][x] = true;
          state2.r_filled[z][y] = true;
        }

        if(!state2.r_filled[z][y]) {
          state2.res[pid] = ++cid2;
          state2.exist[x][y][z] = cid2;
          state2.f_filled[z][x] = true;
          state2.r_filled[z][y] = true;
        }
      }
    }


    all_space(x, y, z, D) {
      if(state1.exist[x][y][z] > 0) {
        state1.h_f[z][x]++;
        state1.h_r[z][y]++;
        state1.f_filled[z][x] = true;
        state1.r_filled[z][y] = true;
      }

      if(state2.exist[x][y][z] > 0) {
        state2.h_f[z][x]++;
        state2.h_r[z][y]++;
        state2.f_filled[z][x] = true;
        state2.r_filled[z][y] = true;
      }
    }
  } 


  void change_state4(State& state1, State& state2, Silhouette& sil1, Silhouette& sil2, int start_id, const std::vector<int>& pid_list1, const std::vector<int>& pid_list2) {
    int maxIter = 20;
    int id = start_id + 1;
    std::vector<int> d_list = {0, 1, 2, 3, 4, 5};
    while(maxIter--) {
      bool changed = false;
      std::shuffle(d_list.begin(), d_list.end(), mt);

      for(const int& i: pid_list1) {
        int s1 = i;
        // int s1 = sil1.pid_list[i];
        int x1, y1, z1;
        std::tie(x1,y1,z1) = id2Pos(s1,D);
        if(state1.exist[x1][y1][z1] > 0) continue;
        if(state1.f_filled[z1][x1] && state1.r_filled[z1][y1]) continue;
        if(!state1.valid[s1]) continue;

        std::tuple<int,int,int> best_pos;
        int max_block_size = -1;
        int best_rot;
        
        // against object2
        for(const int&j : pid_list2) {
          int s2 = j;
          // int s2 = sil2.pid_list[j];
          int x2,y2,z2;
          std::tie(x2,y2,z2) = id2Pos(s2, D);

          if(state2.exist[x2][y2][z2] > 0) continue;
          if(state2.f_filled[z2][x2] && state2.r_filled[z2][y2]) continue;
          if(!state2.valid[s2]) continue;

          std::vector<bool> used1, used2;

          // search all rotation
          for(int r = 0; r < MAX_ROT; r++) {
            used1.clear();
            used2.clear();
            used1.resize(MAX_V, false);
            used2.resize(MAX_V, false);

            std::queue<pii> q;
            int block_size=1;

            q.push({s1,s2});
            used1[s1] = true;
            used2[s2] = true;

            while(!q.empty()) {
              int cx1,cx2,cy1,cy2,cz1,cz2;
              auto&&[uid1,uid2] = q.front(); q.pop();
              std::tie(cx1,cy1,cz1) = id2Pos(uid1,D);
              std::tie(cx2,cy2,cz2) = id2Pos(uid2,D);

              for(int d = 0; d < 6; d++) {
                int nx1,nx2,ny1,ny2,nz1,nz2;
                int vid1,vid2;
                nx1=cx1+dirs[0][0][d]; ny1=cy1+dirs[0][1][d]; nz1=cz1+dirs[0][2][d];
                nx2=cx2+dirs[r][0][d]; ny2=cy2+dirs[r][1][d]; nz2=cz2+dirs[r][2][d];
                vid1=pos2ID(nx1,ny1,nz1,D);
                vid2=pos2ID(nx2,ny2,nz2,D);
                if(is_out(nx1,ny1,nz1,D) || is_out(nx2,ny2,nz2,D)) continue;
                if(state1.exist[nx1][ny1][nz1] || state2.exist[nx2][ny2][nz2]) continue;
                if(used1[vid1] || used2[vid2]) continue;
                if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
                q.push(pii{vid1,vid2});
                used1[vid1] = true;
                used2[vid2] = true;
                block_size++;
              }
            } // while(!q.empty())

            if(block_size > max_block_size) {
              max_block_size = block_size;
              best_rot = r;
              best_pos = {x2,y2,z2};
            }
          } // search_rotation
          if(max_block_size > 5 * D) break;
        } // all_space 2


        if(max_block_size == -1 || max_block_size == 1) continue;
        
        changed = false;
        
        int x2,y2,z2;
        int s2;
        std::tie(x2,y2,z2) = best_pos;
        s2 = pos2ID(x2,y2,z2,D);
        std::vector<bool> used1(MAX_V, false), used2(MAX_V, false);
        std::vector<int> cubes1, cubes2;
        std::vector<pii> used_list;
        used_list.push_back(std::make_pair(s1, s2));
        used1[s1] = true;
        used2[s2] = true;
        cubes1.push_back(s1);
        cubes2.push_back(s2);
        int block_size = 1;
        int block_size_limit; // = std::min(max_block_size, (int)std::sqrt(state1.max_vol));
        if(max_block_size > D ) {
            int diff = max_block_size - D;
            block_size_limit = D + mt() % diff + 1;
        } else {
          block_size_limit = max_block_size;
        }

        while(block_size < block_size_limit) {
          std::vector<pii> candidates;
          for(auto&&[uid1, uid2]: used_list) {
            int cx1, cx2, cy1, cy2, cz1, cz2;
            int nx1, nx2, ny1, ny2, nz1, nz2;
            int vid1, vid2;
            std::tie(cx1, cy1, cz1) = id2Pos(uid1, D);
            std::tie(cx2, cy2, cz2) = id2Pos(uid2, D);

            for(int d = 0; d < 6; d++) {
              nx1 = cx1 + dirs[0][0][d]; ny1 = cy1 + dirs[0][1][d]; nz1 = cz1 + dirs[0][2][d];
              nx2 = cx2 + dirs[best_rot][0][d]; ny2 = cy2 + dirs[best_rot][1][d]; nz2 = cz2 + dirs[best_rot][2][d];
              if(is_out(nx1, ny1, nz1, D) || is_out(nx2, ny2, nz2, D)) continue;
              vid1 = pos2ID(nx1, ny1, nz1, D);
              vid2 = pos2ID(nx2, ny2, nz2, D);
              if(!state1.valid[vid1] || !state2.valid[vid2]) continue;
              if(state1.exist[nx1][ny1][nz1] || state2.exist[nx2][ny2][nz2]) continue;
              if(used1[vid1] || used2[vid2]) continue;
              candidates.push_back(std::make_pair(vid1, vid2));
            }
          }

          int c_size = candidates.size();
          int choosed_idx = mt() % c_size;
          pii p = candidates[choosed_idx];
          cubes1.push_back(p.first);
          cubes2.push_back(p.second);
          used1[p.first] = true;
          used2[p.second] = true;
          used_list.push_back(p);
          block_size++;
        }


        for(int &pid: cubes1) {
          int x, y, z;
          std::tie(x,y,z) = id2Pos(pid, D);
          state1.res[pid] = id;
          state1.exist[x][y][z] = id;
          state1.f_filled[z][x] = true;
          state1.r_filled[z][y] = true;
        }

        for(int &pid: cubes2) {
          int x, y, z;
          std::tie(x, y, z) = id2Pos(pid, D);
          state2.res[pid] = id;
          state2.exist[x][y][z] = id;
          state2.f_filled[z][x] = true;
          state2.r_filled[z][y] = true;
        }

        id++;

      } // all_space 1


      bool ok = true;
      for(int z = 0; z < D && ok; z++) {
        for(int x = 0; x < D && ok; x++) {
          if(sil1.f[z][x] && !state1.f_filled[z][x]) ok = false;
          if(sil2.f[z][x] && !state2.f_filled[z][x]) ok = false;
        }

        for(int y = 0; y < D && ok; y++) {
          if(sil1.r[z][y] && !state1.r_filled[z][y]) ok = false;
          if(sil2.r[z][y] && !state2.r_filled[z][y]) ok = false;
        }
      }
      if(ok || !changed) break;
      
    } // while(true)

    // fill 1*1*1 block
    int sid1 = id;
    int sid2 = id;
    int x, y, z;
    for(int &uid: gen_shuffle_arr(MAX_V)) {
      std::tie(x, y, z) = id2Pos(uid, D);

      if(!state1.f_filled[z][x] && sil1.f[z][x] && sil1.r[z][y] && !sil1.unnecessary[uid]) {
        state1.res[uid] = sid1++;
        state1.f_filled[z][x] = true;
        state1.r_filled[z][y] = true;
      }

      if(!state1.r_filled[z][y] && sil1.f[z][x] && sil1.r[z][y] && !sil1.unnecessary[uid]) {
        state1.res[uid] = sid1++;
        state1.f_filled[z][x] = true;
        state1.r_filled[z][y] = true;
      }

      if(!state2.f_filled[z][x] && sil2.f[z][x] && sil2.r[z][y] && !sil2.unnecessary[uid]) {
        state2.res[uid] = sid2++;
        state2.f_filled[z][x] = true;
        state2.r_filled[z][y] = true;
      }

      if(!state2.r_filled[z][y] && sil2.f[z][x] && sil2.r[z][y] && !sil2.unnecessary[uid]) {
        state2.res[uid] = sid2++;
        state2.f_filled[z][x] = true;
        state2.r_filled[z][y] = true;
      }

    } // while(maxIter--)

  } // change_state4
};



int main(int argv, char* argc[]) {
  start_time = clock();
  Solver solver;

  #ifdef _OPTUNA
  bool ok = arg_parse(argv, argc)
  if(not ok) {
    std::cerr << "Invalid arguments!\n";
    return -1;
  }
  #endif

  #ifdef ONLINE_JUDGE
  solver.init();
  solver.solve4();
  solver.show();

  #elif defined(TEST)
  bool ok = arg_parse(argv, argc);
  if(not ok) {
    std::cerr << "Incvalid arguments!\n";
    return -1;
  }
  solver.init();
  solver.test();
  solver.summary();
  solver.show();

  #else
  solver.init();
  solver.solve4();
  solver.judge_WA();
  solver.show();
  solver.summary();

  #endif

  return 0;
}