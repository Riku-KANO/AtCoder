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
#define FOR(i, s, n) for (int i = (int)s; i < (int)(n); i++)

using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int16_t;
using i64 = int16_t;

// --------------------- constants ----------------------
constexpr int INF = 1 << 30;
constexpr long long LINF = 1LL << 60;
constexpr double EPS = 1e-6;
const double PI = std::acos(-1);
constexpr double LAMBDA = 1e-5;

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
      if (i + 1 > argv)
      {
        std::cerr << "No arguments." << std::endl;
      }
      double _TL = std::stod(argc[i + 1]);
      TL = _TL;
    }
  }
  return true;
}

std::vector<int> generate_weights(int num_item, int num_div, int seed_)
{
  std::mt19937 rng(seed_);
  std::exponential_distribution<> dist(LAMBDA);
  std::vector<int> ret;
  const double THRESHOLD = 1e5 * num_item / num_div;
  while (ret.size() < num_item)
  {
    double w_prime = dist(rng);
    if (w_prime > THRESHOLD)
      continue;
    ret.emplace_back(std::max(1, (int)std::round(w_prime)));
  }
  return ret;
}

// --------------------- classes ----------------------
enum class Comp
{
  LARGER,
  LESS,
  EQUAL
};

struct IOServer
{
  int n, d, q;
  vector<int> w;
  IOServer() {}
  IOServer(int n_, int d_, int q_, const std::vector<int> &w_)
  {
    n = n_;
    d = d_;
    q = q_;
    w = w_;
  }

  Comp query(int nl, int nr, const std::vector<int> &vl, const std::vector<int> &vr) const
  {
    assert(vl.size() != 0 && vr.size() != 0 && nl == vl.size() && nr == vr.size());
    ll lsum = 0LL;
    ll rsum = 0LL;
    for (int idx : vl)
      lsum += w[idx];
    for (int idx : vr)
      rsum += w[idx];
    if (lsum > rsum)
    {
      return Comp::LARGER;
    }
    else if (lsum < rsum)
    {
      return Comp::LESS;
    }
    else
    {
      return Comp::EQUAL;
    }
  }

  long long calc_score(const std::vector<int> &ans)
  {
    long long t_sum = 0LL;
    double var = 0.0;
    std::vector<long long> t(D, 0);
    FOR(i, 0, N)
    t[ans[i]] += w[i];
    FOR(d, 0, D)
    t_sum += t[d];
    double t_mean = (double)t_sum / D;
    FOR(d, 0, D)
    var += pow(t[d] - t_mean, 2);
    var /= D;
    return 1 + std::round(100.0 * std::sqrt(var));
  }
};

struct Operation
{
  int u; // from
  int v; // to
  std::vector<int> present_from_u;
  std::vector<int> present_from_v;
  int last_largest;
  int last_smallest;
};

// ----------------------------------------------------

class Solver
{
public:
  IOServer server;
  int num_search = 0;
  std::vector<std::vector<int>> initial_cluster;

  Solver()
  {
  } // constructor

  void init()
  {
    std::cin >> N >> D >> Q;
    initial_cluster.resize(D);
#ifdef _LOCAL
    std::vector<int> w(N);
    FOR(i, 0, N)
    std::cin >> w[i];
    server = IOServer(N, D, Q, w);
#endif
    /**
     *  貪欲法で理想的な割り当てを求める
     */
    std::vector<ll> weight_pred(N);
    for (int seed = 1000; seed < 2000; seed++)
    {
      std::vector<int> weights = generate_weights(N, D, seed);
      std::sort(weights.begin(), weights.end());
      FOR(i, 0, N)
      weight_pred[i] += weights[i];
    }
    FOR(i, 0, N)
    weight_pred[i] /= 1000;
    std::reverse(weight_pred.begin(), weight_pred.end());
    std::vector<ll> cluster_sum(D, 0LL);
    FOR(i, 0, N)
    {
      int lowest_pos = -1;
      ll lowest_val = LINF;
      FOR(d, 0, D)
      {
        if (cluster_sum[d] < lowest_val)
        {
          lowest_val = cluster_sum[d];
          lowest_pos = d;
        }
      }
      cluster_sum[lowest_pos] += weight_pred[i];
      initial_cluster[lowest_pos].push_back(N - 1 - i);
    }
    FOR(d, 0, D)
    std::reverse(initial_cluster[d].begin(), initial_cluster[d].end());

    fprintf(stderr, "[case] N: %d, D: %d, Q: %d \n", N, D, Q);
    std::cerr << "Initialization done\n";
  }

  //
  void solve()
  {
    const int NUM_SIMULATION = 200;
    int res = simulation(NUM_SIMULATION);
    std::cerr << "simulation done\n";
    std::cerr << "best result is solver" << res << std::endl;
    std::vector<int> ans;
    switch (res)
    {
    case 1:
      ans = solve1(this->server, false);
      break;
    case 2:
      ans = solve2(this->server, false);
      break;
    case 3:
      ans = solve3(this->server, false);
      break;
    case 4:
      ans = solve4(this->server, false);
      break;
    case 5:
      ans = solve5(this->server, false);
      break;
    case 6:
      ans = solve6(this->server, false);
      break;
    case 7:
      ans = solve7(this->server, false);
      break;
    case 8:
      ans = solve8(this->server, false);
      break;
    default:
      std::cerr << "Invalid result\n";
      break;
    }
    std::cout << "# best result is solver " << res << std::endl;

    for (int i = 0; i < N; i++)
      std::cout << ans[i] << " \n"[i + 1 == N];
#ifdef _LOCAL
    std::cerr << "Score: " << server.calc_score(ans) << "\n";
#endif

  } // solve

  void find_largest(
      int &largest_idx,
      bool &suspend,
      const std::vector<std::vector<int>> &cluster,
      const std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation)
  {
    int nl, nr;
    std::vector<int> vl, vr;
    std::vector<int> cluster_indices;
    FOR(d, 0, D)
    if (fixed.find(d) == fixed.end())
      cluster_indices.emplace_back(d);
    largest_idx = cluster_indices.front();
    FOR(d, 0, cluster_indices.size() - 1)
    {
      int u = largest_idx;
      int v = cluster_indices[d + 1];

      nl = cluster[u].size();
      nr = cluster[v].size();
      vl = cluster[u];
      vr = cluster[v];

      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);

      if (res == Comp::LESS)
      {
        largest_idx = v;
      }

      if (num_query == Q && d < D - 1)
      {
        suspend = true;
        return;
      }
      if (num_query == Q)
      {
        return;
      }
    }
  }

  void find_smallest(
      int &smallest_idx,
      bool &suspend,
      const std::vector<std::vector<int>> &cluster,
      const std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation)
  {
    int nl, nr;
    std::vector<int> vl, vr;
    std::vector<int> cluster_indices;
    FOR(d, 0, D)
    if (fixed.find(d) == fixed.end())
      cluster_indices.emplace_back(d);
    smallest_idx = cluster_indices.front();
    FOR(d, 0, cluster_indices.size() - 1)
    {
      int u = smallest_idx;
      int v = cluster_indices[d + 1];

      nl = cluster[u].size();
      nr = cluster[v].size();
      vl = cluster[u];
      vr = cluster[v];

      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LARGER)
      {
        smallest_idx = v;
      }

      if (num_query == Q && d < D - 1)
      {
        suspend = true;
        return;
      }
      if (num_query == Q)
      {
        return;
      }
    }
  }

  /**
   * @brief
   *
   * @param largest_idx 総和が最も大きいクラスターのインデックス
   * @param smallest_idx 総和が最も小さいクラスターのインデックス
   * @param cluster
   * @param fixed
   * @param num_query クエリ回数
   * @return true クエリ回数がオーバーする場合,true
   * @return false クエリ回数がオーバーしていない場合,false
   */
  void find_largest_and_smallest(
      int &largest_idx,
      int &smallest_idx,
      bool &suspend,
      const std::vector<std::vector<int>> &cluster,
      const std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation)
  {
    int nl, nr;
    std::vector<int> vl, vr;
    vector<vector<int>> graph(D);

    /**
     * largest indexを見つける
     *
     */
    std::vector<int> cluster_indices;
    FOR(d, 0, D)
    if (fixed.find(d) == fixed.end())
      cluster_indices.emplace_back(d);
    largest_idx = cluster_indices.front();
    FOR(d, 0, cluster_indices.size() - 1)
    {
      int u = largest_idx;
      int v = cluster_indices[d + 1];

      nl = cluster[u].size();
      nr = cluster[v].size();
      vl = cluster[u];
      vr = cluster[v];

      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LARGER)
      {
        graph[u].push_back(v);
      }
      else if (res == Comp::LESS)
      {
        graph[v].push_back(u);
        largest_idx = v;
      }
      else
      {
        graph[v].push_back(u);
        largest_idx = v;
      }

      if (num_query == Q && d < D - 1)
      {
        suspend = true;
        return;
      }
      if (num_query == Q)
      {
        return;
      }
    }

    /**
     * smallest indexを見つける
     *
     */
    std::vector<int> leaves;
    FOR(d, 0, D)
    if (graph[d].empty() && fixed.find(d) == fixed.end())
      leaves.push_back(d);
    smallest_idx = leaves.front();
    FOR(d, 0, leaves.size() - 1)
    {
      int u = smallest_idx;
      int v = leaves[d + 1];

      nl = cluster[u].size();
      nr = cluster[v].size();
      vl = cluster[u];
      vr = cluster[v];

      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LARGER)
      {
        smallest_idx = v;
      }

      if (num_query == Q && d < leaves.size() - 1)
      {
        suspend = true;
        return;
      }
      if (num_query == Q)
        return;
    }
  }

  Comp query(int nl, int nr, const std::vector<int> &vl, const std::vector<int> &vr, int &num_query, const IOServer &server, bool simulation = false)
  {
    num_query++;
    if (simulation)
      return server.query(nl, nr, vl, vr);
    else
    {
      std::cout << nl << " " << nr << " ";
      for (int idx : vl)
        std::cout << idx << " ";
      for (int idx : vr)
        std::cout << idx << " ";
      std::cout << std::endl;
#ifdef _LOCAL
      return server.query(nl, nr, vl, vr);
#else
      char res;
      std::cin >> res;
      switch (res)
      {
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
  /**
   * @brief 二分探索で交換するアイテムを見つける
   *
   * @param from
   * @param to
   * @param suspend
   * @param cluster
   * @param num_query
   * @param server
   * @param is_simulation
   * @return int
   */
  int find_move_item_by_bs(
      int from,
      int to,
      bool &suspend,
      const std::vector<std::vector<int>> &cluster,
      int &num_query,
      const IOServer &server,
      bool is_simulation)
  {
    int ret = -1;
    int l = 0;
    int r = cluster[from].size();
    int m = (l + r) / 2;
    int nl, nr;
    std::vector<int> vl, vr;
    while (l + 1 != r)
    {
      if (num_query == Q)
      {
        suspend = true;
        return ret;
      }
      m = (l + r) / 2;
      vl = cluster[to];
      vr = cluster[from];
      vl.emplace_back(cluster[from][m]);
      vr.erase(vr.begin() + m);
      nl = vl.size();
      nr = vr.size();
      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LARGER)
      {
        r = m;
      }
      else if (res == Comp::LESS)
      {
        l = m;
      }
      else if (res == Comp::EQUAL)
      {
        return cluster[from][m];
      }
    }

    if (l == 0)
    {
      if (num_query == Q)
      {
        suspend = true;
        return ret;
      }
      vl = cluster[to];
      vr = cluster[from];
      vl.emplace_back(cluster[from][0]);
      vr.erase(vr.begin());
      nl = vl.size();
      nr = vr.size();
      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LESS || res == Comp::EQUAL)
      {
        ret = cluster[from][0];
      }
      else
      {
        ret = -1;
      }
    }
    else if (l == cluster[from].size() - 1)
    {
      return ret;
    }
    else
    {
      ret = cluster[from][l];
    }
    return ret;
  }

  /**
   * @brief 二分探索で挿入すべき箇所を見つけ、挿入する
   *
   * @param to
   * @param item
   * @param cluster
   * @param suspend
   * @param num_query
   * @param server
   * @param is_simulation
   */
  void insert_item_by_bs(
      int to,
      int item,
      std::vector<std::vector<int>> &cluster,
      bool &suspend,
      int &num_query,
      const IOServer &server,
      bool is_simulation)
  {
    int l = 0;
    int r = cluster[to].size();
    int m;
    int nl, nr;
    std::vector<int> vl, vr;
    int pos = -1;
    bool contain_less = false;
    bool contain_larger = false;
    while (l + 1 != r)
    {
      if (num_query == Q)
      {
        suspend = true;
        return;
      }
      m = (l + r) / 2;
      nl = 1;
      nr = 1;
      vl = {item};
      vr = {cluster[to][m]};
      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LESS)
      {
        contain_less = true;
        r = m;
      }
      else if (res == Comp::LARGER)
      {
        contain_larger = true;
        l = m;
      }
      else
      {
        pos = m;
        break;
      }
    }
    if (!contain_less)
    {
      if (num_query == Q)
      {
        suspend = true;
        return;
      }
      nl = 1;
      nr = 1;
      vl = {item};
      vr = {cluster[to].back()};
      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LARGER)
      {
        pos = cluster[to].size();
      }
      else
      {
        pos = cluster[to].size() - 1;
      }
    }
    else if (!contain_larger)
    {
      if (num_query == Q)
      {
        suspend = true;
        return;
      }
      nl = 1;
      nr = 1;
      vl = {item};
      vr = {cluster[to].front()};
      Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res == Comp::LESS)
      {
        pos = 0;
      }
      else
      {
        pos = 1;
      }
    }
    else
    {
      pos = l;
    }
    cluster[to].insert(cluster[to].begin() + pos, item);
  }

  int simulation(int num_simulation = 10)
  {
    assert(num_simulation > 0);
    int num_solvers = 8;
    std::vector<ll> scores(num_solvers, 0);
    for (int sid = 0; sid < num_simulation; sid++)
    {
      std::cerr << "simulation iter: " << sid + 1 << std::endl;
      std::vector<int> W = generate_weights(N, D, sid);
      IOServer sim_server(N, D, Q, W);
      std::vector<int> ans;

      // solve1
      ans = this->solve1(sim_server, true);
      ll score1 = sim_server.calc_score(ans);
      scores[0] += score1;

      // solve2
      ans = this->solve2(sim_server, true);
      ll score2 = sim_server.calc_score(ans);
      scores[1] += score2;

      // solve3
      ans = this->solve3(sim_server, true);
      ll score3 = sim_server.calc_score(ans);
      scores[2] += score3;

      // solve4
      ans = this->solve4(sim_server, true);
      ll score4 = sim_server.calc_score(ans);
      scores[3] += score4;

      // solve5
      ans = this->solve5(sim_server, true);
      ll score5 = sim_server.calc_score(ans);
      scores[4] += score5;

      // solve6
      ans = this->solve6(sim_server, true);
      ll score6 = sim_server.calc_score(ans);
      scores[5] += score6;

      // solve7
      ans = this->solve7(sim_server, true);
      ll score7 = sim_server.calc_score(ans);
      scores[6] += score7;

      // solve8
      ans = this->solve8(sim_server, true);
      ll score8 = sim_server.calc_score(ans);
      scores[7] += score8;
    }

    int best_solver = -1;
    ll best_score = LINF;
    for (int i = 0; i < num_solvers; i++)
    {
      scores[i] /= num_simulation;
      std::cerr << "solver" << i + 1 << " score: " << scores[i] << std::endl;
      if (scores[i] < best_score)
      {
        best_score = scores[i];
        best_solver = i + 1;
      }
    }
    return best_solver;
  }
  /**
   * @brief sub4の解法。
   * 大きいクラスターから小さいクラスターに一つアイテムを移動。
   * 条件を満たさなければひとつ前のクラスターの状態に戻る。
   * 14 ms 程度。
   * @param IOServer シミュレーションではない場合、適当なserverが渡される。
   * @param is_simulation
   * @return std::vector<int> ans
   */
  std::vector<int> solve1(const IOServer &server, bool is_simulation)
  {
    std::vector<int> ans(N, 0);
    std::vector<int> ans_tmp(N, 0);
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), mt);
    FOR(i, 0, N)
    cluster[i % D].push_back(indices[i]);
    best_cluster = cluster;
    std::set<int> fixed_cluster_indices;
    int q = 0;
    while (q < Q)
    {
      int nl, nr;
      std::vector<int> vl, vr;
      int largest_cluster_idx;
      int smallest_cluster_idx;
      bool suspend = false;

      find_largest_and_smallest(
          largest_cluster_idx,
          smallest_cluster_idx,
          suspend,
          cluster,
          fixed_cluster_indices,
          q,
          server,
          is_simulation);

      if (!suspend)
      {
        int from = largest_cluster_idx;
        int to = smallest_cluster_idx;
        int fsize = cluster[from].size();
        int tsize = cluster[to].size();
        if (fsize == 1)
        {
          fixed_cluster_indices.insert(from);
          continue;
        }

        int idx = (u32)mt() % fsize;
        int target_item_idx = cluster[from][idx];
        pre_cluster = cluster;
        cluster[from].erase(cluster[from].begin() + idx);
        cluster[to].push_back(target_item_idx);

        std::vector<int> vl, vr;
        vl = pre_cluster[from];
        vr = cluster[to];
        std::tie(vl, vr) = delete_intersection_set(vl, vr);
        nl = vl.size();
        nr = vr.size();
        Comp res = query(nl, nr, vl, vr, q, server, is_simulation);
        if (res == Comp::LARGER)
        {
          best_cluster = cluster;
          if (!is_simulation)
          {
            output_for_vis(best_cluster);
          }
        }
        else
        {
          cluster = pre_cluster;
        }
      }
    }

    for (int d = 0; d < D; d++)
    {
      for (int idx : best_cluster[d])
        ans[idx] = d;
    }
    return ans;
  }

  /**
   * @brief 大きいクラスターと小さいクラスターで1個交換する
   *
   * @param server IOServer
   * @param is_simulation bool
   * @return std::vector<int>
   */
  std::vector<int> solve2(const IOServer &server, bool is_simulation)
  {
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::vector<Operation> history;
    std::vector<int> ans_tmp(N);
    FOR(i, 0, N)
    {
      cluster[i % D].push_back(i);
    }
    best_cluster = cluster;
    std::set<int> fixed_cluster_indices;
    int q = 0;
    while (q < Q)
    {
      int nl, nr;
      std::vector<int> vl, vr;
      int largest_cluster_idx;
      int smallest_cluster_idx;
      bool suspend = false;

      find_largest_and_smallest(
          largest_cluster_idx,
          smallest_cluster_idx,
          suspend,
          cluster,
          fixed_cluster_indices,
          q,
          server,
          is_simulation);

      if (!suspend)
      {
        int from = largest_cluster_idx;
        int to = smallest_cluster_idx;
        if (cluster[from].size() == 1)
        {
          fixed_cluster_indices.insert(from);
          continue;
        }
        bool bad_pre_operation = false;
        if (!history.empty())
        {
          Operation pre_operation = history.back();
          if (pre_operation.last_largest == from && pre_operation.last_smallest == to)
          {
            best_cluster = cluster;
          }
          else if (pre_operation.last_largest == to && pre_operation.last_smallest == from)
          {
            bad_pre_operation = true;
          }
          else if (pre_operation.v == largest_cluster_idx || pre_operation.u == smallest_cluster_idx)
          {
            bad_pre_operation = true;
          }
          else
          {
            best_cluster = cluster;
          }
        }
        if (bad_pre_operation)
        {
          cluster = pre_cluster;
          from = history.back().last_largest;
          to = history.back().last_smallest;
          history.pop_back();
          bool ok = false;
          while (!ok)
          {
            nl = 1;
            nr = 2;
            int fidx = (u32)mt() % cluster[from].size();
            vl = {cluster[from][fidx]};

            /**
             * この実装を入れるとスコアが悪くなる。
             *
             */
            // to = from;
            // while (to == from)
            // {
            //   to = (u32)mt() % D;
            //   if (cluster[to].size() == 1)
            //   {
            //     to = from;
            //     continue;
            //   }
            // }

            std::vector<int> range(cluster[to].size());
            std::iota(range.begin(), range.end(), 0);
            std::shuffle(range.begin(), range.end(), mt);
            vr.clear();
            int tsize = (u32)mt() % 2 + 1;
            for (int i = 0; i < std::min(tsize, (int)range.size()); i++)
            {
              vr.push_back(cluster[to][range[i]]);
            }

            nr = vr.size();
            Comp res = query(nl, nr, vl, vr, q, server, is_simulation);

            if (res == Comp::LESS)
            {
              ok = true;
              pre_cluster = cluster;
              for (int item_idx : vl)
              {
                FOR(j, 0, cluster[from].size())
                {
                  if (cluster[from][j] == item_idx)
                  {
                    cluster[from].erase(cluster[from].begin() + j);
                    break;
                  }
                }
              }
              for (int item_idx : vr)
              {
                FOR(j, 0, cluster[to].size())
                {
                  if (cluster[to][j] == item_idx)
                  {
                    cluster[to].erase(cluster[to].begin() + j);
                    break;
                  }
                }
              }
              for (int item_idx : vr)
                cluster[from].emplace_back(item_idx);
              for (int item_idx : vl)
                cluster[to].emplace_back(item_idx);

              history.push_back(Operation{from, to, vl, vr, largest_cluster_idx, smallest_cluster_idx});
            }
            if (q == Q)
              break;
          }
        }
        else
        {
          int fsize = cluster[from].size();
          int tsize = cluster[to].size();
          int idx = (u32)mt() % fsize;
          int target_item_idx = cluster[from][idx];
          pre_cluster = cluster;
          cluster[from].erase(cluster[from].begin() + idx);
          cluster[to].push_back(target_item_idx);
          Operation operation = {from, to, {target_item_idx}, {}, largest_cluster_idx, smallest_cluster_idx};
          history.push_back(operation);
        }

        if (!is_simulation)
        {
          output_for_vis(best_cluster);
        }
      }
    }
    std::vector<int> ans(N, 0);
    for (int d = 0; d < D; d++)
    {
      for (int idx : best_cluster[d])
        ans[idx] = d;
    }

    return ans;
  }

  /**
   * @brief 最初にマージソートで順位を求める。一つのアイテムだけ移動
   * (Qが小さい場合にはこのソルバーは無効)
   *
   * @param server
   * @param is_simulation
   * @return std::vector<int>
   */
  std::vector<int> solve3(const IOServer &server, bool is_simulation)
  {
    std::vector<int> items(N);
    std::vector<int> ans(N, 0);
    std::iota(items.begin(), items.end(), 0);
    int q = 0;
    bool suspend = false;
    std::vector<int> sorted_items = merge_sort(items, q, suspend, server, is_simulation);
    std::vector<int> sorted_item_indices(N);
    FOR(i, 0, N)
    sorted_item_indices[sorted_items[i]] = i;
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> pre_cluster;
    std::vector<std::vector<int>> best_cluster;
    std::set<int> fixed_cluster_indices;
    std::vector<int> ans_tmp(N);

    FOR(d, 0, D)
    {
      for (int idx : initial_cluster[d])
      {
        cluster[d].push_back(sorted_items[idx]);
      }
    }

    if (q < Q)
      best_cluster = cluster;
    else
      return ans;

    while (q < Q)
    {
      int nl, nr;
      std::vector<int> vl, vr;
      int largest_cluster_idx;
      int smallest_cluster_idx;
      bool suspend = false;

      find_largest_and_smallest(
          largest_cluster_idx,
          smallest_cluster_idx,
          suspend,
          cluster,
          fixed_cluster_indices,
          q,
          server,
          is_simulation);

      if (suspend)
        break;

      int from = largest_cluster_idx;
      int to = smallest_cluster_idx;
      int fsize = cluster[from].size();
      int tsize = cluster[to].size();
      if (fsize == 1)
      {
        fixed_cluster_indices.insert(from);
        continue;
      }

      int target_item = find_move_item_by_bs(from, to, suspend, cluster, q, server, is_simulation);
      if (suspend)
        break;
      if (target_item == -1)
      {
        continue;
      }

      FOR(i, 0, cluster[from].size())
      {
        if (cluster[from][i] == target_item)
        {
          cluster[from].erase(cluster[from].begin() + i);
          break;
        }
      }

      int target_item_idx = sorted_item_indices[target_item];
      std::vector<int> to_indices;
      for (int item_idx : cluster[to])
        to_indices.push_back(sorted_item_indices[item_idx]);
      int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
      cluster[to].insert(cluster[to].begin() + pos, target_item);

      if (!is_simulation)
      {
        output_for_vis(cluster);
      }

      // insert_item_by_bs(to, target_item, cluster, suspend, q, server, is_simulation);
      // if (suspend)
      //   break;
    }

    FOR(d, 0, D)
    for (int item_idx : cluster[d])
      ans[item_idx] = d;
    return ans;
  }

  /**
   * @brief マージソートを使う。
   * 総和が大きいクラスターを見つけ、その一つの要素をソート済みの配列で隣り合うアイテムと交換する。
   *
   * @param server
   * @param is_simulation
   * @return std::vector<int>
   */
  std::vector<int> solve4(const IOServer &server, bool is_simulation)
  {
    std::vector<int> items(N);
    std::vector<int> ans(N, 0);
    std::iota(items.begin(), items.end(), 0);
    int q = 0;
    bool suspend = false;
    std::vector<int> sorted_items = merge_sort(items, q, suspend, server, is_simulation);
    std::vector<int> sorted_item_indices(N);
    FOR(i, 0, N)
    sorted_item_indices[sorted_items[i]] = i;
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::set<int> fixed_cluster_indices;
    std::vector<int> ans_tmp(N);

    FOR(d, 0, D)
    {
      for (int idx : initial_cluster[d])
      {
        cluster[d].push_back(sorted_items[idx]);
      }
    }

    if (q < Q)
      best_cluster = cluster;
    else
      return ans;

    while (q < Q)
    {
      int nl, nr;
      std::vector<int> vl, vr;
      int largest_cluster_idx;
      int smallest_cluster_idx;
      bool suspend = false;

      find_largest(
          largest_cluster_idx,
          suspend,
          cluster,
          fixed_cluster_indices,
          q,
          server,
          is_simulation);

      if (suspend)
        break;

      int from = largest_cluster_idx;
      int fsize = cluster[from].size();
      if (fsize == 1)
      {
        fixed_cluster_indices.insert(from);
        continue;
      }
      int from_idx = (u32)mt() % (fsize - 1);
      int target_item = cluster[from][from_idx];
      int from_item_idx = sorted_item_indices[target_item];
      if (from_item_idx == 0)
      {
        int to = from;
        while (to == from)
        {
          to = (u32)mt() % D;
        }
        cluster[to].insert(cluster[to].begin(), target_item);
        cluster[from].erase(cluster[from].begin());
      }
      else
      {
        int to_item_idx = from_item_idx;
        int to_item = sorted_items[to_item_idx];
        FOR(d, 0, D)
        FOR(i, 0, cluster[d].size())
        if (cluster[d][i] == to_item)
        {
          std::swap(cluster[d][i], cluster[from][from_idx]);
          break;
        }
      }
      if (!is_simulation)
      {
        output_for_vis(cluster);
      }
    }

    FOR(d, 0, D)
    for (int item_idx : cluster[d])
      ans[item_idx] = d;
    return ans;
  }

  /**
   * @brief ランダムに2つのクラスターを選ぶ。その2つを比較し、アイテムのやり取りを行う
   *
   * @param server
   * @param is_simulation
   * @return std::vector<int>
   */
  std::vector<int> solve5(const IOServer &server, bool is_simulation)
  {
    std::vector<int> ans(N, 0);
    std::vector<int> ans_tmp(N);
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), mt);
    FOR(i, 0, N)
    cluster[i % D].push_back(indices[i]);
    best_cluster = cluster;
    int q = 0;
    int nl, nr;
    std::vector<int> vl, vr;
    while (q < Q)
    {
      int pre_nl, pre_nr;
      std::vector<int> pre_vl, pre_vr;
      int u, v;

      u = (u32)mt() % D;
      v = u;
      while (u == v)
      {
        v = (u32)mt() % D;
      }

      pre_vl = cluster[u];
      pre_vr = cluster[v];
      pre_nl = pre_vl.size();
      pre_nr = pre_vr.size();
      Comp pre_res = query(pre_nl, pre_nr, pre_vl, pre_vr, q, server, is_simulation);
      if (q == Q)
      {
        break;
      }
      int u_num_choose;
      int v_num_choose;
      if (pre_res == Comp::LARGER)
      {
        if (pre_nl == 1)
          continue;
        u_num_choose = std::min((int)((u32)mt() % pre_nl + 1), pre_nl - 1);
        v_num_choose = (u32)mt() % pre_nr;
      }
      else if (pre_res == Comp::LESS)
      {
        if (pre_nr == 1)
          continue;
        u_num_choose = (u32)mt() % pre_nl;
        v_num_choose = std::min((int)((u32)mt() % pre_nr + 1), pre_nr - 1);
      }
      else
      {
        continue;
      }

      std::shuffle(cluster[u].begin(), cluster[u].end(), mt);
      std::shuffle(cluster[v].begin(), cluster[v].end(), mt);
      std::vector<int> present_from_u, present_from_v;
      FOR(i, 0, u_num_choose)
      present_from_u.emplace_back(cluster[u][i]);
      FOR(i, 0, v_num_choose)
      present_from_v.emplace_back(cluster[v][i]);

      pre_cluster = cluster;

      for (int present_item : present_from_u)
      {
        cluster[v].emplace_back(present_item);
        FOR(i, 0, cluster[u].size())
        if (cluster[u][i] == present_item)
        {
          cluster[u].erase(cluster[u].begin() + i);
          break;
        }
      }
      for (int present_item : present_from_v)
      {
        cluster[u].emplace_back(present_item);
        FOR(i, 0, cluster[v].size())
        if (cluster[v][i] == present_item)
        {
          cluster[v].erase(cluster[v].begin() + i);
          break;
        }
      }
      int nx_nl = cluster[u].size();
      int nx_nr = cluster[v].size();
      std::vector<int> nx_vl = cluster[u];
      std::vector<int> nx_vr = cluster[v];

      if (pre_res == Comp::LARGER)
      {
        std::tie(vl, vr) = delete_intersection_set(nx_vl, pre_vl);
        nl = vl.size();
        nr = vr.size();
        if (nl > 0)
        {
          Comp nx_res1 = query(nl, nr, vl, vr, q, server, is_simulation);
          if (nx_res1 == Comp::LARGER)
          {
            cluster = pre_cluster;
            continue;
          }
          if (q == Q)
            break;
        }
        std::tie(vl, vr) = delete_intersection_set(nx_vr, pre_vl);
        nl = vl.size();
        nr = vr.size();
        Comp nx_res2 = query(nl, nr, vl, vr, q, server, is_simulation);
        if (nx_res2 == Comp::LARGER)
        {
          cluster = pre_cluster;
          continue;
        }
      }
      else if (pre_res == Comp::LESS)
      {
        std::tie(vl, vr) = delete_intersection_set(nx_vl, pre_vr);
        nl = vl.size();
        nr = vr.size();
        Comp nx_res1 = query(nl, nr, vl, vr, q, server, is_simulation);
        if (nx_res1 == Comp::LARGER)
        {
          cluster = pre_cluster;
          continue;
        }
        std::tie(vl, vr) = delete_intersection_set(nx_vr, pre_vr);
        nl = vl.size();
        nr = vr.size();
        if (nl > 0)
        {
          if (q == Q)
            break;
          Comp nx_res2 = query(nl, nr, vl, vr, q, server, is_simulation);
          if (nx_res2 == Comp::LARGER)
          {
            cluster = pre_cluster;
            continue;
          }
        }
      }
      best_cluster = cluster;
      if (!is_simulation)
      {
        output_for_vis(best_cluster);
      }
    } // while(q < Q)

    for (int d = 0; d < D; d++)
    {
      for (int idx : best_cluster[d])
        ans[idx] = d;
    }
    return ans;
  }

  /**
   * @brief
   * 1. 最大or最小のターゲットクラスターを見つける
   * 2. それ以外のクラスターをランダムに選ぶ。
   * 3. その2つのクラスター間でアイテムのやり取り
   * 4. やり取り前のターゲットクラスターと値を比較してオーバーしたらundo
   *
   * @param server
   * @param is_simulation
   * @return std::vector<int>
   */
  std::vector<int> solve6(IOServer &server, bool is_simulation)
  {
    std::vector<int> ans(N, 0);
    std::vector<int> ans_tmp(N);
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::set<int> fixed_cluster_indices;
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), mt);
    FOR(i, 0, N)
    cluster[i % D].push_back(indices[i]);
    best_cluster = cluster;
    int q = 0;
    int nl, nr;
    std::vector<int> vl, vr;
    while (q < Q)
    {
      bool large;
      int largest_index;
      int smallest_index;
      bool suspend = false;
      if ((u32)mt() % 2 == 0)
      {
        large = true;
        find_largest(
            largest_index,
            suspend,
            cluster,
            fixed_cluster_indices,
            q,
            server,
            is_simulation);

        if (cluster[largest_index].size() == 1)
        {
          fixed_cluster_indices.insert(largest_index);
          continue;
        }
      }
      else
      {
        large = false;
        find_smallest(
            smallest_index,
            suspend,
            cluster,
            fixed_cluster_indices,
            q,
            server,
            is_simulation);
      }

      if (!suspend)
      {
        int max_iter = 5;
        bool done = false;
        while (max_iter-- && !done)
        {
          int pat = (u32)mt() % 2;
          if (pat == 0)
          {
            std::vector<int> indices_D(D);
            std::iota(indices_D.begin(), indices_D.end(), 0);
            std::shuffle(indices_D.begin(), indices_D.end(), mt);
            int from;
            int to;
            if (large)
            {
              from = largest_index;
              for (int idx : indices_D)
              {
                if (fixed_cluster_indices.find(idx) != fixed_cluster_indices.end() || idx == from)
                {
                  continue;
                }
                to = idx;
                break;
              }
            }
            else
            {
              to = smallest_index;
              for (int idx : indices_D)
              {
                if (fixed_cluster_indices.find(idx) != fixed_cluster_indices.end() || idx == to)
                {
                  continue;
                }
                if (cluster[idx].size() == 1)
                {
                  continue;
                }
                from = idx;
                break;
              }
            }

            try_move_one_item(
                from,
                to,
                cluster, best_cluster,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                large,
                done);
          }
          else if (pat == 1)
          {
            int from;
            if (large)
            {
              from = largest_index;
            }
            else
            {
              from = smallest_index;
            }

            try_exchange_one_item(
                from,
                cluster,
                best_cluster,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                large,
                done);
          }

          if (q == Q)
            break;
        } // while(max_iter--)
      }
    }

    FOR(d, 0, D)
    for (int item_idx : best_cluster[d])
      ans[item_idx] = d;

    return ans;
  }

  /**
   * @brief merge sort ベースの解法
   * solve6に続いてアイテムの移動、アイテムの交換で状態の遷移をする
   *
   *
   * @param server
   * @param is_simulation
   * @return std::vector<int>
   */
  std::vector<int> solve7(IOServer &server, bool is_simulation)
  {
    std::vector<int> items(N);
    std::vector<int> ans(N, 0);
    FOR(i, 0, N)
    ans[i] = i % D;
    std::iota(items.begin(), items.end(), 0);
    int q = 0;
    bool suspend = false;
    std::vector<int> sorted_items = merge_sort(items, q, suspend, server, is_simulation);
    std::vector<int> sorted_item_indices(N);
    FOR(i, 0, N)
    sorted_item_indices[sorted_items[i]] = i;
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::set<int> fixed_cluster_indices;
    std::vector<int> ans_tmp(N);

    FOR(d, 0, D)
    {
      for (int idx : initial_cluster[d])
      {
        cluster[d].push_back(sorted_items[idx]);
      }
    }

    if (q < Q)
      best_cluster = cluster;
    else
      return ans;

    while (q < Q)
    {
      // FOR(d, 0, D) {
      //   std::cerr << "cluster " << d << ": ";
      //   for(int item_idx: cluster[d]) std::cerr << item_idx << " ";
      //   std::cerr << "\n";
      // }
      bool large;
      int largest_index;
      int smallest_index;
      bool suspend = false;
      if ((u32)mt() % 2 == 0)
      {
        large = true;
        find_largest(
            largest_index,
            suspend,
            cluster,
            fixed_cluster_indices,
            q,
            server,
            is_simulation);

        if (cluster[largest_index].size() == 1)
        {
          fixed_cluster_indices.insert(largest_index);
          continue;
        }
      }
      else
      {
        large = false;
        find_smallest(
            smallest_index,
            suspend,
            cluster,
            fixed_cluster_indices,
            q,
            server,
            is_simulation);
      }
      // >---- 最大最小 ----<

      if (!suspend)
      {
        bool done = false;
        int max_iter = 5;
        std::vector<int> indices_D(D);
        std::iota(indices_D.begin(), indices_D.end(), 0);
        std::shuffle(indices_D.begin(), indices_D.end(), mt);

        int from;
        int to;
        int nl, nr;
        std::vector<int> vl, vr;

        while (max_iter-- && !done)
        {
          int pat = (u32)mt() % 10;
          if (q == Q)
          {
            break;
          }
          if (pat < -1)
          { // exchange pattern
            int from;
            if (large)
            {
              from = largest_index;
            }
            else
            {
              from = smallest_index;
            }

            try_exchange_one_item(
                from,
                cluster,
                best_cluster,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                large,
                done);

            if (done)
            {
              best_cluster = cluster;
            }
          }
          else if (pat < 5)
          { // move pattern
            if (large)
            {
              from = largest_index;
              to = -1;
              for (int idxD : indices_D)
              {
                if (idxD == from || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
                {
                  continue;
                }
                to = idxD;
                break;
              }
            }
            else
            {
              from = -1;
              to = smallest_index;
              for (int idxD : indices_D)
              {
                if (idxD == to || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
                {
                  continue;
                }
                if (cluster[idxD].size() == 1)
                {
                  continue;
                }
                from = idxD;
                break;
              }
            }
            if (from == -1 || to == -1)
            {
              continue;
            }

            int target_item = find_move_item_by_bs(from, to, suspend, cluster, q, server, is_simulation);
            if (q == Q)
            {
              break;
            }
            if (target_item == -1)
            {
              continue;
            }
            pre_cluster = cluster;

            auto it = std::find(cluster[from].begin(), cluster[from].end(), target_item);
            cluster[from].erase(it);

            int target_item_idx = sorted_item_indices[target_item];
            std::vector<int> to_indices;
            for (int item_idx : cluster[to])
              to_indices.push_back(sorted_item_indices[item_idx]);
            int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
            cluster[to].insert(cluster[to].begin() + pos, target_item);

            // if (large)
            // {
            //   vl = cluster[from];
            //   vr = pre_cluster[to];
            // }
            // else
            // {
            //   vl = cluster[to];
            //   vr = pre_cluster[from];
            // }
            // std::tie(vl, vr) = delete_intersection_set(vl, vr);
            // nl = vl.size();
            // nr = vr.size();
            // Comp res = query(nl, nr, vl, vr, q, server, is_simulation);
            // if (large && res == Comp::LESS)
            // {
            //   cluster = pre_cluster;
            //   continue;
            // }
            // else if (!large && res == Comp::LARGER)
            // {
            //   cluster = pre_cluster;
            //   continue;
            // }
            best_cluster = cluster;
            done = true;
          }
          else
          { // other pattern
            int target;
            int opponent = -1;
            if (large)
            {
              target = largest_index;
            }
            else
            {
              target = smallest_index;
            }
            try_exchange_one_sorted_item(
                target,
                opponent,
                cluster,
                best_cluster,
                sorted_items,
                sorted_item_indices,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                large,
                done);
            if (done)
            {
              best_cluster = cluster;
            }
          }
        } // while(max_iter--)
        if (!is_simulation)
        {
          output_for_vis(best_cluster);
        }
      }
    }

    FOR(d, 0, D)
    for (int item_idx : best_cluster[d])
      ans[item_idx] = d;

    return ans;
  }

  /**
   * @brief tolerate機能を追加
   * 
   * @param server 
   * @param is_simulation 
   */
  std::vector<int> solve8(IOServer& server, bool is_simulation) {
    std::vector<int> ans(N);
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster(D);
    std::vector<std::vector<int>> pre_cluster;
    std::vector<int> max_cluster, min_cluster;
    std::set<int> fixed_cluster_indices;
    int max_over_cluster_idx, min_over_cluster_idx;
    int q = 0;
    bool suspend = false;
    enum class Status {
      NORMAL,
      MAX_OVER, 
      MIN_OVER,
    };
    Status cur_status = Status::NORMAL;

    // ansの適当な初期化
    FOR(i, 0, N) ans[i] = i % D;

    // merge sortとアイテムのインデックスのソート
    std::vector<int> items(N);
    std::iota(items.begin(), items.end(), 0);
    std::vector<int> sorted_items = merge_sort(items, q, suspend, server, is_simulation);
    std::vector<int> sorted_item_indices(N);
    FOR(i, 0, N) sorted_item_indices[sorted_items[i]] = i;

    FOR(d, 0, D)
    {
      for (int idx : initial_cluster[d])
      {
        cluster[d].push_back(sorted_items[idx]);
      }
    }

    if(suspend) {
      return ans;
    } else {
      best_cluster = cluster;
      if(!is_simulation) {
        output_for_vis(best_cluster);
      }
    }

    while(q < Q) {
      int num_tolerate = 10;
      int largest_index;
      int smallest_index;
      int nl, nr;
      std::vector<int> vl, vr;
      if(!is_simulation) {
        output_for_vis(best_cluster);
      }
      if(cur_status == Status::NORMAL) {
        bool done = false;
        bool target_is_large = (u32)mt() % 2 == 0;
        int from, to;
        if(target_is_large) {
          find_largest(
              largest_index,
              suspend,
              cluster,
              fixed_cluster_indices,
              q,
              server,
              is_simulation);
          if(cluster[largest_index].size() == 1) {
            fixed_cluster_indices.insert(largest_index);
            continue;
          }
        } else {
          find_smallest(
              smallest_index,
              suspend,
              cluster,
              fixed_cluster_indices,
              q,
              server,
              is_simulation);
        }

        if(suspend) {
          break;
        }

        int max_iter = 5;
        std::vector<int> indices_D(D);
        std::iota(indices_D.begin(), indices_D.end(), 0);

        pre_cluster = cluster;

        while(max_iter-- && !done) {
          int pat = (u32)mt() % 50;
          std::shuffle(indices_D.begin(), indices_D.end(), mt);
          if(q == Q)
          {
            break;
          }

          if(pat == 0)
          { // move random item to another cluster
            int idx_tmp = -1;
            int idx_tmp2 = -1;
            if(target_is_large) {
              idx_tmp = largest_index;
            } else {
              idx_tmp = smallest_index;
            }
            for(int idxD: indices_D) {
              if(idxD == idx_tmp || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end()) {
                continue;
              }
              if(!target_is_large && cluster[idxD].size() == 1) {
                continue;
              }
              idx_tmp2 = idxD;
              break;
            }
            if(idx_tmp == -1 || idx_tmp2 == -1) {
              continue;
            }

            if(target_is_large) {
              from = idx_tmp;
              to = idx_tmp2;
            } else {
              from = idx_tmp2;
              to = idx_tmp;
            }

            int idx = (u32)mt() % cluster[from].size();
            int target_item = cluster[from][idx];
            int target_item_idx = sorted_item_indices[target_item];
            std::vector<int> to_indices;
            for(int item_idx: cluster[to]) to_indices.emplace_back(sorted_item_indices[item_idx]);
            cluster[from].erase(cluster[from].begin() + idx);
            int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
            cluster[to].insert(cluster[to].begin() + pos, target_item);
            if(target_is_large) {
              vl = cluster[to];
              vr = pre_cluster[from];
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size(); nr = vr.size();
              assert(nl > 0 && nr > 0);
              Comp res = query(nl, nr, vl, vr, q, server, is_simulation);

              if(res == Comp::LARGER) {
                cur_status = Status::MAX_OVER;
                max_cluster = pre_cluster[from];
                max_over_cluster_idx = to;
                done = true;
                break;
              }
            } else {
              vl = cluster[from];
              vr = pre_cluster[to];
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size(); nr = vr.size();
              assert(nl > 0 && nr > 0);
              Comp res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(res == Comp::LESS) {
                cur_status = Status::MIN_OVER;
                min_cluster = pre_cluster[to];
                min_over_cluster_idx = from;
                done = true;
                break;
              }
            }
            best_cluster = cluster;
            done = true;
          } 
          else if(pat == 1)
          { // exchange ramdom items
            // >-- from to の設定 --<
            int idx_tmp;
            int idx_tmp2;
            if(target_is_large) {
              idx_tmp = largest_index;
            } else {
              idx_tmp = smallest_index;
            }
            idx_tmp2 = -1;
            for(int idxD: indices_D) {
              if(idxD == idx_tmp || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end()) {
                continue;
              }
              idx_tmp2 = idxD;
              break;
            }
            if(idx_tmp2 == -1) {
              continue;
            }

            if(target_is_large) {
              from = idx_tmp;
              to = idx_tmp2;
            } else {
              from = idx_tmp2;
              to = idx_tmp;
            }
            // >-- from toの設定終了 --<
            int u = from;
            int v = to;
            int usize = cluster[from].size();
            int vsize = cluster[to].size();
            int uidx = (u32)mt() % usize;
            int vidx = (u32)mt() % vsize;
            int uitem = cluster[u][uidx];
            int vitem = cluster[v][vidx];
            int uitem_idx = sorted_item_indices[uitem];
            int vitem_idx = sorted_item_indices[vitem];
            std::vector<int> u_indices, v_indices;
            cluster[u].erase(cluster[u].begin() + uidx);
            cluster[v].erase(cluster[v].begin() + vidx);
            for(int item_idx: cluster[from]) u_indices.emplace_back(sorted_item_indices[item_idx]);
            for(int item_idx: cluster[to]) v_indices.emplace_back(sorted_item_indices[item_idx]);
            auto uit = std::lower_bound(u_indices.begin(), u_indices.end(), vitem_idx);
            auto vit = std::lower_bound(v_indices.begin(), v_indices.end(), uitem_idx);
            u_indices.insert(uit, vitem_idx);
            v_indices.insert(vit, uitem_idx);
            std::vector<int> nx_u_cluster, nx_v_cluster;
            for(int idx: u_indices) nx_u_cluster.emplace_back(sorted_items[idx]);
            for(int idx: v_indices) nx_v_cluster.emplace_back(sorted_items[idx]);
            cluster[from] = nx_u_cluster;
            cluster[to] = nx_v_cluster;

            if(target_is_large)
            {
              vl = cluster[from];
              vr = pre_cluster[from];
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size();
              nr = vr.size();
              if(nl > 0 && nr > 0) {
                Comp res = query(nl, nr, vl, vr, q, server, is_simulation); 
                if(res == Comp::LARGER) {
                  cur_status = Status::MAX_OVER;
                  max_cluster = pre_cluster[from];
                  max_over_cluster_idx = from;
                  done = true;
                  break;
                }
              }
              vl = cluster[to];
              vr = pre_cluster[from];
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size();
              nr = vr.size();
              if(nl > 0 && nr > 0) {
                Comp res = query(nl, nr, vl, vr, q, server, is_simulation); 
                if(res == Comp::LARGER) {
                  cur_status = Status::MAX_OVER;
                  max_cluster = pre_cluster[to];
                  max_over_cluster_idx = to;
                  done = true;
                  break;
                }
              }

            }
            else
            {// for smallest cluster
              vl = cluster[from];
              vr = pre_cluster[to]; // ref
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size();
              nr = vr.size();
              if(nl > 0 && nr > 0) {
                Comp res = query(nl, nr, vl, vr, q, server, is_simulation); 
                if(res == Comp::LESS) {
                  cur_status = Status::MIN_OVER;
                  min_cluster = pre_cluster[from];
                  min_over_cluster_idx = from;
                  done = true;
                  break;
                }
              }
              vl = cluster[to];
              vr = pre_cluster[to]; // ref
              std::tie(vl, vr) = delete_intersection_set(vl, vr);
              nl = vl.size();
              nr = vr.size();
              if(nl > 0 && nr > 0) {
                Comp res = query(nl, nr, vl, vr, q, server, is_simulation); 
                if(res == Comp::LESS) {
                  cur_status = Status::MIN_OVER;
                  min_cluster = pre_cluster[to];
                  min_over_cluster_idx = to;
                  done = true;
                  break;
                }
              }
            }

            best_cluster = cluster;
            done = true;
            
          }
          else if(pat < 26) 
          { // move valid item to another
            if (target_is_large)
            {
              from = largest_index;
              to = -1;
              for (int idxD : indices_D)
              {
                if (idxD == from || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
                {
                  continue;
                }
                to = idxD;
                break;
              }
            }
            else
            {
              from = -1;
              to = smallest_index;
              for (int idxD : indices_D)
              {
                if (idxD == to || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
                {
                  continue;
                }
                if (cluster[idxD].size() == 1)
                {
                  continue;
                }
                from = idxD;
                break;
              }
            }
            if (from == -1 || to == -1)
            {
              continue;
            }

            int target_item = find_move_item_by_bs(from, to, suspend, cluster, q, server, is_simulation);
            if (q == Q)
            {
              break;
            }
            if (target_item == -1)
            {
              continue;
            }
            pre_cluster = cluster;

            auto it = std::find(cluster[from].begin(), cluster[from].end(), target_item);
            cluster[from].erase(it);

            int target_item_idx = sorted_item_indices[target_item];
            std::vector<int> to_indices;
            for (int item_idx : cluster[to]) {
              to_indices.push_back(sorted_item_indices[item_idx]);
            }
            int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
            cluster[to].insert(cluster[to].begin() + pos, target_item);

            best_cluster = cluster;
            done = true;

          } // move valid one item
          else 
          { // exchange valid one item
            int target;
            int opponent = -1;
            if (target_is_large)
            {
              target = largest_index;
            }
            else
            {
              target = smallest_index;
            }
            bool success = false;
            try_exchange_one_sorted_item(
                target,
                opponent,
                cluster,
                best_cluster,
                sorted_items,
                sorted_item_indices,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                target_is_large,
                success);
            if (success)
            {
              done = true;
              best_cluster = cluster;
            }

          }
        }

      }
      else if(cur_status == Status::MAX_OVER)
      {
        bool target_is_large = true;
        int from = max_over_cluster_idx;
        int to;
        bool done = false;
        std::vector<int> indices_D(D);
        std::iota(indices_D.begin(), indices_D.end(), 0);
        while(num_tolerate-- && !done) {
          if(q == Q) {
            break;
          }
          find_largest(
            max_over_cluster_idx, 
            suspend, 
            cluster, 
            fixed_cluster_indices, 
            q, 
            server, 
            is_simulation);
          if(suspend) {
            break;
          }
          int pat = (u32)mt() % 10;
          std::shuffle(indices_D.begin(), indices_D.end(), mt);

          if(pat < 5) 
          { // move valid item to another
            from = max_over_cluster_idx;
            to = -1;
            for (int idxD : indices_D)
            {
              if (idxD == from || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
              {
                continue;
              }
              to = idxD;
              break;
            }
            if (from == -1 || to == -1)
            {
              continue;
            }

            int target_item = find_move_item_by_bs(from, to, suspend, cluster, q, server, is_simulation);
            if (q == Q)
            {
              break;
            }
            if (target_item == -1)
            {
              continue;
            }
            pre_cluster = cluster;

            auto it = std::find(cluster[from].begin(), cluster[from].end(), target_item);
            cluster[from].erase(it);

            int target_item_idx = sorted_item_indices[target_item];
            std::vector<int> to_indices;
            for (int item_idx : cluster[to]) {
              to_indices.push_back(sorted_item_indices[item_idx]);
            }
            int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
            cluster[to].insert(cluster[to].begin() + pos, target_item);

            vl = cluster[from];
            vr = max_cluster;  // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LESS) {
                best_cluster = cluster;
                done = true;
              }
            }
          } // move valid one item
          else 
          { // exchange valid one item
            int target = max_over_cluster_idx;
            int opponent = -1;
            bool success = false;
            try_exchange_one_sorted_item(
                target,
                opponent,
                cluster,
                best_cluster,
                sorted_items,
                sorted_item_indices,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                target_is_large,
                success);

            if(q == Q) {
              break;
            }

            if(!success) {
              continue;
            }

            bool res1 = true;
            bool res2 = true;

            vl = cluster[target];
            vr = max_cluster; // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LARGER) {
                res1 = false;
              }
            }
            vl = cluster[opponent];
            vr = max_cluster; // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LARGER) {
                res2 = false;
              }
            }

            if(res1 && res2) {
              best_cluster = cluster;
              done = true;
            }
          }

          cur_status = Status::NORMAL;
        } // while(num_tolerate--)

        if(!done) {
          cluster = best_cluster; // undo
        }
      }
      else if(cur_status == Status::MIN_OVER)
      {
        int to = min_over_cluster_idx;
        int from;
        bool target_is_large = false;
        bool done = false;
        std::vector<int> indices_D(D);
        std::iota(indices_D.begin(), indices_D.end(), 0);
        while(num_tolerate-- && !done) {
          if(q == Q) {
            break;
          }
          find_smallest(
            min_over_cluster_idx, 
            suspend, 
            cluster, 
            fixed_cluster_indices, 
            q, 
            server, 
            is_simulation);
          if(suspend) {
            break;
          }

          int pat = (u32)mt() % 10;
          std::shuffle(indices_D.begin(), indices_D.end(), mt);

          if(pat < 5) 
          { // move valid item to another
            to = min_over_cluster_idx;
            from = to;
            for (int idxD : indices_D)
            {
              if (idxD == to || fixed_cluster_indices.find(idxD) != fixed_cluster_indices.end())
              {
                continue;
              }
              if(cluster[idxD].size() == 1) {
                continue;
              }
              from = idxD;
              break;
            }
            if (from == -1 || to == -1)
            {
              continue;
            }

            int target_item = find_move_item_by_bs(from, to, suspend, cluster, q, server, is_simulation);
            if (q == Q)
            {
              break;
            }
            if (target_item == -1)
            {
              continue;
            }

            auto it = std::find(cluster[from].begin(), cluster[from].end(), target_item);
            cluster[from].erase(it);

            int target_item_idx = sorted_item_indices[target_item];
            std::vector<int> to_indices;
            for (int item_idx : cluster[to]) {
              to_indices.push_back(sorted_item_indices[item_idx]);
            }
            int pos = std::lower_bound(to_indices.begin(), to_indices.end(), target_item_idx) - to_indices.begin();
            cluster[to].insert(cluster[to].begin() + pos, target_item);

            vl = cluster[to];
            vr = min_cluster;  // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LARGER) {
                best_cluster = cluster;
                done = true;
              }
            }
          } // move valid one item
          else 
          { // exchange valid one item
            int target = min_over_cluster_idx;
            int opponent = -1;
            bool success = false;
            try_exchange_one_sorted_item(
                target,
                opponent,
                cluster,
                best_cluster,
                sorted_items,
                sorted_item_indices,
                fixed_cluster_indices,
                q,
                server,
                is_simulation,
                target_is_large,
                success);

            if(q == Q) {
              break;
            }

            if(!success) {
              continue;
            }

            bool res1 = true;
            bool res2 = true;
            vl = cluster[target];
            vr = min_cluster; // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LESS) {
                res1 = false;
              }
            }
            if(q == Q) {
              break;
            }
            vl = cluster[opponent];
            vr = min_cluster; // ref
            std::tie(vl, vr) = delete_intersection_set(vl, vr);
            nl = vl.size();
            nr = vr.size();
            if(nl > 0 && nr > 0) {
              Comp valid_res = query(nl, nr, vl, vr, q, server, is_simulation);
              if(valid_res == Comp::LESS) {
                res2 = false;
              }
            }

            if(res1 && res2) {
              best_cluster = cluster;
              done = true;
            }
          }
        } // while()
        if(!done) {
          cluster = best_cluster; // undo
        }
        cur_status = Status::NORMAL;
      }

    }

    FOR(d, 0, D) for(int item_idx: best_cluster[d]) ans[item_idx] = d;

    return ans;
  }

  void try_move_one_item(
      int from,
      int to,
      std::vector<std::vector<int>> &cluster,
      std::vector<std::vector<int>> &best_cluster,
      std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation,
      bool large,
      bool &done)
  {

    std::vector<std::vector<int>> pre_cluster = cluster;
    int fsize = cluster[from].size();
    int nl, nr;
    std::vector<int> vl, vr;

    int idx = (u32)mt() % fsize;
    int target_item = cluster[from][idx];

    cluster[to].emplace_back(target_item);
    FOR(i, 0, cluster[from].size())
    if (cluster[from][i] == target_item)
    {
      cluster[from].erase(cluster[from].begin() + i);
      break;
    }

    if (large)
    {
      vl = pre_cluster[from];
      vr = cluster[to];
    }
    else
    {
      vl = cluster[from];
      vr = pre_cluster[to];
    }

    std::tie(vl, vr) = delete_intersection_set(vl, vr);
    nl = vl.size();
    nr = vr.size();

    Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
    if (res == Comp::LESS)
    {
      cluster = pre_cluster;
    }
    else
    {
      best_cluster = cluster;
      done = true;
      if (!is_simulation)
      {
        output_for_vis(best_cluster);
      }
    }
  }

  void try_exchange_one_item(
      int from,
      std::vector<std::vector<int>> &cluster,
      std::vector<std::vector<int>> &best_cluster,
      std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation,
      bool large,
      bool &done)
  {
    std::vector<std::vector<int>> pre_cluster = cluster;
    std::vector<int> indices(D);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), mt);
    int u = from;
    int v;
    int nl, nr;
    std::vector<int> vl, vr;
    for (int idx : indices)
    {
      if (fixed.find(idx) != fixed.end() || idx == u)
      {
        continue;
      }
      v = idx;
      break;
    }

    int usize = cluster[u].size();
    int vsize = cluster[v].size();
    int uidx = (u32)mt() % usize;
    int vidx = (u32)mt() % vsize;
    int uitem = cluster[u][uidx];
    int vitem = cluster[v][vidx];
    std::swap(cluster[u][uidx], cluster[v][vidx]);

    vl = pre_cluster[u];
    vr = cluster[u];
    std::tie(vl, vr) = delete_intersection_set(vl, vr);
    nl = vl.size();
    nr = vr.size();
    if (nl > 0 && nr > 0)
    {
      Comp res1 = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res1 == Comp::LESS && large)
      {
        cluster = pre_cluster;
        return;
      }
      else if (res1 == Comp::LARGER && !large)
      {
        cluster = pre_cluster;
        return;
      }
    }

    vl = pre_cluster[u];
    vr = cluster[v];
    std::tie(vl, vr) = delete_intersection_set(vl, vr);
    nl = vl.size();
    nr = vr.size();
    if (nl > 0 && nr > 0)
    {
      if (num_query == Q)
      {
        return;
      }
      Comp res2 = query(nl, nr, vl, vr, num_query, server, is_simulation);
      if (res2 == Comp::LESS && large)
      {
        cluster = pre_cluster;
        return;
      }
      else if (res2 == Comp::LARGER && !large)
      {
        cluster = pre_cluster;
        return;
      }
    }

    best_cluster = cluster;
    done = true;
    if (!is_simulation)
    {
      output_for_vis(best_cluster);
    }
  }

  void try_exchange_one_sorted_item(
      int target_cluster_idx,
      int& opponent_cluster_idx,
      std::vector<std::vector<int>> &cluster,
      std::vector<std::vector<int>> &best_cluster,
      const std::vector<int> &sorted_items,
      const std::vector<int> &sorted_item_indices,
      std::set<int> &fixed,
      int &num_query,
      const IOServer &server,
      bool is_simulation,
      bool large,
      bool &done)
  {
    constexpr int NUM_RANGE = 3;
    std::vector<std::vector<int>> pre_cluster = cluster;
    int tsize = cluster[target_cluster_idx].size();
    int target_idx = (u32)mt() % tsize;
    int target_item = cluster[target_cluster_idx][target_idx];
    int target_item_idx = sorted_item_indices[target_item];
    int nl, nr;
    std::vector<int> vl, vr;
    auto swap_two_item = [&](int u, int v, int uitem, int vitem)
    {
      std::vector<int> u_item_indices;
      std::vector<int> v_item_indices;
      int u_item_idx = sorted_item_indices[uitem];
      int v_item_idx = sorted_item_indices[vitem];
      for (int item_idx : cluster[u])
        u_item_indices.emplace_back(sorted_item_indices[item_idx]);
      for (int item_idx : cluster[v])
        v_item_indices.emplace_back(sorted_item_indices[item_idx]);
      auto it1 = std::lower_bound(u_item_indices.begin(), u_item_indices.end(), v_item_idx);
      auto it2 = std::lower_bound(v_item_indices.begin(), v_item_indices.end(), u_item_idx);
      u_item_indices.insert(it1, v_item_idx);
      v_item_indices.insert(it2, u_item_idx);
      auto uit = std::find(u_item_indices.begin(), u_item_indices.end(), u_item_idx);
      auto vit = std::find(v_item_indices.begin(), v_item_indices.end(), v_item_idx);
      u_item_indices.erase(uit);
      v_item_indices.erase(vit);
      std::vector<int> nx_u_cluster, nx_v_cluster;
      for (int idx : u_item_indices)
        nx_u_cluster.emplace_back(sorted_items[idx]);
      for (int idx : v_item_indices)
        nx_v_cluster.emplace_back(sorted_items[idx]);
      cluster[u] = nx_u_cluster;
      cluster[v] = nx_v_cluster;
    };

    if (large)
    {
      if (target_item_idx == 0)
      { // first index
        int from = target_cluster_idx;
        int to = target_cluster_idx;
        while (to == target_cluster_idx)
        {
          to = (u32)mt() % D;
        }
        cluster[to].insert(cluster[to].begin(), target_item);
        cluster[from].erase(cluster[from].begin());

        vl = cluster[to];
        vr = pre_cluster[from];
        std::tie(vl, vr) = delete_intersection_set(vl, vr);
        nl = vl.size();
        nr = vr.size();

        Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
        if (res == Comp::LARGER)
        {
          cluster = pre_cluster;
          return;
        }
        opponent_cluster_idx = to;
        done = true;
      }
      else
      { // not first
        int diff = (u32)mt() % NUM_RANGE + 1;
        int to_item_idx = std::max(0, target_item_idx - diff);
        int to_item = sorted_items[to_item_idx];
        int to_cluster_idx = -1;
        FOR(d, 0, D)
        {
          if (to_cluster_idx != -1)
          {
            break;
          }
          FOR(i, 0, cluster[d].size())
          if (cluster[d][i] == to_item)
          {
            if (d == target_cluster_idx)
              return;
            to_cluster_idx = d;
            break;
          }
        }
        swap_two_item(target_cluster_idx, to_cluster_idx, target_item, to_item);

        vl = cluster[to_cluster_idx];
        vr = pre_cluster[target_cluster_idx];
        std::tie(vl, vr) = delete_intersection_set(vl, vr);
        nl = vl.size();
        nr = vr.size();
        if (nl > 0 && nr > 0)
        {
          Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
          if (res == Comp::LARGER)
          {
            cluster = pre_cluster;
            return;
          }
          opponent_cluster_idx = to_cluster_idx;
          done = true;
        }
      }
    }
    else
    { // for smallest index
      if (target_item_idx == N - 1)
      {
        return;
      }
      else
      {
        int diff = (u32)mt() % NUM_RANGE + 1;
        int to_item_idx = std::min(N - 1, target_item_idx + diff);
        int to_item = sorted_items[to_item_idx];
        int to_cluster_idx = -1;
        FOR(d, 0, D)
        {
          if (to_cluster_idx != -1)
          {
            break;
          }
          FOR(i, 0, cluster[d].size())
          {
            if (cluster[d][i] == to_item)
            {
              if (d == target_cluster_idx)
                return;
              to_cluster_idx = d;
              break;
            }
          }
        }
        swap_two_item(target_cluster_idx, to_cluster_idx, target_item, to_item);

        vl = cluster[to_cluster_idx];
        vr = pre_cluster[target_cluster_idx];
        std::tie(vl, vr) = delete_intersection_set(vl, vr);
        nl = vl.size();
        nr = vr.size();
        if (nl > 0 && nr > 0)
        {
          Comp res = query(nl, nr, vl, vr, num_query, server, is_simulation);
          if (res == Comp::LESS)
          {
            cluster = pre_cluster;
            return;
          }
          opponent_cluster_idx = to_cluster_idx;
          done = true;
        }
      }
    }
  }

  std::vector<int> merge_sort(
      const std::vector<int> &items,
      int &num_query,
      bool &suspend,
      const IOServer &server,
      bool simulation)
  {
    std::vector<int> ret = items;
    int num_item = items.size();
    std::vector<int> buf(ret.size(), 0);
    int mid = num_item / 2;
    sub_merge_sort(0, mid, ret, buf, num_query, suspend, server, simulation);
    sub_merge_sort(mid, num_item, ret, buf, num_query, suspend, server, simulation);
    merge(0, mid, num_item, ret, buf, num_query, suspend, server, simulation);
    return ret;
  }

  void sub_merge_sort(
      int l,
      int r,
      std::vector<int> &arr,
      std::vector<int> &buf,
      int &num_query,
      bool &suspend,
      const IOServer &server,
      bool simulation)
  {
    if (num_query == Q)
    {
      suspend = true;
      return;
    }
    if (l + 1 == r)
    {
      return;
    }
    int mid = (l + r) / 2;
    sub_merge_sort(l, mid, arr, buf, num_query, suspend, server, simulation);
    sub_merge_sort(mid, r, arr, buf, num_query, suspend, server, simulation);
    merge(l, mid, r, arr, buf, num_query, suspend, server, simulation);
  }

  void merge(
      int l,
      int m,
      int r,
      std::vector<int> &arr,
      std::vector<int> &buf,
      int &num_query,
      bool &suspend,
      const IOServer &server,
      bool simulation)
  {
    if (num_query == Q)
    {
      return;
    }
    int i1 = l;
    int i2 = m;
    for (int idx = l; idx < r; idx++)
    {
      if (i1 == m)
      {
        buf[idx] = arr[i2++];
        continue;
      }
      else if (i2 == r)
      {
        buf[idx] = arr[i1++];
        continue;
      }
      else
      {
        Comp res = query(1, 1, {arr[i1]}, {arr[i2]}, num_query, server, simulation);
        if (res == Comp::LESS)
        {
          buf[idx] = arr[i1++];
        }
        else
        {
          buf[idx] = arr[i2++];
        }

        if (num_query == Q)
        {
          suspend = true;
          return;
        }
      }
    }
    for (int idx = l; idx < r; idx++)
      arr[idx] = buf[idx];
  }

  std::pair<std::vector<int>, std::vector<int>> delete_intersection_set(const std::vector<int> &vx, const std::vector<int> &vy)
  {
    std::vector<int> cnt(N, 0);
    for (int item_idx : vx)
      cnt[item_idx]++;
    for (int item_idx : vy)
      cnt[item_idx]++;
    std::vector<int> retx, rety;
    for (int item_idx : vx)
      if (cnt[item_idx] == 1)
        retx.emplace_back(item_idx);
    for (int item_idx : vy)
      if (cnt[item_idx] == 1)
        rety.emplace_back(item_idx);
    return {retx, rety};
  }

  void output_for_vis(const std::vector<std::vector<int>> &cluster)
  {
    std::vector<int> ans(N);
    FOR(d, 0, D)
    for (int item_idx : cluster[d])
      ans[item_idx] = d;
    std::cout << "#c ";
    FOR(i, 0, N)
    std::cout << ans[i] << " \n"[i + 1 == N];
  }
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