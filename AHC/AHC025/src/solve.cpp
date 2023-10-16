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

  Solver()
  {
  } // constructor

  void init()
  {
    std::cin >> N >> D >> Q;

#ifdef _LOCAL
    std::vector<int> w(N);
    FOR(i, 0, N)
    std::cin >> w[i];
    server = IOServer(N, D, Q, w);
#endif

    fprintf(stderr, "[case] N: %d, D: %d, Q: %d \n", N, D, Q);
    std::cerr << "Initialization done\n";
  }

  //
  void solve()
  {
    const int NUM_SIMULATION = 20;
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
    default:
      std::cerr << "Invalid result\n";
      break;
    }

    for (int i = 0; i < N; i++)
      std::cout << ans[i] << " \n"[i + 1 == N];
#ifdef _LOCAL
    std::cerr << "Score: " << server.calc_score(ans) << "\n";
#endif

  } // solve

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

  void merge_sort_test()
  {
    std::vector<std::vector<int>> cluster(D);
    FOR(i, 0, N)
    cluster[i % D].push_back(i);
    int q = 0;
    bool suspend = false;
    while (q < Q)
    {
      if (q == 0)
        cluster[0] = this->merge_sort(cluster[0], q, suspend, false);
      else
        auto tmp = query(1, 1, {0}, {1}, q, server, false);
      if (!suspend)
      {
      }
    }
    std::vector<int> ans(N, 0);
    FOR(d, 0, D)
    for (int item_idx : cluster[d]) ans[item_idx] = d;
    FOR(i, 0, N)
        std::cout << ans[i] << " \n"[i + 1 == N];
    for (int item_idx : cluster[0])
    {
      std::cerr << item_idx << ": " << server.w[item_idx] << std::endl;
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
  int simulation(int num_simulation = 10)
  {
    int num_solvers = 2;
    std::vector<ll> scores(num_solvers, 0);
    for (int sid = 0; sid < num_simulation; sid++)
    {
      std::cerr << "simulation iter: " << sid + 1 << std::endl;
      std::vector<int> W = generate_virtual_weights(N, D, sid);
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
    std::vector<std::vector<int>> cluster(D);
    std::vector<std::vector<int>> best_cluster;
    std::vector<std::vector<int>> pre_cluster;
    std::vector<std::pair<int, int>> history;
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
        if (!history.empty())
        {
          std::pair<int, int> pre = history.back();
          if (pre.first == from && pre.second == to)
          {
            best_cluster = cluster;
          }
          else if (!(pre.first == to && pre.second == from))
          {
            best_cluster = cluster;
          }
          else
          {
            cluster = pre_cluster;
            history.pop_back();
            continue;
          }
        }
        int idx = (u32)mt() % fsize;
        int target_item_idx = cluster[from][idx];
        pre_cluster = cluster;
        cluster[from].erase(cluster[from].begin() + idx);
        cluster[to].push_back(target_item_idx);
        history.emplace_back(from, to);
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

  std::vector<int> generate_virtual_weights(int num_item, int num_div, int seed_)
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

  std::vector<int> merge_sort(const std::vector<int> &items, int &num_query, bool &suspend, bool simulation)
  {
    std::vector<int> ret = items;
    int num_item = items.size();
    std::vector<int> buf(ret.size(), 0);
    int mid = num_item / 2;
    sub_merge_sort(0, mid, ret, buf, num_query, suspend, simulation);
    sub_merge_sort(mid, num_item, ret, buf, num_query, suspend, simulation);
    merge(0, mid, num_item, ret, buf, num_query, suspend, simulation);
    return ret;
  }

  void sub_merge_sort(int l, int r, std::vector<int> &arr, std::vector<int> &buf, int &num_query, bool &suspend, bool simulation)
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
    sub_merge_sort(l, mid, arr, buf, num_query, suspend, simulation);
    sub_merge_sort(mid, r, arr, buf, num_query, suspend, simulation);
    merge(l, mid, r, arr, buf, num_query, suspend, simulation);
  }

  void merge(int l, int m, int r, std::vector<int> &arr, std::vector<int> &buf, int &num_query, bool &suspend, bool simulation)
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
        buf[idx] = arr[i2];
        i2++;
        continue;
      }
      else if (i2 == r)
      {
        buf[idx] = arr[i1];
        i1++;
        continue;
      }
      else
      {
        Comp res = query(1, 1, {arr[i1]}, {arr[i2]}, num_query, server, simulation);
        if (res == Comp::LESS)
        {
          buf[idx] = arr[i1];
          i1++;
        }
        else
        {
          buf[idx] = arr[i2];
          i2++;
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
};

int main(int argv, char *argc[])
{
  start_time = clock();
  Solver solver;

  solver.init();
  // solver.solve();
  solver.merge_sort_test();
  solver.summary();
  return 0;
}