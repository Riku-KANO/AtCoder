// includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <queue>
#include <deque>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <stack>
#include <bitset>
#include <complex>

#include <random>
#include <algorithm>
#include <limits>
#include <utility>
#include <iomanip>
#include <ios>
#include <iterator>
#include <numeric>
#include <chrono>
#include <memory>
#include <optional>
#include <utility>
#include <functional>
#include <thread>

#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <ctime>

// --------------------- macros ----------------------
#define pii std::pair<int, int>
#define Vec std::vector
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define FOR(i, s, n) for (int i = (int)s; i < (int)(n); i++)

#define LOG_INFO(message, ...) fprintf(stderr, "[INFO] " message "\n", ##__VA_ARGS__)
#define LOG_ERROR(message, ...) fprintf(stderr, "[ERROR] " message "\n", ##__VA_ARGS__)
#define LOG_WARNING(message, ...) fprintf(stderr, "[WARNING] " message "\n", ##__VA_ARGS__)

using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int16_t;
using i64 = int16_t;

typedef std::string Commands;
typedef char Command;

// # 考察メモ
//
// 1 → 1 → 4 → 4
// ↑       ↑   ↓      左の経路の1ターンあたりのゴミの獲得量: (1 * 4 * 16 + 4 * 12 * 16) / 16 = 52
// 1   1 ← 4   4      1周のスコア : 300   RRRDDDLLLURRULLU                                     avg a:  8, 30
// ↑   ↓   ↑   ↓      4の地点を2周: 289   RRRDDDLUUURDDDLUULDDLUUU                             avg a: 12, 25
// 1   1   4   4      4の地点を3周: 300   RRRDDDLUUURDDDLUUURDDDLUULDDLUUU                     avg a: 16, 20
// ↑   ↓   ↑   ↓      4の地点を4周: 319   RRRDDDLUUURDDDLUUURDDDLUUURDDDLUULDDLUUU             avg a: 20, 20
// 1 ← 1   4 ← 4      avg aの値が均一になっているときが最適な状態とは限らない。
//
// 1 → 1 → 4 → 4
// ↑       ↑   ↓      左の経路の1ターンあたりのゴミの獲得量: (1 * 4 * 16 + 4 * 12 * 16) / 16 = 52
// 1 ← 1 ← 4   4      1周のスコア : 390   RRRDDDLLLURRULLU                                     avg a:  8, 30
//         ↑   ↓      4の地点を2周: 373   RRRDDDLLLURRUURDDDLLLURRULLU                         avg a: 14, 27
// 4 → 4 → 4   4      4の地点を3周: 380   RRRDDDLLLURRUURDDDLLLURRUURDDDLLLURRULLU             avg a: 20, 25
// ↑           ↓      4の地点を4周: 396   RRRDDDLLLURRUURDDDLLLURRUURDDDLLLURRUURDDDLLLURRULLU avg a: 26, 24
// 4 ← 4 ← 4 ← 4      avg aの値が均一になっているときが最適な状態とは限らない。

// 1 → 1 →100→100     1周のスコア    : 6060   RRRDDDLLLURRULLU                                     avg a:  8, 750
// ↑       ↑   ↓      100の地点を 2周: 5025   RRRDDDLUUURDDDLUULDDLUUU                             avg a: 12, 617
// 1   1 ←100 100     100の地点を 3周: 4524   RRRDDDLUUURDDDLUUURDDDLUULDDLUUU                     avg a: 16, 550
// ↑   ↓   ↑   ↓      100の地点を 4周: 4236   RRRDDDLUUURDDDLUUURDDDLUUURDDDLUULDDLUUU             avg a: 20, 510
// 1   1  100 100     100の地点を 5周: 4055   RRRDDDLUUURDDDLUUURDDDLUUURDDDLUUURDDDLUULDDLUUU     avg a: 24, 483
// ↑   ↓   ↑   ↓      100の地点を 6周: 3934   RRRDDDLUUURDDDLUUURDDDLUUURDDDLUUURDDDLUUURDDDLU...         28, 464
// 1 ← 1  100←100     100の地点を 7周: 3852   RRRDDDLUUURDDDLUUURDDDLUUURDDDLUUURDDDLUUURDDDLU...         32, 450
//                    100の地点を 8周: 3795   略                                                          36, 439
//                    100の地点を 9周: 3756   略                                                          40, 430
//                    100の地点を10周: 3730   略                                                          44, 423
//                    100の地点を11周: 3713   略                                                          48, 417
//                    100の地点を12周: 3704   略                                                          52, 412
//                    100の地点を13周: 3701   略                                                          56, 407
//                    100の地点を14周: 3703   略                                                          60, 403
//                    100の地点を15周: 3708   略                                                          64, 400
//
// 1 1 ......... 1 1    30 * 30グリッド. 右下2*2だけ100
// 1 1 ......... 1 1    1周... (1 * 896 + 100 * 4) * 900 - 1 * 896 - 100 * 4 = 1165104 / 2 = 582552
// | | \                2周... (0 + 900) * 904 +
// | |  \ 
// 1 1 --------100 100
// 1 1 ------- 100 100

// 以下のスコアの通り、マップが広くて汚くなりやすい箇所が限られたスペースの場合はさっさとスタート地点に帰った方がスコアが高くなる。
// 周: score
//  1: 582552.0
//  2: 585146.6548672566
//  3: 587741.2863436124
//  4: 590335.8947368421
//  5: 592930.480349345
//  6: 595525.0434782609
//  7: 598119.5844155845
//  8: 600714.1034482758
//  9: 603308.6008583691
// 10: 605903.0769230769
// 11: 608497.5319148937
// 12: 611091.966101695
// 13: 613686.3797468354
// 14: 616280.7731092437
// 15: 618875.1464435146
// 16: 621469.5
// 17: 624063.8340248963
// 18: 626658.1487603306
// 19: 629252.4444444445

// 40*40マップの10分の１が'128'でそれ以外は'1'である場合、'128'を25回程度ループすると最適値っぽくなる。
//             100分の1が'128'の場合は14回ループ
// N:  40, lv:    1, uv:  128, loop:   25, numL: 1440, score:       9881181.18

// --------------------- constants ----------------------
constexpr int INF = 1 << 30;
constexpr long long LINF = 1LL << 60;
constexpr double EPS = 1e-6;
const double PI = std::acos(-1);

double TL = 1.90;

constexpr int DI[4] = {0, 1, 0, -1};
constexpr int DJ[4] = {1, 0, -1, 0};

constexpr int MAX_N = 40;
// --------------------- global variables ----------------------
std::random_device seedGen;
std::mt19937 mt(seedGen());
std::uniform_int_distribution<> rand01(0, 1);
std::uniform_real_distribution<> randReal(0, 1);
clock_t startTime;

// --------------------- functionss ----------------------
inline bool is_out(int r, int c, int N)
{
  return (r < 0 || r >= N || c < 0 || c >= N);
}

inline double getTime(clock_t startTime)
{
  return (double)(clock() - startTime) / CLOCKS_PER_SEC;
}

bool argParse(int argv, char *argc[])
{
  for (int i = 0; i < argv; i++)
  {
    if (std::strcmp(argc[i], "--seed") == 0)
    {
      if (i + 1 > argv)
      {
        LOG_ERROR("No arguments.");
        return false;
      }
      int _seed = std::stoi(argc[i + 1]);
      mt = std::mt19937(_seed);
    }
    if (std::strcmp(argc[i], "--TL") == 0)
    {
      if (i + 1 > argv)
      {
        LOG_ERROR("No arguments.");
      }
      double _TL = std::stod(argc[i + 1]);
      TL = _TL;
    }
  }
  return true;
}

// --------------------- classes ----------------------
/**
 * @brief Input for Application
 *
 */
struct Input
{
  int N;
  std::array<std::string, MAX_N> horizontal;
  std::array<std::string, MAX_N> vertical;
  std::array<int, MAX_N * MAX_N> dirty;

  // factory method
  static Input getInput();
};

Input Input::getInput()
{
  Input ret;
  std::cin >> ret.N;
  for (int i = 0; i < ret.N - 1; i++)
  {
    std::cin >> ret.horizontal[i];
  }
  for (int i = 0; i < ret.N; i++)
  {
    std::cin >> ret.vertical[i];
  }
  for (int i = 0; i < ret.N; i++)
  {
    for (int j = 0; j < ret.N; j++)
    {
      int idx = i * ret.N + j;
      std::cin >> ret.dirty[idx];
    }
  }
  return ret;
}

// ----------------------------------------------------
struct Field
{
  const int N;
  const int numCell;
  const std::array<std::string, MAX_N> horizontals;
  const std::array<std::string, MAX_N> verticals;
  const std::array<int, MAX_N * MAX_N> dirty;
  Vec<Vec<int>> dist; // [from: (r,c)][to: (r, c)]
  Vec<Vec<int>> graph;
  Vec<Vec<int>> graphDir;       // [from: (r, c)][dir]
  Vec<Vec<Vec<int>>> cellsDist; // [from][dist] -> cells
  Vec<double> meanDistSum;      // [pos (r, c)] -> mean of distance summation
  Vec<int> cellIDs;

  Field(const Input &input); // init from input
  inline int getRCIdx(int r, int c) const;
  inline int getDirVal(int from, int to) const;

private:
  inline void initGraph();
  inline void initDistances();
  inline void initStats();
  inline int _getRCIdx(int r, int c) const;
};

Field::Field(const Input &input)
    : N(input.N),
      numCell(input.N * input.N),
      horizontals(input.horizontal),
      verticals(input.vertical),
      dirty(input.dirty)
{
  // fill cellIDs
  cellIDs.resize(N * N);
  std::iota(cellIDs.begin(), cellIDs.end(), 0);

  // Initialize graph data
  this->initGraph();

  // Initialize distances data
  this->initDistances();

  // Initialize statistic data
  this->initStats();

  LOG_INFO("Field Initialization done at %6.4f s", getTime(startTime));
}

inline int Field::getRCIdx(int r, int c) const
{
  return this->_getRCIdx(r, c);
}

inline int Field::getDirVal(int from, int to) const
{
  int diff = to - from;
  if (diff == 1)
  {
    return 0;
  }
  else if (diff == this->N)
  {
    return 1;
  }
  else if (diff == -1)
  {
    return 2;
  }
  else if (diff == -this->N)
  {
    return 3;
  }
  else
  {
    LOG_ERROR("Invalid positions: %4d, %4d at Line %d", from, to, __LINE__);
    std::exit(-1);
  }
}

inline void Field::initGraph()
{
  const char EMPTY = '0';
  const char WALL = '1';
  this->graph.resize(numCell);
  this->graphDir.resize(numCell, Vec<int>(4, -1));
  int from, to;
  int r, c;
  int nr, nc;
  // horizontal bars
  for (r = 0; r < N - 1; r++)
  {
    for (c = 0; c < N; c++)
    {
      from = _getRCIdx(r, c);
      nr = r + 1;
      nc = c;
      to = _getRCIdx(nr, nc);
      if (horizontals[r][c] == EMPTY)
      {
        graph[from].emplace_back(to);
        graph[to].emplace_back(from);
        graphDir[from][1] = to;
        graphDir[to][3] = from;
      }
    }
  }
  // vertical bars
  for (r = 0; r < N; r++)
  {
    for (c = 0; c < N - 1; c++)
    {
      from = _getRCIdx(r, c);
      nr = r;
      nc = c + 1;
      to = _getRCIdx(nr, nc);
      if (verticals[r][c] == EMPTY)
      {
        graph[from].emplace_back(to);
        graph[to].emplace_back(from);
        graphDir[from][0] = to;
        graphDir[to][2] = from;
      }
    }
  }
} // initGraph()

inline void Field::initDistances()
{
  this->dist.resize(numCell);
  int maxDist = 0;

  // initialize distances
  for (int from : cellIDs)
  {
    this->dist[from].resize(numCell, INF);
    std::queue<int> q;
    q.push(from);
    dist[from][from] = 0;

    while (!q.empty())
    {
      int u = q.front();
      q.pop();
      for (int v : graph[u])
      {
        if (dist[from][v] > dist[from][u] + 1)
        {
          dist[from][v] = dist[from][u] + 1;
          maxDist = std::max(maxDist, dist[from][v]);
          q.push(v);
        }
      }
    }
  }

  LOG_INFO("Max distance: %3d", maxDist);
  this->cellsDist.resize(numCell, Vec<Vec<int>>(maxDist + 1));
  for (int from : cellIDs)
  {
    for (int to : cellIDs)
    {
      int d = dist[from][to];
      cellsDist[from][d].emplace_back(to);
    }
  }
} // Field::initDistances()

inline void Field::initStats()
{
  this->meanDistSum.resize(numCell);

  for (int from : cellIDs)
  {
    long long totalDist = 0LL;
    for (int to : cellIDs)
    {
      totalDist += dist[from][to];
    }

    meanDistSum[from] = (double)totalDist / numCell;
  }
}

inline int Field::_getRCIdx(int r, int c) const
{
  return r * this->N + c;
}

/**
 * @brief Mode enum for bot.
 *
 */
enum class Mode
{
  NOT_VISIT,
  HIGH_VALUE,
};

/**
 * @brief Abstract class of cleaner bots. This is like an interface, not contain field variables.
 *
 */
class CleanerBot
{
public:
  virtual void init(std::shared_ptr<Field> _field) = 0;
  virtual void execute(const Commands &commands) = 0;
};

/**
 * @brief Cleaner bot class
 *
 * @todo `TakahashikunCleanerNo2`か`OsoujiTakahashikunMk2`のどちらの名前がいいだろうか。
 */
class TakahashikunCleanerNo2 : CleanerBot
{
public:
  std::shared_ptr<Field> field;
  int curPos;
  int curTurn;
  Vec<Vec<int>> timeStamp; // [position (r, c)] -> (Turns when you visited the position. 0 value is default.)
  Vec<int> history;        // history of location
  Vec<bool> visited;
  TakahashikunCleanerNo2() {}
  void init(std::shared_ptr<Field> _field) override;
  bool allVisited() const;
  void setCurPos(int pos);
  Commands findPath(bool excludeVisited) const;
  Commands findPath(int to, Mode mode) const;
  Commands findPathCustom() const;
  Commands findTSPPath();
  Commands findTreePath();
  Commands findPathBySA();
  Commands findLargeCyclePath();
  void execute(const Commands &commands) override;
  void revert(int turn);
  Commands generateCommandsFromHistory() const;
  Commands shiftCommands(const Commands &commands) const;
  inline int getNextPosByCommand(int cur, Command command) const;

private:
  int visitCount;
  std::mt19937 rng;
  Commands findPathByNotVisit(int to) const;
  Commands findPathByHighValue(int to) const;
  void generateRandomDoubleBondSpanningTree(Vec<Vec<u16>> *graph);
  void compressDoubleBondTree(Vec<Vec<u16>> *graph);
  inline void processDoubleBond(int from, int to, Vec<Vec<u16>> *graph);
  inline Commands generateCommandsFromDoubleBondSpanningTree(Vec<Vec<u16>> *graph);
  inline Vec<u16> generateInitialSASolution();
  inline Vec<u16> reconstructHistory(const Vec<u16> &cur);
  inline Command getAction(int cur, int nx) const;
  inline int _getNextPosByCommand(int cur, Command command) const;
  inline int getRCIdx(int r, int c) const;
}; // TakahashikunCleanerNo2

void TakahashikunCleanerNo2::init(std::shared_ptr<Field> _field)
{
  this->curPos = 0;
  this->curTurn = 0;
  this->visitCount = 1;
  this->rng = std::mt19937(seedGen());
  this->field = _field;
  this->timeStamp.resize(field->numCell, {0});
  this->history = {0};
  this->visited.resize(field->numCell, false);
  this->visited[0] = true;

  LOG_INFO("Bot Initialization done at %6.4f s", getTime(startTime));
}

bool TakahashikunCleanerNo2::allVisited() const
{
  return this->visitCount == field->N * field->N;
}

void TakahashikunCleanerNo2::setCurPos(int pos)
{
  this->curPos = pos;
}

/**
 * @brief 今いるポジションから平均的に価値の高い最短パスを提案する。
 *
 * @param excludeVisited 既に訪れているポジションを除外するかどうか。default: false
 * @return Commands
 */
Commands TakahashikunCleanerNo2::findPath(bool excludeVisited = false) const
{
  Commands ret;
  int from = curPos;
  Vec<long long> dp(field->numCell, LINF); // TODO: check whether 32 bit or 64 bit
  Vec<int> pre(field->numCell, -1);
  dp[from] = 0;
  pre[from] = from;
  std::queue<std::pair<long long, int>> pq;
  pq.push({0LL, from});

  // BFS search
  while (!pq.empty())
  {
    auto [cost, u] = pq.front();
    pq.pop();

    if (cost > dp[u])
    {
      continue;
    }

    for (int v : field->graph[u])
    {
      int timeStep = field->dist[from][v];
      long long value = -field->dirty[v] * pow((timeStep + curTurn - timeStamp[v].back()), 2);
      if (field->dist[from][v] == field->dist[from][u] + 1 && dp[v] > dp[u] + value)
      {
        dp[v] = dp[u] + value;
        pre[v] = u;
        pq.push({dp[v], v});
      }
    }
  }

  // find best average value
  double highest = -1e9;
  int to = -1;
  for (int pos : field->cellIDs)
  {
    if (pos == from)
    {
      continue;
    }
    if (excludeVisited && visited[pos])
    {
      continue;
    }

    double avgValue = (double)(-dp[pos]) / field->dist[from][pos] * pow(0.99, (double)field->dist[from][pos]);
    if (avgValue > highest)
    {
      highest = avgValue;
      to = pos;
    }
  }

  // Reconstruction
  int cur = to;
  while (cur != from)
  {
    ret.push_back(getAction(pre[cur], cur));
    cur = pre[cur];
  }

  std::reverse(ret.begin(), ret.end());

  return ret;
}

Commands TakahashikunCleanerNo2::findPath(int to, Mode mode) const
{
  switch (mode)
  {
  case Mode::HIGH_VALUE:
    return this->findPathByHighValue(to);
  case Mode::NOT_VISIT:
    return this->findPathByNotVisit(to);
  default:
    LOG_ERROR("Invalid Mode Argument at Line %d", __LINE__);
    std::exit(-1);
  }
}

/**
 * @brief
 *
 * @return Commands
 */
Commands TakahashikunCleanerNo2::findPathCustom() const
{
  typedef std::pair<double, int> pdi;
  Commands ret;
  int from = curPos;
  Vec<double> dp(field->numCell, 1e9); // TODO: check whether 32 bit or 64 bit
  Vec<int> pre(field->numCell, -1);
  dp[from] = 0;
  pre[from] = from;
  std::priority_queue<pdi, Vec<pdi>, std::greater<pdi>> pq;
  pq.push({0.0, from});

  // dijkstra search
  while (!pq.empty())
  {
    auto [cost, u] = pq.top();
    pq.pop();

    if (cost > dp[u])
    {
      continue;
    }

    for (int v : field->graph[u])
    {
      int timeStep = field->dist[from][v];
      double value = -field->meanDistSum[v] * field->dirty[v] * (timeStep + curTurn - timeStamp[v].back());
      if (field->dist[from][v] == field->dist[from][u] + 1 && dp[v] > dp[u] + value)
      {
        dp[v] = dp[u] + value;
        pre[v] = u;
        pq.push({dp[v], v});
      }
    }
  }

  // find best average value
  double highest = -1e9;
  int to = -1;
  for (int pos : field->cellIDs)
  {
    if (pos == from)
    {
      continue;
    }

    double avgValue = -dp[pos] / field->dist[from][pos];
    if (avgValue > highest)
    {
      highest = avgValue;
      to = pos;
    }
  }

  // Reconstruction
  int cur = to;
  while (cur != from)
  {
    ret.push_back(getAction(pre[cur], cur));
    cur = pre[cur];
  }

  std::reverse(ret.begin(), ret.end());

  return ret;
}

/**
 * @brief
 *
 * @return Commands
 */
Commands TakahashikunCleanerNo2::findTSPPath()
{
  Commands ret;
  const int ORIGIN = 0;
  Vec<Vec<u16>> spanningTree;
  this->generateRandomDoubleBondSpanningTree(&spanningTree);

  this->compressDoubleBondTree(&spanningTree);
  int cur = ORIGIN;

#ifdef DEBUG
  LOG_INFO("Started tsp search");
#endif
  while (true)
  {
    int singleBondDir = -1;
    int doubleBondDir = -1;
    int nx = -1;
    for (int nei : field->graph[cur])
    {
      int dir = field->getDirVal(cur, nei);
      if (spanningTree[cur][dir] == 2)
      {
        doubleBondDir = dir;
        nx = nei;
        break;
      }
      else if (spanningTree[cur][dir] == 1)
      {
        singleBondDir = dir;
        nx = nei;
      }
    }

    // Not found
    if (nx == -1)
    {
      break;
    }

    // when double bond exist
    int dir;
    int invDir;
    if (doubleBondDir != -1)
    {
      dir = doubleBondDir;
    }
    else
    {
      dir = singleBondDir;
    }
    invDir = (dir + 2) % 4;
    spanningTree[cur][dir]--;
    spanningTree[nx][invDir]--;

    ret += getAction(cur, nx);
    cur = nx;
  }

#ifdef DEBUG
  LOG_INFO("Ended TSP search");
#endif
  return ret;
}

/**
 * @brief
 *
 * @return Commands
 * @todo 次数のグラフから有向グラフの作成
 */
Commands TakahashikunCleanerNo2::findTreePath()
{
  Commands ret;
  const int ORIGIN = 0;
  Vec<Vec<u16>> spanningTree;
  this->generateRandomDoubleBondSpanningTree(&spanningTree);
  return this->generateCommandsFromDoubleBondSpanningTree(&spanningTree);
}

Commands TakahashikunCleanerNo2::findPathBySA()
{
  Commands commands;
  const double penalty = 1e9;
  LOG_INFO("Started finding path by SA");
  auto evalHistory = [&](const Vec<u16> &his) -> double
  {
    double ret = 0.0;
    int size = his.size();
    Vec<Vec<int>> ts(field->numCell);
    for (int cellID : field->cellIDs)
    {
      ts[cellID].resize(1, 0);
    }
    int turn = 0;
    for (int i = 0; i < size; i++)
    {
      int cellID = his[i];
      ts[cellID].emplace_back(turn);
      turn++;
    }
    int curPos = 0;
    int nxPos;
    int length = his.size();
    turn = his.size() - 1;
    long long deltaS = 0;
    long long totalS = 0;
    long long curS = 0;
    for (int cellID : field->cellIDs)
    {
      curS += (long long)field->dirty[cellID] * (turn - ts[cellID].back());
      deltaS += field->dirty[cellID];
      if (ts[cellID].back() == 0)
      {
        ret += penalty;
      }
    }

    for (int i = 1; i < size; i++)
    {
      totalS += curS;
      nxPos = his[i];
      turn++;
      curS += deltaS;
      curS -= (long long)field->dirty[nxPos] * (turn - ts[nxPos].back());
      ts[nxPos].emplace_back(turn);
      curPos = nxPos;
    }
    ret += (double)totalS / length;
    return ret;
  }; // evalHistory

  Vec<u16> curHistory = this->generateInitialSASolution();
  LOG_INFO("generated initial solution");
  Vec<u16> bestHistory = curHistory;
  double curScore = evalHistory(curHistory);
  double bestScore = curScore;
  LOG_INFO("Initial score: %12.2f", bestScore);

  while (getTime(startTime) < 1.0)
  {
    Vec<u16> nxHistory = reconstructHistory(curHistory);
    LOG_INFO("reconstruct done");
    double nxScore = evalHistory(nxHistory);
    LOG_INFO("calculation done");
    if (nxScore < curScore)
    {
      curHistory = nxHistory;
      curScore = nxScore;
      if (nxScore < bestScore)
      {
        LOG_INFO("Best score is updated! %12.2f -> %12.2f", bestScore, nxScore);
        bestHistory = nxHistory;
        bestScore = nxScore;
        Commands bestCommands;
        int size = bestHistory.size();
        for (int i = 0; i < size - 1; i++)
        {
          int from = bestHistory[i];
          int to = bestHistory[i + 1];
          bestCommands += getAction(from, to);
        }
        std::cout << bestCommands << "\n";
      }
    }
  }

  // generate commands from history
  int size = bestHistory.size();
  for (int i = 0; i < size - 1; i++)
  {
    int from = curHistory[i];
    int to = curHistory[i + 1];
    commands += getAction(from, to);
  }

  return commands;
}

Commands TakahashikunCleanerNo2::findLargeCyclePath()
{
  Commands ret;
  const int ORIGIN = 0;
  const int EAST = 0;
  const int SOUTH = 1;
  const int WEST = 2;
  const int NORTH = 3;
  int curDir = 0;
  int curPos = ORIGIN;
  Vec<u16> his = {ORIGIN};
  Vec<bool> curVisited(field->numCell, false);
  auto turnRight = [](int dir) -> int
  {
    dir++;
    if (dir == 4)
      dir = 0;
    return dir;
  };
  auto turnLeft = [](int dir) -> int
  {
    dir--;
    if (dir == -1)
      dir = 3;
    return dir;
  };
  auto turnBack = [](int dir) -> int
  {
    dir += 2;
    if (dir >= 4)
      dir -= 4;
    return dir;
  };
  auto getPos = [&](int pos, int dir) -> int
  {
    switch (dir)
    {
    case 0:
      return pos + 1;
    case 1:
      return pos + field->N;
    case 2:
      return pos - 1;
    case 3:
      return pos - field->N;
    default:
      LOG_ERROR("Invalid position at Line %d", __LINE__);
      std::exit(-1);
    }
  };

  while (true)
  {
    if (curPos == ORIGIN && curDir == NORTH)
    {
      break;
    }
    int leftDir = turnLeft(curDir);
    int leftPos = getPos(curPos, leftDir);
    int forwardPos = getPos(curPos, curDir);
    bool isLeftAdjacent = false;
    bool isForwardAdjacent = false;
    for (int nei : field->graph[curPos])
    {
      if (nei == leftPos)
      {
        isLeftAdjacent = true;
        break;
      }
      else if (nei == forwardPos)
      {
        isForwardAdjacent = true;
      }
    }
    if (isLeftAdjacent)
    {
      curDir = leftDir;
      curPos = leftPos;
      his.emplace_back(curPos);
      curVisited[curPos] = true;
    }
    else if (isForwardAdjacent)
    {
      curPos = forwardPos;
      his.emplace_back(curPos);
      curVisited[curPos] = true;
    }
    else
    {
      curDir = turnRight(curDir);
    }
  }

  while (true)
  {
    bool exist = false;

    int hisSize = his.size();
    for (int i = 0; i < hisSize - 1; i++)
    {
      int from = his[i];
      int to = his[i + 1];
      int dir = field->getDirVal(from, to);

      if (dir % 2 == 0)
      { // N / S
        if (field->graphDir[from][SOUTH] != -1 && field->graphDir[to][SOUTH] != -1)
        {
          u16 posf = getPos(from, SOUTH);
          u16 post = getPos(to, SOUTH);
          if (!curVisited[posf] && !curVisited[post] && field->graphDir[posf][dir] != -1)
          {
            Vec<u16> newPath = {posf, post};
            his.insert(his.begin() + i + 1, newPath.begin(), newPath.end());
            exist = true;
            curVisited[posf] = curVisited[post] = true;
            break;
          }
        }
        if (field->graphDir[from][NORTH] != -1 && field->graphDir[to][NORTH] != -1)
        {
          int posf = getPos(from, NORTH);
          int post = getPos(to, NORTH);
          if (!curVisited[posf] && !curVisited[post] && field->graphDir[posf][dir] != -1)
          {
            Vec<u16> newPath = {(u16)posf, (u16)post};
            his.insert(his.begin() + i + 1, newPath.begin(), newPath.end());
            exist = true;
            curVisited[posf] = curVisited[post] = true;
            break;
          }
        }
      }
      else
      { // E / W
        if (field->graphDir[from][EAST] != -1 && field->graphDir[to][EAST] != -1)
        {
          int posf = getPos(from, EAST);
          int post = getPos(to, EAST);
          if (!curVisited[posf] && !curVisited[post] && field->graphDir[posf][dir] != -1)
          {
            Vec<u16> newPath = {(u16)posf, (u16)post};
            his.insert(his.begin() + i + 1, newPath.begin(), newPath.end());
            exist = true;
            curVisited[posf] = curVisited[post] = true;
            break;
          }
        }
        if (field->graphDir[from][WEST] != -1 && field->graphDir[to][WEST] != -1)
        {
          int posf = getPos(from, WEST);
          int post = getPos(to, WEST);
          if (!curVisited[posf] && !curVisited[post] && field->graphDir[posf][dir] != -1)
          {
            Vec<u16> newPath = {(u16)posf, (u16)post};
            his.insert(his.begin() + i + 1, newPath.begin(), newPath.end());
            exist = true;
            curVisited[posf] = curVisited[post] = true;
            break;
          }
        }
      }
    }

    if (!exist)
    {
      break;
    }
  }

  // handle not visit
  while (true)
  {
    bool exist = false;
    int size = his.size();
    for (int i = 0; i < size - 1 && !exist; i++)
    {
      int u = his[i];
      for (int nei : field->graph[u])
      {
        if (!curVisited[nei])
        {
          Vec<u16> newPath = {(u16)nei, (u16)u};
          his.insert(his.begin() + i + 1, newPath.begin(), newPath.end());
          curVisited[nei] = true;
          exist = true;
          break;
        }
      }
    }
    if (!exist)
    {
      break;
    }
  }

  int hisSize = his.size();
  for (int i = 0; i < hisSize - 1; i++)
  {
    int from = his[i];
    int to = his[i + 1];
    ret.push_back(getAction(from, to));
  }

  return ret;
}

void TakahashikunCleanerNo2::execute(const Commands &commands)
{
  int nxPos;
  for (Command command : commands)
  {
    nxPos = _getNextPosByCommand(curPos, command);
    curPos = nxPos;
    curTurn++;
    timeStamp[curPos].emplace_back(curTurn);
    if (!visited[curPos])
    {
      visitCount++;
      visited[curPos] = true;
    }
    history.emplace_back(curPos);
  }
}

void TakahashikunCleanerNo2::revert(int turn)
{
  while (history.size() != turn + 1)
  {
    int lastPos = history.back();
    history.pop_back();
    timeStamp[lastPos].pop_back();
  }
  curTurn = turn;
  curPos = history.back();
  assert(history.size() >= 1);
}

Commands TakahashikunCleanerNo2::generateCommandsFromHistory() const
{
  Commands ret;
  const int size = history.size() - 1;
  int cur, nx;
  for (int i = 0; i < size; i++)
  {
    cur = history[i];
    nx = history[i + 1];
    ret += getAction(cur, nx);
  }
  return ret;
}

Commands TakahashikunCleanerNo2::shiftCommands(const Commands &commands) const
{
  Vec<int> his = this->history;
  Commands ret = commands;
  his.pop_back();
  while (his.front() != 0)
  {
    std::rotate(his.begin(), his.begin() + 1, his.end());
    std::rotate(ret.begin(), ret.begin() + 1, ret.end());
  }
  return ret;
}

int TakahashikunCleanerNo2::getNextPosByCommand(int cur, Command command) const
{
  return this->_getNextPosByCommand(cur, command);
}

/**
 * @brief 価値の高い最短パスを通るようにパスを提案する。ダイクストラで見つける。
 *
 * @param to 行き先
 * @return Commands
 */
Commands TakahashikunCleanerNo2::findPathByHighValue(int to) const
{
  Commands ret;
  int from = curPos;
  Vec<int> dp(field->numCell, INF); // TODO: check whether 32 bit or 64 bit
  Vec<int> pre(field->numCell, -1);
  dp[from] = 0;
  pre[from] = from;
  std::priority_queue<pii, Vec<pii>, std::greater<pii>> pq;
  pq.push({0, from});

  // dijkstra search
  while (!pq.empty())
  {
    auto [cost, u] = pq.top();
    pq.pop();

    if (cost > dp[u])
    {
      continue;
    }

    for (int v : field->graph[u])
    {
      int timeStep = field->dist[from][v];
      int value = -field->dirty[v] * (timeStep + curTurn - timeStamp[v].back());
      if (field->dist[from][v] == field->dist[from][u] + 1 && dp[v] > dp[u] + value)
      {
        dp[v] = dp[u] + value;
        pre[v] = u;
        if (v != to)
        {
          pq.push({dp[v], v});
        }
      }
    }
  }

  // Reconstruction
  int cur = to;
  while (cur != from)
  {
    ret.push_back(getAction(pre[cur], cur));
    cur = pre[cur];
  }

  std::reverse(ret.begin(), ret.end());

  return ret;
}

/**
 * @brief 「訪れていない箇所が多い最短パス」かつ「価値の高い最短パス」をみたすようにパスを提案する。
 *
 * @param to 行き先
 * @return Commands
 *
 * @todo 実装する?
 */
Commands TakahashikunCleanerNo2::findPathByNotVisit(int to) const
{
  Commands ret;
  // todo
  return ret;
}

void TakahashikunCleanerNo2::generateRandomDoubleBondSpanningTree(Vec<Vec<u16>> *graph)
{
#ifdef DEBUG
  LOG_INFO("Started generating random double bond spanning tree");
#endif
  graph->resize(field->numCell, Vec<u16>(4, 0));
  // int startPos = (u16)this->rng() % field->numCell;
  int startPos = 0;
  std::queue<int> q;
  Vec<bool> visited(field->numCell, false);
  q.push(startPos);
  visited[startPos] = true;

  int u;
  while (!q.empty())
  {
    u = q.front();
    q.pop();
    for (int v : field->graph[u]) // TODO: randomize
    {
      if (visited[v])
      {
        continue;
      }
      int dir = field->getDirVal(u, v);
      int invDir = (dir + 2) % 4;
      q.push(v);
      visited[v] = true;
      graph->at(u)[dir] = 2;    // set double bond
      graph->at(v)[invDir] = 2; // set double bond
    }
  }
#ifdef DEBUG
  LOG_INFO("Ended generating random double bond spanning tree");
#endif
}

void TakahashikunCleanerNo2::compressDoubleBondTree(Vec<Vec<u16>> *graph)
{
  Vec<int> cellIDs(field->numCell);
  std::iota(cellIDs.begin(), cellIDs.end(), 0);
#ifdef DEBUG
  int searchCount = 0;
#endif
  int iter = 1000;
  while (iter--)
  {
#ifdef DEBUG
    searchCount++;
#endif
    int from = (u32)rng() % field->numCell;
    for (int nei : field->graph[from])
    {
      LOG_INFO("from: %4d, nei: %4d", from, nei);
      int dir = field->getDirVal(from, nei);
      if (graph->at(from)[dir] > 0)
      {
        continue;
      }
      this->processDoubleBond(from, nei, graph);
    }
  }

#ifdef DEBUG
  LOG_INFO("search count of compressing double bond: %d", searchCount);
#endif
}

inline void TakahashikunCleanerNo2::processDoubleBond(int from, int to, Vec<Vec<u16>> *graph)
{
  Vec<int> dist(field->numCell, INF);
  Vec<int> pre(field->numCell, -1);
  std::queue<int> q;
  dist[from] = 0;
  pre[from] = from;
  q.push(from);
  while (!q.empty() && dist[to] == INF)
  {
    int u = q.front();
    q.pop();
    for (int v : field->graph[u])
    {
      int dir = field->getDirVal(u, v);
      if (graph->at(u)[dir] == 0)
      {
        continue;
      }
      if (dist[v] > dist[u] + 1)
      {
        dist[v] = dist[u] + 1;
        pre[v] = u;
        q.push(v);
      }
    }
  }
#ifdef DEBUG
  LOG_INFO("dist[from][to] = %d", dist[to]);
#endif
  if (dist[to] == INF)
    return;
  int doubleBond = 0;
  int singleBond = 0;
  int cur = to;

  Vec<int> path = {to};
  while (cur != from)
  {

    LOG_INFO("pre: %4d, cur: %4d", pre[cur], cur);
    int dir = field->getDirVal(pre[cur], cur);
    if (graph->at(pre[cur])[dir] == 1)
    {
      singleBond++;
    }
    else
    {
      doubleBond++;
    }
    cur = pre[cur];
    path.emplace_back(cur);
  }

  if (doubleBond <= singleBond)
  {
    return;
  }

  // increment degree of edge between from and to
  int dir = field->getDirVal(from, to);
  int invDir = (dir + 2) % 4;
  graph->at(from)[dir]++;
  graph->at(to)[invDir]++;

  // decrement degree of edges on path
  int pathSize = path.size();
  for (int i = 0; i < pathSize - 1; i++)
  {
    int u = path[i];
    int v = path[i + 1];
    LOG_INFO("u: %4d, v: %4d", u, v);
    int dir = field->getDirVal(u, v);
    int invDir = (dir + 2) % 4;
    graph->at(u)[dir]--;
    graph->at(v)[invDir]--;
  }

  // if the cell is isorated, connect it to neighbor
  // exclude cells; 'from' and 'to'
  for (int i = 1; i < pathSize; i++)
  {
    int cellID = path[i];
    bool isIsorated = true;
    for (int d = 0; d < 4; d++)
    {
      if (graph->at(cellID)[d] > 0)
      {
        isIsorated = false;
        break;
      }
    }
    if (!isIsorated)
    {
      continue;
    }

    int targetNei = -1;

    for (int nei : field->graph[cellID])
    {
      bool isNeiIsorated = true;
      for (int d = 0; d < 4; d++)
      {
        if (graph->at(nei)[d] > 0)
        {
          isNeiIsorated = false;
          break;
        }
      }

      if (!isNeiIsorated)
      {
        targetNei = nei;
        break;
      }
    }

#ifdef DEBUG
    if (targetNei == -1)
    {
      LOG_ERROR("Not found target neighbor at Line %d", __LINE__);
    }
#endif

    int dir = field->getDirVal(cellID, targetNei);
    int invDir = (dir + 2) % 4;
    graph->at(cellID)[dir] += 2;
    graph->at(targetNei)[invDir] += 2;
  }
}

inline Commands TakahashikunCleanerNo2::generateCommandsFromDoubleBondSpanningTree(Vec<Vec<u16>> *graph)
{
  Commands ret;
  Vec<Vec<u16>> edgePair(field->numCell, Vec<u16>(4)); // + 字のセル、かつそれぞれの辺の次数が1の場合をチェックする。
  for (int cellID : field->cellIDs)
  {
    int numSingleBond = 0;
    for (int nei : field->graph[cellID])
    {
      int dir = field->getDirVal(cellID, nei);
      if (graph->at(cellID)[dir] == 1)
      {
        numSingleBond++;
      }
    }
    if (numSingleBond == 2)
    {
      for (int nei : field->graph[cellID])
      {
        int dir = field->getDirVal(cellID, nei);
        if (graph->at(cellID)[dir] == 1)
        {
          edgePair[cellID][dir] = 0;
        }
      }
    }

    else if (numSingleBond == 4)
    {
      Vec<int> labels(field->numCell, -1);
      std::queue<std::pair<int, int>> q; //(cellID, label)
      bool done = false;
      for (int nei : field->graph[cellID])
      {
        int dirLabel = field->getDirVal(cellID, nei);
        q.push({nei, dirLabel});
        labels[nei] = dirLabel;
      }

      while (!q.empty() && !done)
      {
        auto [u, label] = q.front();
        q.pop();
        for (int v : field->graph[u])
        {
          int dir = field->getDirVal(u, v);
          if (graph->at(u)[dir] == 0 || graph->at(u)[dir] == 2)
          {
            continue;
          }
          if (v == cellID || labels[v] == label)
          {
            continue;
          }

          // if adjacent cell label is already filled by other label, search ends
          if (labels[v] != -1)
          {
            edgePair[cellID][label] = 0;
            edgePair[cellID][labels[v]] = 0;
            for (int d = 0; d < 4; d++)
            {
              if (edgePair[cellID][d] == -1)
              {
                edgePair[cellID][d] = 1;
              }
            }
            done = true;
            break;
          }
          labels[v] = label;
          q.push({v, label});
        }
      }
    }
  }

  auto dfs = [&](auto dfs, int start, int pre) -> void
  {
    int numDoubleBond = 0;
    int numSingleBond = 0;
    int nx;
    for (int nei : field->graph[start])
    {
      int dir = field->getDirVal(start, nei);
      if (graph->at(start)[dir] == 1)
      {
        numSingleBond++;
      }
      else if (graph->at(start)[dir] == 2)
      {
        numDoubleBond++;
      }
    }

    if (numSingleBond == 3)
    { // ペアが二つある方向へ行く。
      Vec<int> label0;
      Vec<int> label1;
      for (int nei : field->graph[start])
      {
        int dir = field->getDirVal(start, nei);
        if (graph->at(start)[dir] == 1)
        {
          if (edgePair[start][dir] == 0)
          {
            label0.emplace_back(nei);
          }
          else
          {
            label1.emplace_back(nei);
          }
        }
      }

      if (label0.size() == 2)
      {
        nx = label0.front();
      }
      else
      {
        nx = label1.front();
      }
    }
    else if (numDoubleBond > 0)
    { // ダブルボンドがあればそこを優先。
      for (int nei : field->graph[start])
      {
        int dir = field->getDirVal(start, nei);
        if (graph->at(start)[dir] == 2)
        {
          nx = nei;
          break;
        }
      }
    }
    else if (numSingleBond == 1)
    {
      for (int nei : field->graph[start])
      {
        int dir = field->getDirVal(start, nei);
        if (graph->at(start)[dir] == 1)
        {
          nx = nei;
          break;
        }
      }
    }
    else
    {
      return;
    }

    int dir = field->getDirVal(start, nx);
    int invDir = (dir + 2) % 4;
    graph->at(start)[dir]--;
    graph->at(nx)[invDir]--;
    ret += getAction(start, nx);
    dfs(dfs, nx, start);
  };

  dfs(dfs, 0, -1);
  return ret;
}

inline Vec<u16> TakahashikunCleanerNo2::generateInitialSASolution()
{
  Vec<u16> ret;

  Commands commands = this->findLargeCyclePath();
  u16 cur = 0;
  ret.emplace_back(cur);
  for (Command command : commands)
  {
    switch (command)
    {
    case 'R':
      cur += 1;
      break;
    case 'L':
      cur -= 1;
      break;
    case 'U':
      cur -= field->N;
      break;
    case 'D':
      cur += field->N;
      break;
    default:
      break;
    }
    ret.emplace_back(cur);
  }
  return ret;
}

inline Vec<u16> TakahashikunCleanerNo2::reconstructHistory(const Vec<u16> &cur)
{
  Vec<u16> ret = cur;

  int curSize = cur.size();
  int length = (u32)rng() % 5 + 2;
  int start = (u32)rng() % (curSize - length);
  u16 startPos = ret[start];
  u16 endPos = ret[start + length - 1];
  ret.erase(ret.begin() + start, ret.begin() + start + length);
  Vec<u16> path;
  Vec<int> dist(field->numCell, INF);
  Vec<int> pre(field->numCell, -1);
  dist[start] = 0;
  pre[start] = start;
  auto dfs = [&](auto dfs, int s) -> void
  {
    if (dist[endPos] != INF)
    {
      return;
    }
    Vec<int> indices(field->graph[s].size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    for (int idx : indices)
    {
      int nei = field->graph[s][idx];
      if (dist[nei] != INF)
      {
        continue;
      }
      dist[nei] = dist[s] + 1;
      pre[nei] = s;
      dfs(dfs, nei);
    }
  };

  dfs(dfs, startPos);
  LOG_INFO("DFS DONE");
  LOG_INFO("start pos: %d, end pos: %d", startPos, endPos);
  u16 pos = endPos;
  path = {endPos};
  while (pos != startPos)
  {
    pos = pre[pos];
    path.emplace_back(pos);
  }
  for (int p : path)
  {
    std::cerr << p << " ";
  }
  std::cerr << "\n";
  std::reverse(path.begin(), path.end());
  ret.insert(ret.begin() + start, path.begin(), path.end());

  return ret;
}

inline Command TakahashikunCleanerNo2::getAction(int cur, int nx) const
{
  if (nx - cur == field->N)
  {
    return 'D';
  }
  else if (nx - cur == 1)
  {
    return 'R';
  }
  else if (nx - cur == -1)
  {
    return 'L';
  }
  else if (nx - cur == -field->N)
  {
    return 'U';
  }
  else
  {
    LOG_ERROR("Invalid Positions!, %d, %d at Line %d", nx, cur, __LINE__);
    std::exit(-1);
  }
}

inline int TakahashikunCleanerNo2::_getNextPosByCommand(int cur, Command command) const
{
  switch (command)
  {
  case 'L':
    return cur - 1;
  case 'R':
    return cur + 1;
  case 'U':
    return cur - field->N;
  case 'D':
    return cur + field->N;
  default:
    LOG_ERROR("Invalid command '%c' at Line %d", command, __LINE__);
    std::exit(-1);
  }
}

inline int TakahashikunCleanerNo2::getRCIdx(int r, int c) const
{
  return this->field->getRCIdx(r, c);
}

/**
 * @brief Main App class
 *
 */
class App
{
public:
  App(const Input &input);
  void init();
  void run();
  void test();
  void show() const;
  void summary() const;

private:
  const Input initialInput;
  TakahashikunCleanerNo2 takahashi;
  std::shared_ptr<Field> field; // bot の cleanerと共有する
  Commands bestCommands;
  long long bestScore;
  std::mt19937 rng;
  long long calcScore() const;
  long long calcScore(const Commands &commands) const;
};

/**
 * @brief Construct a new App:: App object.
 *
 * @param input
 */
App::App(const Input &input) : initialInput(input)
{
  this->field = std::make_shared<Field>(input);
}

/**
 * @brief 初期化メソッド
 *
 */
void App::init()
{
  this->takahashi.init(this->field);
  this->rng = std::mt19937(seedGen());
  LOG_INFO("Application Initialization done at %6.4f s", getTime(startTime));
  LOG_INFO("Takahashi-kun is Ready!");
}

/**
 * @brief 問題を解くメインメソッド
 *
 */
void App::run()
{
  Commands commands1 = this->takahashi.findLargeCyclePath();
  Commands commands2;
  Commands commands;
  long long score1 = this->calcScore(commands1);
  long long score2;
  bestScore = score1;
  bestCommands = commands1;
  bool excludeVisited = true;
  int startPos = 0;
  int dirtiest = -1;
  for (int cellID : field->cellIDs)
  {
    if (field->dirty[cellID] > dirtiest)
    {
      dirtiest = field->dirty[cellID];
      int startPos = cellID;
    }
  }
  this->takahashi.setCurPos(startPos);

  while (!takahashi.allVisited())
  {
    if ((u32)rng() % 10 == 0)
    {
      commands = this->takahashi.findPath(excludeVisited);
    }
    else
    {
      commands = this->takahashi.findPath();
    }
    commands2 += commands;
    this->takahashi.execute(commands);
  }
  commands = this->takahashi.findPath(startPos, Mode::HIGH_VALUE);
  commands2 += commands;
  score2 = this->calcScore(commands2);
  LOG_INFO("score1: %12lld, score2: %12lld", score1, score2);

  if (score2 < bestScore)
  {
    this->bestCommands = takahashi.shiftCommands(commands2);
    bestScore = score2;
  }
}

void App::test()
{
}

void App::show() const
{
  std::cout << bestCommands << "\n";
  std::cout.flush();
}

void App::summary() const
{
  fprintf(stderr, "\n########## SUMMARY #########\n");
  fprintf(stderr, "Score       : %10lld pt\n", bestScore);
  fprintf(stderr, "Elapsed Time: %10.2f s\n", getTime(startTime));
  fprintf(stderr, "############################\n");
}

long long App::calcScore() const
{
  Vec<Vec<int>> timeStamp = this->takahashi.timeStamp;
  int curTurn = this->takahashi.curTurn;
  int curPos = 0;
  int nxPos;
  int length = bestCommands.size();
  long long deltaS = 0;
  long long totalS = 0;
  long long curS = 0;
  for (int cellID : field->cellIDs)
  {
    curS += (long long)field->dirty[cellID] * (curTurn - timeStamp[cellID].back());
    deltaS += field->dirty[cellID];
  }

  for (Command command : bestCommands)
  {
    totalS += curS;
    nxPos = this->takahashi.getNextPosByCommand(curPos, command);
    curTurn++;
    curS += deltaS;
    curS -= (long long)field->dirty[nxPos] * (curTurn - timeStamp[nxPos].back());
    timeStamp[nxPos].emplace_back(curTurn);
    curPos = nxPos;
  }

  return std::round((double)totalS / length);
}

long long App::calcScore(const Commands &commands) const
{
  Vec<Vec<int>> timeStamp(field->numCell);
  int curTurn = 0;
  int curPos = 0;
  for (int cellID : field->cellIDs)
  {
    timeStamp[cellID] = {0};
  }
  for (Command command : commands)
  {
    curPos = this->takahashi.getNextPosByCommand(curPos, command);
    curTurn++;
    timeStamp[curPos].emplace_back(curTurn);
  }

  int nxPos;
  int length = commands.size();
  long long deltaS = 0;
  long long totalS = 0;
  long long curS = 0;
  for (int cellID : field->cellIDs)
  {
    curS += (long long)field->dirty[cellID] * (curTurn - timeStamp[cellID].back());
    deltaS += field->dirty[cellID];
  }

  for (Command command : commands)
  {
    totalS += curS;
    nxPos = this->takahashi.getNextPosByCommand(curPos, command);
    curTurn++;
    curS += deltaS;
    curS -= (long long)field->dirty[nxPos] * (curTurn - timeStamp[nxPos].back());
    timeStamp[nxPos].emplace_back(curTurn);
    curPos = nxPos;
  }

  return std::round((double)totalS / length);
}

int main(int argv, char *argc[])
{
  startTime = clock();

  // I/O optimization
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  std::cout.tie(nullptr);
  std::cin.rdbuf()->pubsetbuf(nullptr, 0);
  std::cout.rdbuf()->pubsetbuf(nullptr, 0);

  // Read input
  Input input = Input::getInput();

  // Run Application
  App app(input);
  app.init();
  app.run();
  // app.test();
#ifdef LOCAL
  app.summary();
#endif
  app.show();

  return 0;
}