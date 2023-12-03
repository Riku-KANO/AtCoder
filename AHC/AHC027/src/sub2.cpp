// 近距離のものを重点的に探索

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
  const std::array<std::string, MAX_N> horizontals;
  const std::array<std::string, MAX_N> verticals;
  const std::array<int, MAX_N * MAX_N> dirty;
  Vec<Vec<int>> dist; // [from: (r,c)][to: (r, c)]
  Vec<Vec<int>> graph;
  Vec<Vec<Vec<int>>> cellsDist; // [from][dist] -> cells

  Field(const Input &input); // init from input

private:
  inline int getRCIdx(int r, int c) const;
  inline void initGraph();
  inline void initDistances();
};

Field::Field(const Input &input) : N(input.N), horizontals(input.horizontal), verticals(input.vertical), dirty(input.dirty)
{
  // Initialize graph data
  this->initGraph();

  // Initialize distances data
  this->initDistances();

  LOG_INFO("Field Initialization done at %6.4f s", getTime(startTime));
}

inline void Field::initGraph()
{
  const char EMPTY = '0';
  const char WALL = '1';
  this->graph.resize(N * N);
  // horizontal bars
  for (int r = 0; r < N - 1; r++)
  {
    for (int c = 0; c < N; c++)
    {
      int from = getRCIdx(r, c);
      int nr = r + 1;
      int nc = c;
      int to = getRCIdx(nr, nc);
      if (horizontals[r][c] == EMPTY)
      {
        graph[from].emplace_back(to);
        graph[to].emplace_back(from);
      }
    }
  }
  // vertical bars
  for (int r = 0; r < N; r++)
  {
    for (int c = 0; c < N - 1; c++)
    {
      int from = getRCIdx(r, c);
      int nr = r;
      int nc = c + 1;
      int to = getRCIdx(nr, nc);
      if (verticals[r][c] == EMPTY)
      {
        graph[from].emplace_back(to);
        graph[to].emplace_back(from);
      }
    }
  }
} // initGraph()

inline void Field::initDistances()
{
  const int NUM_CELL = N * N;
  this->dist.resize(NUM_CELL);
  int maxDist = 0;

  // initialize distances
  for (int from = 0; from < NUM_CELL; from++)
  {
    this->dist[from].resize(NUM_CELL, INF);
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
  this->cellsDist.resize(NUM_CELL, Vec<Vec<int>>(maxDist + 1));
  for (int from = 0; from < NUM_CELL; from++)
  {
    for (int to = 0; to < NUM_CELL; to++)
    {
      int d = dist[from][to];
      cellsDist[from][d].emplace_back(to);
    }
  }
} // Field::initDistances()

inline int Field::getRCIdx(int r, int c) const
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
  virtual void execute(const Commands &commands, int limit) = 0;
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
  Commands findPath(bool excludeVisited) const;
  Commands findPath(int to, Mode mode) const;
  void execute(const Commands &commands, int limit) override;
  void revert(int turn);
  Commands generateCommandsFromHistory() const;

private:
  int visitCount;
  Commands findPathByNotVisit(int to) const;
  Commands findPathByHighValue(int to) const;
  inline int getRCIdx(int r, int c) const;
  inline char getAction(int cur, int nx) const;
}; // TakahashikunCleanerNo2

void TakahashikunCleanerNo2::init(std::shared_ptr<Field> _field)
{
  const int NUM_CELL = _field->N * _field->N;
  this->curPos = 0;
  this->curTurn = 0;
  this->visitCount = 1;
  this->field = _field;
  this->timeStamp.resize(NUM_CELL, {0});
  this->history = {0};
  this->visited.resize(NUM_CELL, false);
  this->visited[0] = true;

  LOG_INFO("Bot Initialization done at %6.4f s", getTime(startTime));
}

bool TakahashikunCleanerNo2::allVisited() const
{
  return this->visitCount == field->N * field->N;
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
  const int NUM_CELL = field->N * field->N;
  int from = curPos;
  Vec<int> dp(NUM_CELL, INF); // TODO: check whether 32 bit or 64 bit
  Vec<int> pre(NUM_CELL, -1);
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
        pq.push({dp[v], v});
      }
    }
  }

  // find best average value
  double highest = -1e9;
  int to = -1;
  for (int pos = 0; pos < NUM_CELL; pos++)
  {
    if (pos == from)
    {
      continue;
    }
    if (excludeVisited && visited[pos])
    {
      continue;
    }

    
    double avgValue = (double)(-dp[pos]) / field->dist[from][pos] * pow(0.95, field->dist[from][pos]);
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
  {
    return this->findPathByHighValue(to);
  }
  case Mode::NOT_VISIT:
  {
    return this->findPathByNotVisit(to);
  }
  default:
  {
    LOG_ERROR("Invalid Mode Argument at Line %d", __LINE__);
    std::exit(-1);
    break;
  }
  }
}

void TakahashikunCleanerNo2::execute(const Commands &commands, int limit=10) 
{
  typedef char Command;
  int nxPos;
  int cnt = 0;
  for (Command command : commands)
  {
    if(cnt == limit)
    {
      return;
    }
    cnt++;
    switch (command)
    {
    case 'L':
    {
      nxPos = curPos - 1;
      break;
    }
    case 'R':
    {
      nxPos = curPos + 1;
      break;
    }
    case 'U':
    {
      nxPos = curPos - field->N;
      break;
    }
    case 'D':
    {
      nxPos = curPos + field->N;
      break;
    }
    default:
    {
      LOG_ERROR("Invalid Command '%c' at Line %d", command, __LINE__);
      std::exit(-1);
      break;
    }
    }
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

/**
 * @brief 価値の高い最短パスを通るようにパスを提案する。ダイクストラで見つける。
 *
 * @param to 行き先
 * @return Commands
 */
Commands TakahashikunCleanerNo2::findPathByHighValue(int to) const
{
  Commands ret;
  const int NUM_CELL = field->N * field->N;
  int from = curPos;
  Vec<int> dp(NUM_CELL, INF); // TODO: check whether 32 bit or 64 bit
  Vec<int> pre(NUM_CELL, -1);
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
 */
Commands TakahashikunCleanerNo2::findPathByNotVisit(int to) const
{
  Commands ret;

  return ret;
}

inline int TakahashikunCleanerNo2::getRCIdx(int r, int c) const
{
  return field->N * r + c;
}

inline char TakahashikunCleanerNo2::getAction(int cur, int nx) const
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

/**
 * @brief Main App class
 *
 */
class App
{
public:
  Input initialInput;
  TakahashikunCleanerNo2 takahashi;
  std::shared_ptr<Field> field;
  Commands bestCommands;
  std::mt19937 rng;

  App() {}
  App(const Input &input);
  void init();
  void run();
  void show() const;
  void summary() const;
  long long calcScore() const;

private:
  inline int getRCIdx(int r, int c) const;
};

App::App(const Input &input) : initialInput(input)
{
  this->field = std::make_shared<Field>(input);
}

void App::init()
{
  this->takahashi.init(this->field);
  this->rng = std::mt19937(seedGen());
  LOG_INFO("Application Initialization done at %6.4f s", getTime(startTime));
  LOG_INFO("Takahashi-kun is Ready!");
}

void App::run()
{
  const int NUM_CELL = field->N * field->N;
  Vec<int> cellIDs(NUM_CELL - 1);
  std::iota(cellIDs.begin(), cellIDs.end(), 1); // fill from 1 because 0 is origin
  std::shuffle(cellIDs.begin(), cellIDs.end(), rng);
  cellIDs.emplace_back(0); // add origin to end of the vec;

  bool excludeVisited = true;
  while (!this->takahashi.allVisited())
  {
    Commands commands;
    if ((u32)rng() % 10 == 0)
    {
      commands = this->takahashi.findPath(excludeVisited);
    }
    else
    {
      commands = this->takahashi.findPath();
    }
    this->takahashi.execute(commands, 10000);
  }

  // return to origin
  Commands commands = this->takahashi.findPath(0, Mode::HIGH_VALUE);
  this->takahashi.execute(commands, 10000);

  this->bestCommands = this->takahashi.generateCommandsFromHistory();
}

void App::show() const
{
  std::cout << bestCommands << "\n";
  std::cout.flush();
}

void App::summary() const
{
  fprintf(stderr, "\n######### SUMMARY ########\n");
  fprintf(stderr, "Score       : %10lld\n", calcScore());
  fprintf(stderr, "Elapsed Time: %10.2f s\n", getTime(startTime));
  fprintf(stderr, "##########################\n");
}

long long App::calcScore() const
{
  typedef char Command;
  long long ret = 0;
  const int NUM_CELL = field->N * field->N;

  Vec<Vec<int>> timeStamp = this->takahashi.timeStamp;
  int curTurn = bestCommands.size();
  int curPos = 0;
  int length = bestCommands.size();
  long long deltaS = 0;
  long long totalS = 0;
  long long curS = 0;
  for (int cellID = 0; cellID < NUM_CELL; cellID++)
  {
    curS += (long long)field->dirty[cellID] * (curTurn - timeStamp[cellID].back());
    deltaS += field->dirty[cellID];
  }

  for (Command command : bestCommands)
  {
    totalS += curS;
    int nxPos;
    switch (command)
    {
    case 'R':
    {
      nxPos = curPos + 1;
      break;
    }
    case 'L':
    {
      nxPos = curPos - 1;
      break;
    }
    case 'U':
    {
      nxPos = curPos - field->N;
      break;
    }
    case 'D':
    {
      nxPos = curPos + field->N;
      break;
    }
    default:
    {
      LOG_ERROR("Invalid command '%c' at Line %d\n", command, __LINE__);
      std::exit(-1);
    }
    }
    curTurn++;
    curS += deltaS;
    curS -= (long long)field->dirty[nxPos] * (curTurn - timeStamp[nxPos].back());
    timeStamp[nxPos].emplace_back(curTurn);
    totalS += curS;
    curPos = nxPos;
  }

  return std::round((double)totalS / length / 2);
}

inline int App::getRCIdx(int r, int c) const
{
  return r * field->N + c;
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
#ifdef LOCAL
  app.summary();
#endif
  app.show();

  return 0;
}