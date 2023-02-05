#include<iostream>
#include<random>
#include<vector>
#include<utility>
#include<queue>
#include<algorithm>
#include<cmath>

#define pii std::pair<int,int>

clock_t start_time, cur_time;
std::random_device seed_gen;
std::mt19937 mt(seed_gen());

const double TL = 6.0;
const int INF = 1e9;

int N, M, D, K;
int U[3001], V[3001], W[3001];
int X[1001], Y[1001];
std::vector<std::vector<int>> min_dist(1001);
std::vector<std::vector<std::vector<int>>> dist(31, std::vector<std::vector<int>>(1001));

std::vector<std::vector<std::pair<int,int>>> G(1001);

struct Result {
  std::vector<int> plan;
  int score;
  Result() {
    this->score = 1e9;
  }
  void show() {
    for(int p: this->plan) std::cout << p << "\n";
  }
};

struct Edge {
  int u, v, dist;
};

bool operator<(const Edge& lhs, const Edge&rhs) {
  return lhs.dist < rhs.dist;
}

bool operator>(const Edge& lhs, const Edge&rhs) {
  return lhs.dist > rhs.dist;
}

std::vector<int> dijkstra(int start, const std::vector<std::vector<std::pair<int,int>>> &graph)
{

    int n = N;
    std::vector<int> dist(n + 1, INF);
    dist[start] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{-1, start, 0LL});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.v;
        int cost = edge.dist;
        if(dist[from] != cost) continue;
        for (std::pair<int, int> ne : graph[from])
        {
            int to = ne.first;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pq.push(Edge{from, to, dist[to]});
            }
        }
    }
    return dist;
}

std::vector<int> make_plan_at_random() {
  std::vector<int> plan(M);
  std::vector<int> remain(D, K);
  std::vector<int> cum(D + 1, 0);

  for(int i = 0; i < M; i++) {
    for(int j = 0; j < D; j++) cum[j+1] = cum[j]+remain[j];
    int v = mt() % cum[D];
    int d = std::upper_bound(cum.begin(), cum.end(), v) - cum.begin();
    plan[i] = d;
  }

  return plan;
}

void calc_dist(std::vector<int> plan) {
  std::vector<std::vector<std::pair<int,int>>> graph = G;
  for(int i = 0; i < D; i++) {
    std::cerr << "DAY: " << i + 1 << std::endl;
    std::vector<int> buf;
    for(int j = 0; j < M; j++) {
      if(plan[j] == i + 1) {
        int u = U[j];
        int v = V[j];
        int w = W[j];
        auto it = std::find(graph[u].begin(), graph[u].end(), pii{v, w});
        graph[u].erase(it);
        it = std::find(graph[v].begin(), graph[v].end(), pii{u, w});
        graph[v].erase(it);
        buf.push_back(j);
      }
    }
    for(int s = 0; s < N; s++) dist[i][s] = dijkstra(s, graph);
    for(int b: buf) {
      int u = U[b];
      int v = V[b];
      int w = W[b];
      graph[u].push_back({v, w});
      graph[v].push_back({u, w});
    }
    // graph.clear();
    // graph.resize(N);
  }
}

long long calc_plan_score() {
  std::vector<double> score_per_day(D);
  double ret = 0;
  for(int k = 0; k < D; k++) {
    for(int i = 0; i < N; i++) {
      for(int j = i + 1; j < N; j++) {
        score_per_day[k] += dist[k][i][j] - min_dist[i][j];
      }
    }
    score_per_day[k] /= N * (N - 1);
  }

  for(int i = 0; i < D; i++) ret += score_per_day[i];

  return (long long)std::round(ret * 1000 / D);
}

double get_time() {
  cur_time = clock();
  return (double) (cur_time - start_time) / CLOCKS_PER_SEC;
}

int main(int argv, char** argc) {
  start_time = clock();
  Result best_result;

  std::cin >> N >> M >> D >> K;

  for(int i = 0; i < M; i++) {
    std::cin >> U[i] >> V[i] >> W[i];
    U[i]--;
    V[i]--;
    G[U[i]].push_back({V[i], W[i]});
    G[V[i]].push_back({U[i], W[i]});
  }

  // calc min dist
  for(int i = 0; i < N; i++) min_dist[i] = dijkstra(i, G);
  


  for(int i = 0; i < N; i++) {
    std::cin >> X[i] >> Y[i];
  }
  int iter = 0;
  while(get_time() < TL * 0.9) {
    std::cerr << "ITERATION: " << ++iter << std::endl;
    std::vector<int> plan = make_plan_at_random();
    std::cerr << "calc" << std::endl;
    calc_dist(plan);
    long long score = calc_plan_score();

    if(score < best_result.score) {
      best_result.score = score;
      best_result.plan = plan;
    }
  }
  best_result.show();
  std::cerr << get_time() << std::endl;
  std::cerr << "SCORE: " << best_result.score << std::endl;
  return 0;
}