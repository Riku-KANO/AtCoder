#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <bits/stdc++.h>
#define pii std::pair<int,int>
using namespace std;

const int INF = 1e9;
const int NUM_NODE = 200;

int N, M, D, K;
int U[3001], V[3001], W[3001];
int X[1001], Y[1001];
std::vector<std::vector<int>> dist(1001);
vector<vector<pair<int,int>>> G(1001);

std::random_device engine;
std::mt19937 mt(engine());
clock_t cur_time, start_time;

std::vector<int> dijkstra(int start, const vector<vector<pair<int,int>>>& graph) {
  typedef pair<int, pair<int,int>> Edge;
  std::vector<int> d(N, INF);
  d[start] = 0;
  priority_queue<Edge, vector<Edge>, greater<Edge>> pq;
  pq.push({0, {-1, start}});
  while(!pq.empty()) {
    Edge e = pq.top(); pq.pop();
    int u = e.second.second;
    int cost = e.first;
    if(cost > d[u]) continue;
    for(pair<int,int> p: graph[u]) {
      int v = p.first;
      int w = p.second;
      if(d[v] > d[u] + w) {
        d[v] = d[u] + w;
        pq.push({d[v], {u, v}});
      }
    }
  }
  return d;
}

void load_input() {
  cin >> N >> M >> D >> K;
  for(int i = 0; i < M; i++) {
    int u, v, w;
    cin >> u >> v >> w;
    u--; v--;
    U[i] = u;
    V[i] = v;
    W[i] = w;
    G[u].push_back({v, w});
    G[v].push_back({u, w});
    assert(u < N and v < N);
  }
  for(int i = 0; i < N; i++) {
    cin >> X[i] >> Y[i];
  }
  cerr << "INPUT ENDED" << endl;
}



long long calc_assign() {
  std::vector<int> assign(M);
  long long ans = 0;
  for(int i = 0; i < M; i++) {
    cin >> assign[i];
  }
  for(int d = 1; d <= D; d++) {
    std::vector<std::vector<int>> distK(N);
    std::vector<std::vector<pair<int,int>>> graph(N);
    for(int i = 0; i < M; i++) {
      if(assign[i] == d) continue;
      int u = U[i];
      int v = V[i];
      int w = W[i];
      graph[u].push_back({v, w});
      graph[v].push_back({u, w});
    }
    for(int i = 0; i < N; i++) {
      distK[i] = dijkstra(i, graph);
      for(int j = 0; j < N; j++) ans += distK[i][j] - dist[i][j];
    } 
  }
  ans = std::round((double)ans / N / (N - 1) / D * 1e3);
  return ans;
}

long long calc_r2(const std::vector<int>& nodes) {
  cerr << "CALC" << endl;
  long long ret = 1e18;
  std::vector<pii> points;
  for(int nid: nodes) points.push_back({X[nid], Y[nid]});
  int n = points.size();
  auto pow2 = [](long long a) {
    return a * a;
  };  
  for(int i = 0; i < n; i++) {
    for(int j = i + 1; j < n; j++) {
      long long r2 = pow2(points[i].first-points[j].first) + pow2(points[i].second-points[j].second);
      ret = std::min(ret, r2);
    }
  }
  return ret;
}

double get_time() {
  cur_time = clock();
  return (double)(cur_time-start_time)/CLOCKS_PER_SEC;
}

std::vector<int> move_random_node(const std::vector<int>& nodes) {
  int n = nodes.size();
  std::vector<int> ret = nodes;
  std::vector<bool> checked(n, false);
  for(int nid: nodes) checked[nid] = true;
  cerr << n << endl;
  std::uniform_int_distribution<int> d1(0, n - 1);
  int idx = d1(mt);
  cerr << "OK1" << endl;
  int nid = nodes[idx];
  int maxIter = 10;
  int n_adj = G[nid].size();
  std::uniform_int_distribution<int> d2(0, n_adj-1);
  cerr << "OK2" << endl;
  assert(nid < N);
  for(pii nx: G[nid]) cerr << nx.first << " ";
  cerr << endl;
  while(maxIter--) {
    cerr << "OK3" << endl;
    int to_idx = d2(mt);
    int to = G[nid][to_idx].first;

    cerr << n_adj << " " << to_idx << " " << to << endl;
    cerr << "OK4" << endl;
    if(!checked[to]) {
      ret[idx] = to;
      cerr << "OK5" << endl;
      break;
    }
  }
  return ret;
}


int main() {
  start_time = clock();
  load_input();
  for(int i = 0; i < N; i++) {
    dist[i] = dijkstra(i, G);
  }

  std::vector<int> range;
  for(int i= 0; i < N; i++) range.push_back(i);
  std::vector<int> cur_nodes;
  long long best_score, cur_score;
  std::sample(range.begin(), range.end(), std::back_inserter(cur_nodes), NUM_NODE, mt);
  cur_score = calc_r2(cur_nodes);
  best_score = cur_score;
  std::vector<int> best_nodes = cur_nodes;
  cerr << "INIT NODE: " << endl;
  for(int i = 0; i < NUM_NODE; i++) cerr << cur_nodes[i] << " ";
  cerr << endl;
  cerr << "INIT SCORE: " << cur_score << endl;
  while(get_time() < 2.0) {
    cerr << "STA" << endl;
    for(int i = 0; i < N; i++) {
      cerr << "I: " << i << ", ";
      for(auto p: G[i]) cerr << p.first << " ";
      cerr << endl;
    }
    std::vector<int> next_nodes(NUM_NODE);
    next_nodes = move_random_node(cur_nodes);
    for(int nid: next_nodes) cerr << nid << " ";
    cerr << endl;
    long long next_score = calc_r2(next_nodes);
    cerr << "NEXT: "<< next_score << endl;
    if(next_score > best_score) {
      cerr << "Hi" << endl;
      best_nodes = next_nodes;
      best_score = next_score;
      cerr << best_score << endl;
    } else {
      if(mt() % 10 < 9) {
        cerr << "Oopps" << endl;
        cur_nodes = next_nodes;
      } 
    }
  }
  return 0;
}