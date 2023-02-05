#include<iostream>
#include <queue>
#include <vector>
#include <cmath>
using namespace std;

const int INF = 1e9;

int N, M, D, K;
int U[3001], V[3001], W[3001];
int X[1001], Y[1001];
std::vector<std::vector<int>> dist(1001);
vector<vector<pair<int,int>>> G(1001);

std::vector<int> dijkstra(int start, const vector<vector<pair<int,int>>>& graph) {
  typedef pair<int, pair<int,int>> Edge;
  std::vector<int> dist(N, INF);
  dist[start] = 0;
  priority_queue<Edge, vector<Edge>, greater<Edge>> pq;
  pq.push({0, {-1, start}});
  while(!pq.empty()) {
    Edge e = pq.top(); pq.pop();
    int u = e.second.second;
    int cost = e.first;
    if(cost > dist[u]) continue;
    for(pair<int,int> p: graph[u]) {
      int v = p.first;
      int w = p.second;
      if(dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push({dist[v], {u, v}});
      }
    }
  }
  return dist;
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

int main() {
  load_input();
  for(int i = 0; i < N; i++) {
    dist[i] = dijkstra(i, G);
  }
  long long score = calc_assign();
  cout << score << endl;
  return 0;
}