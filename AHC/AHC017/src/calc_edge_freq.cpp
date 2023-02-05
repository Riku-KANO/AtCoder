#include <iostream>
#include <queue>
#include <vector>
#include <tuple>
#include <cmath>
#include <utility>
#include <algorithm>

using namespace std;
#define pii std::pair<int,int>

const int INF = 1e9;

int N, M, D, K;
int U[3001], V[3001], W[3001];
int X[1001], Y[1001];
std::vector<std::vector<int>> dist(1001);
vector<vector<pair<int,int>>> G(1001);

bool operator==(const pii& lhs, const pii& rhs) {
  return ((lhs.first == rhs.first and lhs.second == rhs.second) or (lhs.first == rhs.second and lhs.second == rhs.first));
}

std::tuple<std::vector<int>, std::vector<int>> dijkstra(int start, const vector<vector<pair<int,int>>>& graph) {
  typedef pair<int, pair<int,int>> Edge;
  std::vector<int> dist(N, INF);
  std::vector<int> pre(N, -1);
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
        pre[v] = u;
        pq.push({dist[v], {u, v}});
      }
    }
  }
  return {dist, pre};
}

std::tuple<std::vector<int>, std::vector<int>> dijkstra(int s, int e1)
{
    typedef pair<int, pii> Edge;
    pii e1p = {U[e1], V[e1]};
    int n = N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);
    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[from])
            continue;

        for (pii ne : G[from])
        {
            int to = ne.first;
            if (pii{from, to} == e1p)
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist, pre};
}

std::tuple<std::vector<int>, std::vector<int>> dijkstra(int s, int e1, int e2)
{
    typedef pair<int, pii> Edge;
    pii e1p = {U[e1], V[e1]};
    pii e2p = {U[e2], V[e2]};
    int n = N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);
    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[from])
            continue;

        for (pii ne : G[from])
        {
            int to = ne.first;
            if (pii{from, to} == e1p or pii{from, to} == e2p)
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist, pre};
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


int main() {
  load_input();
  std::vector<pii> freq(M, {0, 0});
  std::vector<std::vector<int>> edge_id_matrix(N, std::vector<int>(N, -1));
  for(int i = 0; i < M; i++) {
    freq[i].second = i;
    int u = U[i];
    int v = V[i];
    edge_id_matrix[u][v] = edge_id_matrix[v][u] = i;
  }
  for(int i = 0; i < N; i++) {
    std::tuple<std::vector<int>, std::vector<int>> shortest = dijkstra(i, G);
    dist[i] = std::get<0>(shortest);
    std::vector<int> pre = std::get<1>(shortest);
    std::vector<bool> used(M, false);
    for(int j = 0; j < N; j++) {
      if(i == j) continue;
      int v = j;
      while(pre[v] != -1) {
        int u = pre[v];
        int eid = edge_id_matrix[u][v];
        v = u;
        used[eid] = true;
      }
    }
    for(int j = 0; j < M; j++) if(used[j]) freq[j].first++;
  }

  std::sort(freq.rbegin(), freq.rend());
  for(int i = 0; i < 10; i++) {
    int eid = freq[i].second;
    long long f = 0, ff = 0;
    std::vector<std::vector<int>> dist2(N);
    std::vector<pii> freq2(M, {0, 0});
    for(int j = 0; j < M; j++) freq2[j].second = j;
    for(int j = 0; j < N; j++) {
      std::vector<bool> used2(M, false);
      std::tuple<std::vector<int>, std::vector<int>> shortest2 = dijkstra(j, eid);
      dist2[j] = std::get<0>(shortest2);
      std::vector<int> pre2 = std::get<1>(shortest2);
      for(int k = 0; k < N; k++) f += dist2[j][k] - dist[j][k];
      for(int k = 0; k < N; k++) {
        if(j == k) continue;
        int v = k;
        while(pre2[v] != -1) {
          int u = pre2[v];
          int eid2 = edge_id_matrix[u][v];
          v = u;
          used2[eid2] = true;
        }
      }
      for(int k = 0; k < M; k++) if(used2[k]) freq2[k].first++;
    }
    f = std::round(1e3 * f / N / (N - 1) / D);
    cerr << "EID: "<< eid <<  ", INCREASED F: " << f << endl;


    std::sort(freq2.rbegin(), freq2.rend());
    int eid2 = freq2[0].second;
    for(int j = 0; j < N; j++) {
      std::tuple<std::vector<int>, std::vector<int>> shortest3 = dijkstra(j, eid, eid2);
      std::vector<int> dist3 = std::get<0>(shortest3);
      std::vector<int> pre3 = std::get<1>(shortest3);
      for(int k = 0; k < N; k++) ff += dist3[k] - dist[j][k];
    }
    ff = std::round(1e3 * ff / N / (N - 1) / D);
    std::cerr << "Freq EID2: " << freq2[0].first << endl;
    cerr << "EID: ("<< eid << "," << eid2 <<  "), INCREASED F: " << ff << endl;

  }
  // output
  std::sort(freq.begin(), freq.end(),
    [&](const pii&lhs, const pii&rhs) {return lhs.second < rhs.second;}
  );
  for(int i = 0; i < M; i++) {
    cout << freq[i].first << endl;
  }

  return 0;
}