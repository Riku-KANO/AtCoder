#include <iostream>
#include <vector>

using ll = long long;

struct Road {
  int id;
  int from, to;
  int len;
};

template<typename T> void chmin(T& a, T b) {a = std::min(a, b);}

bool operator<(const Road& lhs, const Road& rhs) {
  return lhs.len < rhs.len;
}
bool operator>(const Road& lhs, const Road& rhs) {
  return lhs.len > rhs.len;
}

static const int MAXN = 2e5+5;
static const ll LINF = 1e18 + 10;

std::vector<std::vector<Road>> G(MAXN);

signed main() {
  int N, M, K;
  std::cin >> N >> M >> K;
  std::vector<int> E(K);
  std::vector<Road> roads(M + 1); 
  for(int i = 1; i <= M; i++) {
    int a, b, c;
    std::cin >> a >> b >> c;
    roads[i].id = i;
    roads[i].from = a;
    roads[i].to = b;
    roads[i].len = c;
    G[a].push_back(roads[i]);
  }
  std::vector<ll> dist(N + 1, LINF);
  dist[1] = 0;
  for(int i = 0; i < K; i++) {
    int e;
    std::cin >> e;
    if(dist[roads[e].from] != LINF) {
      chmin(dist[roads[e].to], dist[roads[e].from]+roads[e].len);
    }
  }
  std::cout << ((dist[N]==LINF)? -1 : dist[N]) << std::endl;
  return 0;
}