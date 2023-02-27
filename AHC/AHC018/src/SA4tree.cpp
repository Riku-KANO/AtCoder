//baseline
#include <iostream>
#include <vector>
#include <algorithm>
#include <atcoder/dsu>
#define pii std::pair<int,int>
using namespace std;
using namespace atcoder;

int N, W, K, C;
int wi[5];
int wj[5];
int hi[15];
int hj[15];
int pi[20];
int pj[20];



long long pow2(int x, int y){
  return (long long)(x-y)*(x-y);
}


// normal excavation
int excavation(int x, int y, int power) {
  assert(10 <= power <= 5000);
  std::cout << x << " " << y << " " << power << std::endl;
  int res;
  std::cin >> res;
  return res;
}

int main(int argc, char* argc[]) {
  std::cin >> N >> W >> K >> C;
  int hoge;
  for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) std::cin >> hoge;
  for(int i = 0; i < W; i++) {
    std::cin >> wi[i] >> wj[i];
    pi[i] = wi[i];
    pj[i] = wj[i];
  }
  for(int i = 0; i < K; i++) {
    std::cin >> hi[i] >> hj[i];
    pi[i+W] = hi[i];
    pj[i+W] = hj[i];
  }
  std::vector<std::pair<long long, pii>> dists;
  for(int i = 0; i < W + K; i++) {
    for(int j = i + 1; j < W + K; j++) {
      long long dist = pow2(pi[i], pi[j]) + pow2(pj[i], pj[j]);
      dists.push_back({dist, pii{i, j}});
    }
  }
  std::sort(dists.begin(), dists.end());
  std::vector<std::vector<int>> G(20);
  std::vector<pii> edges;
  std::vector<int> boss(20, -1);
  dsu uf(W+K);
  int cnt = 0;
  for(auto p: dists) {
    long long dist = p.first;
    int i = p.second.first;
    int j = p.second.second;
    if(uf.same(i, j) or (i < W and j < W)) continue;
    if(i < W and boss[j] != -1) continue;
    if(boss[i] != -1 and boss[j] != -1 and boss[i] != boss[j]) continue;
    G[i].push_back(j);
    G[j].push_back(i);
    uf.merge(i, j);
    cnt++;
    edges.push_back({i, j});
    for(int w = 0; w < W; w++) {
      for(int k = 0; k < K; k++) if(uf.same(w, k+W)) boss[W+k] = w;
    }
  }
  int M = edges.size();
  cerr << M << endl;
  for(auto e: edges) {
    std::cerr << e.first << " " << e.second << std::endl;
  }
  long long L = 1e18;
  std::vector<int> ans(M);
  std::vector<std::vector<int>> ans_grid;
  auto fill = [](int si, int sj, int ti, int tj, std::vector<std::vector<int>>& v) {
    if(si == ti) for(int j = std::min(sj, tj); j <= std::max(sj, tj); j++) v[si][j] = 1;
    else for(int i = std::min(si, ti); i <= std::max(si, ti); i++) v[i][sj] = 1;
  };
  std::vector<std::vector<int>> grid(N, std::vector<int>(N, 0));
  for(int i = 0; i < (1<<M); i++) {
    grid.clear();
    grid.resize(N, std::vector<int>(N, 0));
    long long l = 0;
    for(int j = 0; j < M; j++) {
      int s = edges[j].first;
      int t = edges[j].second;
      
      int psi, psj;
      int pti, ptj;
      if(s < W) psi = pi[s], psj = pj[s];
      else psi = pi[s], psj = pj[s];
      if(t < W) pti = pi[t], ptj = pj[t];
      else pti = pi[t], ptj = pj[t];
      // std::cerr << "("<<psx << ","<<psy<<"), ("<<ptx<<","<<pty<<")" << std::endl;
      if((1<<j)&i) {
        int pui = pti;
        int puj = psj;
        fill(psi, psj, pui, puj, grid);
        fill(pui, puj, pti, ptj, grid);
      } else {
        int pui = psi;
        int puj = ptj;
        fill(psi, psj, pui, puj, grid);
        fill(pui, puj, pti, ptj, grid);
      }
    }
    for(int r = 0; r < N; r++) for(int c = 0; c < N; c++) if(grid[r][c]==1) l++;
    if(l < L) {
      L = l;
      ans_grid = grid;
      for(int j = 0; j < M; j++) ans[j] = i&(1<<j) ? 1 : 0;
    }
  }
  std::cerr << "L: " << L << std::endl;


  int power = C * 5;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      if(ans_grid[i][j] == 1) {
        for(int w = 0; w < W; w++) if(i == wi[w] and j == wj[w]) continue;
        int res = -1;
        while(res != 1) {
          // res = excavation(i, j, power);
          if(res == -1) {
            std::cerr << "Invalid output!" << std::endl;
            std::exit(1);
          } else if(res == 2) {
            std::cerr << "All done!!" << std::endl;
            return 0;
          }
        }
      }
    }
  }
  return 0;
}