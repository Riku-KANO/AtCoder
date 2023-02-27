#include <iostream>
#include <vector>
#include <queue>

using namespace std;

bool is_out(int ni, int nj, int N) {
  return (ni < 0 or nj < 0 or ni >= N or nj >= N);
}

static const int di[4] = {0, 1, 0, -1};
static const int dj[4] = {1, 0, -1, 0};

int N, W, K, C;
int wi[5], wj[5];
int hi[10], hj[10];
int S[404][404];
int total_cost = 0;
std::vector<std::vector<bool>> road(404, std::vector<bool>(404, false));

int main() {
  std::cerr << "Judge start" << std::endl;
  std::cin >> N >> W >> K >> C;
  for(int i = 0; i < N; i++) for(int j = 0; j < N; j++) std::cin >> S[i][j];
  for(int i = 0; i < W; i++) std::cin >> wi[i] >> wj[i];
  for(int i = 0; i < K; i++) std::cin >> hi[i] >> hj[i];

  std::cout << N << " " << W << " " << K << " " << C << std::endl;
  for(int i = 0; i < W; i++) std::cout << wi[i] << " " << wj[i] << std::endl;
  for(int i = 0; i < K; i++) std::cout << hi[i] << " " << hj[i] << std::endl;

  
  std::vector<std::vector<bool>> visited(404, std::vector<bool>(404, false));
  while(true) {
    int i, j;
    int power;
    std::cin >> i >> j >> power;
    if(i < 0 or j < 0 or i >= N or j >= N or power < 0 or power > 5000) {
      std::cerr << "Invalid output!!" << std::endl;
      std::cerr << i << " " << j << " " << power << std::endl;
      std::cout << -1 << std::endl;
      return 0;
    } else if(S[i][j] <= 0) {
      std::cerr << "position (" << i << "," << j << ") already excavated." << std::endl;
      std::cout << -1 << std::endl;
      return 0;
    }
    for(int w = 0; w < W; w++) {
      if(wi[w] == i and wj[w] == j) {
        std::cerr << "Invalid Position" << std::endl;
        std::cout << -1 << std::endl;
      }
    }

    S[i][j] -= power;
    total_cost += (power + C);
    if(S[i][j] > 0) {
      std::cout << 0 << std::endl;
    } else {
      road[i][j] = true;

      // connected check
      for(int w = 0; w < W; w++) {
        std::queue<std::pair<int,int>> q;
        q.push({wi[w], wj[w]});
        while(!q.empty()) {
          auto p = q.front();
          int ci = p.first;
          int cj = p.second;
          if(visited[ci][cj]) continue;
          visited[ci][cj] = true;
          for(int d = 0; d < 4; d++) {
            int ni = ci + di[d];
            int nj = cj + dj[d];
            if(is_out(ni, nj, N)) continue;
            if(visited[ni][nj] or !road[ni][nj]) continue;
            q.push({ni, nj});
          }
        }
      }

      bool connected_flag = true;
      for(int k = 0; k < K; k++) {
        if(!visited[hi[k]][hj[k]]) connected_flag = false;
      }

      if(connected_flag) {
        std::cout << 2 << std::endl;
        std::cerr << "SCORE: " << total_cost << std::endl;
        return 0;
      } else {
        std::cout << 1 << std::endl;
      }
      visited.clear();
      visited.resize(N, std::vector<bool>(N, false));
    }
  }
  return 0;
}