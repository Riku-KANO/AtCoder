#include <iostream>
#include <map>
#include <vector>

using ll = long long;

signed main() {
  int N;
  std::cin >> N;
  std::vector<std::vector<int>> cell(N,std::vector<int>(22));
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      std::cin >> cell[i][j];
    }
  }
  
  std::vector<std::vector<std::map<int,ll>>> dp1(N, std::vector<std::map<int,ll>>(N));
  std::vector<std::vector<std::map<int,ll>>> dp2(N, std::vector<std::map<int,ll>>(N));
  ll ans = 0;
  dp1[0][0][cell[0][0]] = 1;
  dp2[N-1][N-1][cell[N-1][N-1]] = 1;

  for(int i = 0; i < N; i++) {
    for(int j = 0; i + j <= N - 1; j++) {
      if(i > 0) {
        for(auto p: dp1[i-1][j]) {
          int v = p.first ^ cell[i][j];
          dp1[i][j][v] += p.second;
        }
      }
      if(j > 0) {
        for(auto p: dp1[i][j-1]) {
          int v = p.first ^ cell[i][j];
          dp1[i][j][v] += p.second;
        }
      }
    }
  }
  for(int i = N-1; i >= 0; i--) {
    for(int j = N-1; i + j >= N - 1; j--) {
      if(i < N - 1) {
        if(i + j == N - 1) {
          for(auto p: dp1[i][j]) {
            int v = p.first;
            ans += p.second * dp2[i+1][j][v];
          }
        } else {
          for(auto p: dp2[i+1][j]) {
            int v = p.first ^ cell[i][j];
            dp2[i][j][v] += p.second;
          }
        }
      }
      if(j < N - 1) {
        if(i + j == N - 1) {
          for(auto p: dp1[i][j]) {
            int v = p.first;
            ans += p.second * dp2[i][j+1][v];
          }
        } else {
          for(auto p: dp2[i][j+1]) {
            int v = p.first ^ cell[i][j];
            dp2[i][j][v] += p.second;
          }
        }
      }
    }
  }
  std::cout << ans << std::endl;
  return 0;
}