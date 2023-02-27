#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

using ll = long long;

std::vector<int> v;

// 1-20, 21-40, ...., 4981-5000
ll hist[11][255][255];
int N, W, K, C;
int S[204][204];
int L[254][254];
int label;

bool is_out(int i, int j, int N) {
  return (i < 0 or j < 0 or i >= N or j >= N);
}



bool solve() {
  std::cin >> N >> W >> K >> C;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      std::cin >> S[i][j];
      assert(S[i][j]>=10 and S[i][j] <= 5000);
      label = std::lower_bound(v.begin(), v.end(), S[i][j]) - v.begin();
      L[i][j] = label;
      assert(label < 250 and label >= 0);
    }
  }
  int hoge;
  for(int i = 0; i < W; i++) std::cin >> hoge >> hoge;
  for(int i = 0; i < K; i++) std::cin >> hoge >> hoge;
  int d;

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      for(int r = -10; r <= 10; r++) {
        for(int c = std::abs(r)-10; c <= 10-std::abs(r); c++) {
          if((r == 0 and c == 0)) continue;
          if(is_out(i+r, j+c,N)) continue;
          d = std::abs(r) + std::abs(c);
          hist[d][L[i][j]][L[i+r][j+c]]++;
        }
      }
    }
  }
  return false;
}

int main() {
  // init
  for(int d = 0; d <= 10; d++) {
    for(int i = 0; i <= 250; i++) {
      for(int j = 0; j <= 250; j++) {
        hist[d][i][j] = 0;
      }
    }
  }
  for(int i = 1; i <= 250; i++) v.push_back(20 * i);

  // solve
  for(int i = 0; i < 5000; i++) {
    if((i + 1) % 50 == 0) {
      std::cerr << i + 1 << " files done." << std::endl;
    }
    if(solve()) {
      std::cerr << "file " << i << " has problem" << std::endl;
      return 0;
    }
  }


  // output
  for(int d = 1; d <= 10; d++) {
    std::cout << "dist " << d << std::endl;
    for(int i = 0; i < 250; i++) {
      for(int j = 0; j < 250; j++) {
        std::cout << hist[d][i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
  return 0;
}