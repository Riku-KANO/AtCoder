#include <iostream>
using ll = long long;

ll cnt[5001];
int N, W, K, C;
int hoge;

void solve() {
  std::cin >> N >> W >> K >> C;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      std::cin >> hoge;
      cnt[hoge]++;
    }
  }
  for(int i = 0; i < W; i++) std::cin >> hoge >> hoge;
  for(int i = 0; i < K; i++) std::cin >> hoge >> hoge;
}


int main() {
  for(int i = 0; i <=5000; i++) cnt[i] = 0;


  for(int i = 0; i < 5000; i++) {
    if((i + 1) % 50 == 0) {
      std::cerr << i + 1 << " file done." << std::endl; 
    }
    solve();
  }

  for(int i = 0; i <= 5000; i++) {
    std::cout << cnt[i] << std::endl;
  }
  return 0;
}