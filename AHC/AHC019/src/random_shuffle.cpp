#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>

std::random_device seed_gen;
std::mt19937 mt(seed_gen());

std::vector<int> gen_shuffle_arr(int n) {
  std::vector<int> ret(n);
  std::iota(ret.begin(), ret.end(), 0);
  std::shuffle(ret.begin(), ret.end(), mt);
  return ret;
}

int main() {
  for(int v: gen_shuffle_arr(10)) {
    std::cout << v << " ";
  }
  std::cout << "\n";
  return 0;
}