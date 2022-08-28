#include<iostream>
using namespace std;
using ll = long long;

int main() {
  ll n; cin >> n;
  ll MOD = 998244353;
  cout << (n%MOD+MOD)%MOD << "\n";
  return 0;
}