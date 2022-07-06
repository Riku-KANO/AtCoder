#include<bits/stdc++.h>
#include<atcoder/all>
using namespace std;
using namespace atcoder;
using ll = long long;
void solve() {
  int N; cin >> N;
  scc_graph graph(N+1);
  int X[N+1];
  ll C[N+1];
  for(int i = 1; i <= N; ++i) {
    cin >> X[i];
    graph.add_edge(i, X[i]);
  }
  for(int i = 1; i <= N; ++i) cin >> C[i];
  ll ans = 0;
  for(auto s: graph.scc()) {
	if(s.size()==1) continue;
    ll tmp = 1LL<<60;
    for(auto v: s)tmp = min(tmp, C[v]);
    ans += tmp;
  }
  cout << ans << endl;
}

int main() {
  solve();
  return 0;
}