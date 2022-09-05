#include<iostream>
#include<vector>
#include<algorithm>
#include<atcoder/all>
using namespace std;
using namespace atcoder;
#define rep(i, n) for(int i = 0; i < n; ++i)

vector<pair<int,int>> wire(500006);

void solve() {
	int n, m, e; cin >> n >> m >> e;
	// 200001発電所の親
	dsu d(200005);
	for(int i = n; i < n+m; ++i) {
		d.merge(200001, i);
	}
	rep(i, e) {
		int u, v; cin >> u >> v;
		u--; v--;
		wire[i]={u,v};
	}
	int q; cin >> q;
	int x[q], xt[q];
	for(int i = 0; i < q; ++i){
		cin >> x[i]; x[i]--;
		xt[i]=x[i];
	}
	sort(xt, xt+q);
	int tmp=0;
	for(int i = 0; i < e; ++i) {
		if(i==xt[tmp]) {
			tmp++;
		} else {
			int u = wire[i].first;
			int v = wire[i].second;
			d.merge(u, v);
		}
	}
	
	vector<int> ans(q);
	for(int i = 0; i < q; ++i) {
		int wire_id =x[q-1-i];
		ans[q-1-i]=d.size(200001)-1-m;
		d.merge(wire[wire_id].first, wire[wire_id].second);
	}

	for(int i = 0; i < q; ++i) cout << ans[i] << "\n";
}

int main(void){
  solve();
  return 0;
}
