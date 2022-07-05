#include<bits/stdc++.h>
using namespace std;
using ll = long long;
#define rep(i, n) for(int i = 0; i < n; i++)

const int INF=1<<29;
vector<vector<int>> G(300'001);

void solve() {
    int N, M; cin >> N >> M;
    rep(i, M) {
        int u, v; cin >> u >> v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    auto bfs=[&](int start) {
        queue<int> q;
        vector<int> ret(N+1, INF);
        ret[start] = 0;
        q.push(start);
        while(!q.empty()){
            int st = q.front(); q.pop();
            for(auto nx: G[st]) {
                if(ret[nx] > ret[st]+1) {
                    ret[nx]=ret[st]+1;
                    q.push(nx);
                }
            }
        }
        return ret;
    };
    vector<int> d1=bfs(1);
    vector<int> dN=bfs(N);
    for(int i = 1; i <= N; i++) {
        ll ans;
        ans = min({d1[N], d1[i]+dN[0], d1[0]+dN[i]});
        if(ans>=INF) cout << -1 << endl;
        else cout << ans << endl;
    }
}

int main () {
    solve();
    return 0;
}