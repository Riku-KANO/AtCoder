#include <bits/stdc++.h>
using namespace std;
using ll=long long;
#define rep(i, n) for(int i = 0; i < n; i ++)

int main(){
    int N, Q, X; cin >> N >> Q >> X;
    int W[N];
    ll cum[N+1]={};
    rep(i, N) {
        cin >> W[i];
        cum[i+1]=cum[i]+W[i];
    }
    int nxt[N];
    ll cycle[N];
    for (int i = 0; i < N; i++){
        ll a = cum[i] + X;
        cycle[i] = a / cum[N];
        a %= cum[N];
        nxt[i] = lower_bound(cum, cum+N+1, a)-cum;
        if (nxt[i] == N){
            cycle[i]++;
            nxt[i] = 0;
        }
    }
    int bl[50][N];
    rep(i, N) bl[0][i] = nxt[i];
    for (int i = 0; i < 49; i++){
        for (int j = 0; j < N; j++){
            bl[i + 1][j] = bl[i][bl[i][j]];
        }
    }

    while(Q--) {
        ll K; cin >> K;
        K--;
        int p = 0;
        for (int j = 0; j < 50; j++) if((K >> j & 1) == 1) p = bl[j][p];
        ll ans = cycle[p] * N + nxt[p] - p;
        cout << ans << endl;
    }
}