#include<bits/stdc++.h>
using namespace std;
using ll = long long;
#define rep(i,n) for(int i = 0; i < n; ++i)

struct Jumper {
    ll x, y;
    ll p;
};

void solve() {
    int N; cin >> N;
    Jumper J[N];
    rep(i, N) cin >> J[i].x >> J[i].y >> J[i].p;
    ll ls = 0, rs = 4e9;
    ll ms;
    while(ls + 1 != rs) {
        ms = (ls + rs) / 2;
        int dist[200][200] = {};
        rep(i,N)dist[i][i]=1;
        rep(i, N) rep(j, N) {
            if(i == j) continue;
            ll dx = abs(J[i].x-J[j].x);
            ll dy = abs(J[i].y-J[j].y);
            if(ms * J[i].p >= dx + dy) dist[i][j]=1;;
        }
        
        rep(k, N)rep(i, N)rep(j, N) dist[i][j] |= dist[i][k]&dist[k][j]; 
        bool flag=false;
        rep(i, N) {
            int sum = 0;
            rep(j, N) if(dist[i][j])sum++;
            if(sum==N) {
                flag = true;
                break;
            }
        }
        if(flag) rs=ms;
        else ls=ms;
    }
    cout << rs << endl;
}

int main() {
    solve();
    return 0;
}