#include<bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    int N, M; cin >> N >> M;
    map<ll, ll> mp;
    ll A[N];
    ll S[N-1];
    ll X[M];
    A[0] = 0;
    for(int i = 0; i < N-1; ++i) {
        cin >> S[i];
        A[i+1]=S[i]-A[i];
    }
    for(int i = 0; i < M; ++i) cin >> X[i];
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            ll d = A[i]-X[j];
            if(i%2)mp[d]++;
            else mp[-d]++;
        }
    }
    ll ans = -1;
    for(auto p: mp) ans=max(ans, p.second);
    cout << ans << endl;
    return 0;
}