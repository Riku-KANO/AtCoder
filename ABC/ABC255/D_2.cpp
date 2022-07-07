#include<bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    ll N, Q; cin >> N >> Q;
    ll A[N];
    ll cum[N + 1] = {};
    for(int i = 0; i < N; ++i) cin >> A[i];
    sort(A, A + N);
    for(int i = 0; i < N; i++) cum[i+1]=cum[i]+A[i];
    while(Q--) {
        ll X; cin >> X;
        ll d = lower_bound(A, A + N, X)-A;
        ll ans = X*d-cum[d] + (cum[N]-cum[d]) - X * (N-d);
        cout << ans << endl;
    }
    return 0;
}
