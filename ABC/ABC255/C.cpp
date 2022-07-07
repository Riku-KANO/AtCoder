#include<bits/stdc++.h>
using namespace std;
using ll = long long;

int main() {
    ll X, A, D, N; cin >> X >> A >> D >>N;
    if(D<0) {
        X = -X;
        A = -A;
        D = -D;
    } else if(D==0) {
        cout << abs(X-A) << endl;
        return 0;
    }
    ll ans;
    ll last = A + (N-1)*D;
    if(X<=A) ans = A-X;
    else if(last<=X) ans = X-last;
    else {
        ans = min((X-A)%D, D-(X-A)%D);
    }
    cout << ans << endl;
    return 0;
}