#include<iostream>
using namespace std;
using ll = long long;
int main() {
    int N; cin >> N;
    int A[N + 1];
    for(int i = 0; i < N; ++i) cin >> A[i + 1];
    ll n_same = 0;
    ll ans = 0;
    for(int i = 1; i <= N; ++i) {
        if(A[i] == i) n_same++;
        else if(A[A[i]] == i) ans++;
    }
    ans /= 2;
    ans += n_same * (n_same-1) / 2;
    cout << ans << endl;
    return 0;
}