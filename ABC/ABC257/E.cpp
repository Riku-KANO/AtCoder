#include<bits/stdc++.h>
using namespace std;
#define rep(i, n) for(int i = 0; i < n; i ++)

void solve() {
    int N; cin >> N;
    int C[10];
    int d;
    int m=1e7;
    rep(i, 9) {
        cin >> C[i+1];
        if(m>=C[i+1]) {
            d = i+1;
            m=C[i+1];
        }
    }
    int max_digit = N/m;
    string ans = "";
    int remain = N;
    rep(i, max_digit) {
        for(int c = 9; c >= 1; c--) {
            if((remain-C[c])/m>=max_digit-1-i&&remain-C[c]>=0) {
                ans+= c+'0';
                remain-=C[c];
                break;
            }
        }
    }
    cout << ans << endl;
}

int main() {
    solve();
    return 0;
}