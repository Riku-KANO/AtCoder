#include<iostream>
using namespace std;
using ll = long long;
ll dp[5005][5005] = {};
int main(){
    int N, M; cin >> N >> M;
    ll X[N];
    ll B[5005] = {};
    for(int i = 0; i < N; ++i) cin >> X[i];
    for(int i = 0; i < M; ++i) {
        int C, Y; cin >> C >> Y;
        B[C] = Y;
    }
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j <= i; ++j) {
            dp[i+1][j+1]=max(dp[i+1][j+1], dp[i][j]+X[i]+B[j+1]); //表が出る
            dp[i+1][0]=max(dp[i+1][0], dp[i][j]); //裏が出る
        }
    }
    ll ans = 0;
    for(int i = 0; i <= N; ++i) {
        ans = max(ans, dp[N][i]);
    }
    cout << ans << "\n";
    return 0;
}