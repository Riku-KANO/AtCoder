#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;
#define ll long long
#define rep(i, n) for(int i = 1; i < n; i++)

void solve() {
    int h1, h2, h3, w1, w2, w3;
    cin >> h1>>h2>>h3>>w1>>w2>>w3;
    int ans = 0;
    rep(i, 29)rep(j, 29)rep(k, 29)rep(l, 29) {
        if(i+j>=h1||i+k>=w1||j+l>=w2||k+l>=h2) continue;
        // printf("%d %d %d %d\n", i, j, k, l);
        int a[3][3];
        a[0][0]=i; a[0][1]=j; a[1][0]=k; a[1][1]=l;
        a[0][2]=h1-i-j;
        a[1][2]=h2-k-l;
        a[2][0]=w1-i-k;
        a[2][1]=w2-j-l;
        a[2][2]=h3-a[2][0]-a[2][1];
        bool flag=false; 
        if(min({a[0][2], a[1][2], a[2][0], a[2][1], a[2][2]})>0 && w3==a[0][2]+a[1][2]+a[2][2])ans++;
    }
    cout << ans << "\n";
}

int main() {
    solve();
    return 0;
}