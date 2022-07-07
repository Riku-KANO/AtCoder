#include<bits/stdc++.h>
using namespace std;

struct Position {
    int x, y;
};

int main() {
    int N, K; cin >> N >> K;
    cout << fixed << setprecision(15);
    double ans=-1;
    int a[K];
    Position pos[N];
    for(int i = 0; i < K; ++i) cin >> a[i], a[i]--;
    for(int i = 0; i < N; i++) cin >> pos[i].x >> pos[i].y;
    for(int i = 0; i < N; ++i) {
        double tmp = 1e20;
        for(int j = 0; j < K; ++j) {
            double dx = pos[i].x-pos[a[j]].x;
            double dy = pos[i].y-pos[a[j]].y;
            double r = sqrt(dx * dx + dy * dy);
            tmp = min(r, tmp);
        }
        ans = max(ans, tmp);
    }
    cout << ans << "\n";
    return 0;
}