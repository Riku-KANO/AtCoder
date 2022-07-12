#include<bits/stdc++.h>
using namespace std;
using ll = long long;

struct Point {
    ll x, y;
};

struct Circle {
    Point p;
    ll r;
};

vector<vector<int>> G(3005);
bool visited[3005];

bool dfs(int start, int goal) {
    if(start == goal) return true;
    bool ret = false;
    if(visited[start]) return false;
    visited[start] = true;
    for(int nx: G[start]) ret |= dfs(nx, goal);
    return ret;
}

int main() {
    int N; cin >> N;
    Point s, t;
    cin >> s.x >> s.y >> t.x >> t.y;
    Circle crc[N];
    int sc, tc;
    auto pow2 = [](ll x) {return x * x;};
    for(int i = 0; i < N; i++) {
        cin >> crc[i].p.x >> crc[i].p.y >> crc[i].r;
        ll dxs = s.x-crc[i].p.x;
        ll dys = s.y-crc[i].p.y; 
        ll dxt = t.x-crc[i].p.x;
        ll dyt = t.y-crc[i].p.y;
        if(pow2(crc[i].r) == pow2(dxs) + pow2(dys)) sc = i;
        if(pow2(crc[i].r) == pow2(dxt) + pow2(dyt)) tc = i;
    }
    for(int i = 0; i < N; i++) {
        for(int j = i + 1; j < N; j++) {
            ll dx = crc[i].p.x - crc[j].p.x;
            ll dy = crc[i].p.y - crc[j].p.y;
            ll d2 = dx * dx + dy * dy;
            if(d2 <= pow2(crc[i].r + crc[j].r) && d2 >= pow2(crc[i].r-crc[j].r)) {
                G[i].push_back(j);
                G[j].push_back(i);
            }
        }
    }
    memset(visited, false, sizeof(visited));
    bool res = dfs(sc, tc);
    if(res) puts("Yes");
    else puts("No");
    return 0;
}