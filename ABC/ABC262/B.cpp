#include<iostream>
using namespace std;

int main() {
    int N, M; cin >> N >> M;
    int A[N + 1][N + 1] = {};
    for(int i = 0; i < M; ++i) {
        int u, v; cin >> u >> v;
        A[u][v] = A[v][u] = 1;
    }
    int ans = 0;
    for(int i = 1; i <= N; ++i) {
        for(int j = i + 1; j <= N; ++j) {
            for(int k = j + 1; k <= N; ++k) {
                if(A[i][j]&A[j][k]&A[k][i]) ans++;
            }
        }
    }
    cout << ans << "\n";
    return 0;
}