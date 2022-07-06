#include<iostream>
using namespace std;
using ll = long long;
#define rep(i, n) for(int i = 0; i < n; i++) 

void solve() {
    int N; cin >> N;
    int A[200003] = {};
    rep(i, N) {
        int L, R; cin >> L >> R;
        A[L]++;
        A[R]--;
    }
    rep(i, 200001) A[i+1]+=A[i];
    bool flag = false;
    rep(i, 200'001) {
        if(A[i]>0 && !flag) {
            flag = true;
            cout << i << " ";
        } else if(A[i]==0&&flag) {
            cout << i << "\n";
            flag=false;
        }
    }
}

int main() {
    solve();
    return 0;
}