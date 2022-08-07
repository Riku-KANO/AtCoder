#include<iostream>
using namespace std;
using ll = long long;

ll dpr[10]={}, dpb[10]={};

int main() {
    int N, X, Y; cin >> N >> X >> Y;
    dpr[N]=1;
    for(int i = N; i >= 2; --i) {
        dpr[i-1] += dpr[i];
        dpb[i] += dpr[i]*X;

        dpr[i-1] += dpb[i];
        dpb[i-1] += dpb[i]*Y;
    }
    cout << dpb[1] << endl;
    return 0;
}