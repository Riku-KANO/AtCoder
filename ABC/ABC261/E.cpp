#include<iostream>
using namespace std;

int main() {
    int N, C; cin >> N >> C;
    int X = C;
    int M = (1<<30)-1;
    int s0 = 0;
    int s1 = M;
    for(int i = 0; i < N; ++i) {
        int t, a; cin >> t >> a;
        if(t == 1) s0&=a, s1 &= a;
        else if(t == 2) s0|= a, s1 |= a;
        else if(t == 3) s0 ^=a, s1 ^= a;
        X = (X&s1)|((~X)&s0);
        cout << X << "\n";
    }
    return 0;
}