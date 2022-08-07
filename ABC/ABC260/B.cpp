#include<iostream>
#include<algorithm>
using namespace std;
#define rep(i, n) for(int i = 0; i < n; ++i)

struct Applicant {
    int id;
    int m, e;
    bool pass = false;
};

int main() {
    int N, X, Y, Z; cin >> N >> X >> Y >> Z;
    Applicant a[N];
    auto f = [](Applicant &a, Applicant &b) {
        if(a.m == b.m) return a.id < b.id;
        else return a.m > b.m;
    };
    auto g = [](Applicant &a, Applicant &b) {
        if(a.e == b.e) return a.id < b.id;
        else return a.e > b.e;
    };
    auto h = [](Applicant &a, Applicant &b) {
        if(a.m + a.e == b.m + b.e) return a.id < b.id;
        else return a.m + a.e > b.m + b.e;
    };
    auto k = [](Applicant &a, Applicant &b) {
        return a.id < b.id;
    };
    rep(i, N) cin >> a[i].m;
    rep(i, N) {
        cin >> a[i].e;
        a[i].id = i + 1;
    }
    sort(a, a + N, f);
    for(int i = 0; i < N && X; ++i) {
        if(a[i].pass) continue;
        else {
            a[i].pass = true;
            X--;
        }
    }
    sort(a, a + N, g);
    for(int i = 0; i < N && Y; ++i) {
        if(a[i].pass) continue;
        else {
            a[i].pass = true;
            Y--;
        }
    }
    sort(a, a + N, h);
    for(int i = 0; i < N && Z; ++i) {
        if(a[i].pass) continue;
        else {
            a[i].pass = true;
            Z--;
        }
    }
    sort(a, a + N, k);
    rep(i, N) {
        if(a[i].pass) cout << a[i].id << "\n";
    }
    return 0;
}