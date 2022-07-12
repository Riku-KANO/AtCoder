#include<iostream>
#include<map>
#include<set>
using namespace std;
using ll = long long;

int main() {
    int N; cin >> N;
    map<int, int> prime;
    map<int, int> mp;
    for(int i = 0; i < N; i++) {
        int m; cin >> m;
        for(int j = 0; j < m; j++) {
            int p, e; cin >> p >> e;
            if(prime[p] < e)  mp[p]=i;
            if(prime[p] == e) mp[p]=-1;
            prime[p]=max(e, prime[p]);
        }
    }
    set<int> s;
    for(auto p: mp) {
        if(p.second == -1) continue;
        s.insert(p.second);
    }
    int ans;
    if(s.size() == N) ans = N;
    else ans = s.size() + 1;
    cout << ans << endl;
    return 0;
}