#include<iostream>
#include<algorithm>
#include<vector>
#include<set>
#include<map>
using namespace std;

int main() {
    int N, K; cin >> N >> K;
    vector<int> parent(N + 1, -1);
    vector<int> ans(N + 1, -1);
    set<int> s;
    s.insert(-1);
    map<int, vector<int>> mp;
    for(int i = 0; i < N; ++i) {
        int p; cin >> p;
        auto it = s.lower_bound(p);
        if(it==s.end()) {
            s.insert(p);
            parent[p] = p;
            mp[parent[p]].push_back(p);
        } else {
            parent[p] = parent[*it];
            mp[parent[p]].push_back(p);
            s.erase(it);
            s.insert(p);
            cout << p << endl;
        }

        if(mp[parent[p]].size()==K) {
            for(int num: mp[parent[p]]) ans[num] = i+1;
            mp.erase(parent[p]);
            it = s.lower_bound(p);
            s.erase(it);
        }
    }
    for(int i = 1; i <= N; ++i) {
        cout << ans[i] << "\n";
    }
    return 0;   
}