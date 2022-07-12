#include<iostream>
#include<string>
#include<vector>
using namespace std;

vector<pair<char, int>> rle(string s) {
    vector<pair<char, int>> ret;
    int n = s.size();
    char cur;
    for(int i = 0; i < n; i++) {
        if(i == 0 || cur != s[i]) {
            cur = s[i];
            ret.push_back({s[i], 1});
        } else if(cur == s[i]) {
            ret.back().second++;
        }
    }
    return ret;
}

int main() {
    string s, t; cin >> s >> t;
    vector<pair<char, int>> vs, vt;
    vs = rle(s); vt = rle(t);
    if(vs.size() != vt.size()) {
        puts("No");
        exit(0);
    }
    int n = vs.size();
    bool judge = false;
    for(int i = 0; i < n; i++) {
        if(vs[i].first != vt[i].first) judge=true;
        else if(vs[i].first == vt[i].first) {
            if(vs[i].second > vt[i].second) judge = true;
            if(vs[i].second == 1 && vt[i].second > 1) judge = true;
        }
    }
    if(judge) puts("No");
    else puts("Yes");
    return 0;
}