#include<iostream>
#include<map>
#include<string>
using namespace std;

int main() {
    int n; cin >> n;
    map<string, int> cnt;
    for(int i = 0; i < n; ++i) {
        string s; cin >> s;
        if(cnt[s] == 0) printf("%s\n", s.c_str());
        else printf("%s(%d)\n", s.c_str(), cnt[s]);
        cnt[s]++;
    }
    return 0;
}