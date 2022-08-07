#include<iostream>
#include<string>
#include<map>
using namespace std;

int main() {
    string s; cin >> s;
    map<char, int> cnt;
    for(int i = 0; i < (int)s.size(); ++i ) cnt[s[i]]++;
    for(auto p: cnt) {
        if(p.second == 1) {
            cout << p.first << endl;
            return 0;
        } 
    }
    cout << -1 << endl;
    return 0;
}