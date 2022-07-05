#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<cmath>
using namespace std;

struct Person {
    int w;
    int label;
};

int main() {
    int n; cin >> n;
    string s; cin >> s;
    vector<Person> p(n);
    int correct = 0;
    for(int i = 0; i < n; ++i) {
        cin >> p[i].w;
        p[i].label=s[i]-'0';
        if(p[i].label)correct++;
    }
    auto compare = [&](const Person &a, const Person &b) {
        return a.w < b.w;
    };
    sort(p.begin(), p.end(), compare);
    int ans = max(correct, n-correct);
    for(int i = 0; i < n; ++i) {
        int base = p[i].w;
        for(int j = i; j < n; j++) {
            if(p[j].w != base) {
                i = j-1;
                break;
            } else {
                if(p[j].label) correct--;
                else correct++;
            }
            if(j == n-1)i=n-1;
        }
        ans = max(ans, correct);
    }
    cout << ans << endl;
    return 0;
}