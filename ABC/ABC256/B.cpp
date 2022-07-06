#include<iostream>
using namespace std;

int main() {
    int n; cin >> n;
    long long s=0;
    for(int i = 0; i < n; i++) {
        int a; cin >> a;
        s++;
        s<<=a;
        s%=(1<<4);
    }
    cout << n - __builtin_popcount(s) << endl;
    return 0;
}