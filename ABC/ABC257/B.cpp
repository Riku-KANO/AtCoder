#include<iostream>
using namespace std;

int main() {
    int n, k, q; cin >> n >> k >> q;
    int a[n + 1] = {};
    a[n]=1;
    for(int i = 0; i < k; ++i) {
        int tmp; cin >> tmp;
        a[--tmp]=1;
    }
    for(int i = 0; i < q; i++) {
        int l; cin >> l;
        int cnt = 0;
        int idx;
        for(int j = 0; j < n; ++j) {
            if(a[j])cnt++;
            if(cnt==l){
                idx=j;
                break;
            }
        }
        if(a[idx+1]!=1){
            a[idx]=0;
            a[idx+1]=1;
        }
    }
    // output
    for(int i = 0; i < n; ++i) {
        if(a[i]) cout << i + 1 << endl;
    }
    return 0;
}