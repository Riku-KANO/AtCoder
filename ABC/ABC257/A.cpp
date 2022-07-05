#include<iostream>
using namespace std;
int main(){
    int n, x; cin >> n >> x;
    char ans = 'A' + (x-1)/n ;
    cout << ans << endl;
    return 0;
}