#include<bits/stdc++.h>
using namespace std;

int main() {
    int R, C; cin >> R >> C;
    int A[2][2];
    cin >> A[0][0] >> A[0][1] >> A[1][0] >> A[1][1];
    cout << A[--R][--C] << endl;
    return 0;
}