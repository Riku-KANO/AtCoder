#include<iostream>
using namespace std;

#define rep(i, n) for(int i = 0; i < n; ++i)

void solve() {
	int n; cin >> n;
	char a[n][n];
	rep(i, n)rep(j, n)cin >> a[i][j];
	bool flag = false;
	rep(i, n) {
		for(int j = i + 1; j < n; ++j) {
			if(a[i][j]=='D' && a[j][i] != 'D') flag=true;
			if(a[i][j]=='W' && a[j][i]!='L') flag=true;
			if(a[i][j]=='L' && a[j][i]!='W') flag=true;			
		}
	}
	if(flag) cout << "incorrect" << endl;
	else cout << "correct" << endl;
}

int main(void){
  solve();
  return 0;
}
