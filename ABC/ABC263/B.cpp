#include<iostream>
using namespace std;

void solve() {
	int n; cin >> n;
	int a[n + 1] = {};
	a[1] = 1;
	for(int i = 2; i <= n; ++i) {
		int b; cin >> b;
		a[i]=b;
	}
	int ans = 0;
	int now = n;
	while(1) {
		ans++;
		if(a[now]==1) {
			cout << ans << endl;
			return;
		} else now = a[now];
	}
}

int main(void){
  solve();
  return 0;
}
