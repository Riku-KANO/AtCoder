#include<iostream>
#include<string>

using namespace std;
#define rep(i, n) for(int i = 0; i < n; ++i)
#define ll long long

void solve() {
	int n, q; cin >> n >> q;
	string s; cin >> s;
	int id=0;
	while(q--) {
		int t, x; cin >> t >> x;
		if(t == 1) {
			id+=x;
			id%=n;
		} else if(t==2) {
			int wh=(n*10-id+x-1)%n;
			cout << s[wh] << endl;
		}
	}
}

int main(void){
  solve();
  return 0;
}
