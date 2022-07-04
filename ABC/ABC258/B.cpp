#include<iostream>
#include<string>
using namespace std;
using ll=long long;
#define rep(i, n) for(int i = 0; i < n; ++i)

int dx[8]={-1, -1, 0, 1, 1, 1, 0, -1};
int dy[8]={0, 1, 1, 1, 0, -1, -1, -1};

void solve() {
	int n; cin >> n;
	string a[n];
	for(int i = 0; i < n; ++i) {
		cin >> a[i];
	}
	ll ans = -1;
	rep(i, n) {
		rep(j, n) {
			rep(l, 8) {
				string tmp="";
				rep(k, n) {
					int x=(n*10+j+dx[l]*k)%n;
					int y=(n*10+i+dy[l]*k)%n;
					tmp+=a[y][x];
				}
				ans=max(ans, stoll(tmp));
			}
		}
	}
	cout << ans << endl;
}

int main(void){
  solve();
  return 0;
}
