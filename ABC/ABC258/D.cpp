#include<iostream>
using namespace std;
#define rep(i, n) for(int i = 0; i < n; ++i)
#define ll long long

const ll LINF = 1001002003004005006ll;

void solve() {
	int n; cin >> n;
	int x; cin >> x;
	ll cum[n+1]={};
	ll a[n], b[n];
	rep(i, n){
		cin >> a[i] >> b[i];
		cum[i+1]=cum[i]+a[i]+b[i];
	}
	ll ans = LINF;
	rep(i, n) {
		ll remain = x-i-1;
		if(remain<0)break;
		ans = min(ans, cum[i+1]+b[i]*remain);
	}
	cout << ans << endl;
}

int main(void){
  solve();
  return 0;
}
