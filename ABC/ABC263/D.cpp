#include<iostream>
using namespace std;
using ll = long long;
#define rep(i, n) for(int i = 0; i < (n); i++)

void solve() {
	int n; cin >> n;
	ll l, r; cin >> l >> r;

	ll cum[n + 1] = {};
	ll rcum[n+1] = {};
	ll a[n];
	ll sum = 0;
	ll lm[n+1]={};
	ll rm[n+1]={};
	rep(i, n) {
		cin >> a[i];
		sum += a[i];
		cum[i+1] = cum[i]+a[i];
	}
	rep(i, n) rcum[i+1] = rcum[i]+a[n-1-i];
	for(ll i = 0; i <= n; ++i) {
		lm[i] = cum[i]-i*l;
		rm[i] = rcum[i]-i*r;
	}
	for(int i = 1; i <= n; ++i) {
		lm[i] = max(lm[i], lm[i-1]);
		rm[i] = max(rm[i], rm[i-1]);
	}
	ll merit = 0;
	ll ans = sum;
	for(int i = 0; i <= n; ++i) {
		if(merit < lm[i]+rm[n-i]) {
			merit = lm[i]+rm[n-i];
			ans = sum-merit;
		}
	}
	cout << ans << endl;
}

int main(void){
  solve();
  return 0;
}
