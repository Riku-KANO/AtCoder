#include<iostream>
using namespace std;

void solve() {
	int k; cin >> k;
	int d = k/60;
	int r = k%60;
	int h = 21 + d;
	printf("%2d:%02d\n", h, r);
}

int main(void){
  solve();
  return 0;
}
