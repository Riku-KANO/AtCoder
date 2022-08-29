#include<iostream>
#include<map>
using namespace std;

void solve() {
	map<int, int> cnt;
	for(int i = 0; i < 5; ++i) {
		int a; cin >> a;
		cnt[a]++;
	}
	int b = 0, c = 0;
	for(auto p: cnt) {
		if(p.second == 2) b = 1;
		if(p.second == 3) c = 1;
	}
	if(b && c) {
		cout << "Yes" << endl;
	} else {
		cout << "No" << endl;
	}
}

int main(void){
  solve();
  return 0;
}
