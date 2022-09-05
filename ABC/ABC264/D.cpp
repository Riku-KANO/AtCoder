#include<iostream>
#include<string>
using namespace std;
#define ll long long

void solve() {
	string target ="atcoder";
	int ans = 0;
	string s; cin >> s;
	for(int i = 0; i < target.size(); ++i) {
		char ch = target[i];
		int idx = 0;
		while(1) {
			if(s[idx]==ch) break;
			idx++;
		}
		ans+=idx;
		s.erase(s.begin()+idx);
	}
	cout << ans << endl;
}

int main(void){
  solve();
  return 0;
}
