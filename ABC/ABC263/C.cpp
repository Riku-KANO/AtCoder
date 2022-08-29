#include<iostream>
#include<map>
#include<vector>

using namespace std;

int n, m;
map<vector<int>, bool> memo;

void rec(int s, int r, vector<int> v) {
	if(memo[v])return;
	memo[v] = true;
	if(r >= 1 && s > m) return;
	if(r == 0) {
		for(auto num: v) cout << num << " ";
		cout << endl;
		return;
	}
	for(int i = s; i <= m; ++i) {
		vector<int> tmp = v;
		tmp.push_back(i);
		rec(i + 1, r-1, tmp);
	}
}

void solve() {
	cin >> n >> m;
	vector<int> v;
	for(int i = 1; i <= m+1; ++i) {
		rec(i, n, v);
	}
}

int main(void){
  solve();
  return 0;
}
