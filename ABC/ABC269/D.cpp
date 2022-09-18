#include<iostream>
#include<map>
#include<atcoder/all>
using namespace std;
using namespace atcoder;
#define rep(i, n) for(int i = 0; i < (n); ++i)

const int di[6] = {-1, -1, 0, 0, 1, 1};
const int dj[6] = {-1, 0, -1, 1, 0, 1};
map<pair<int,int>, int> memo;
void solve() {
	int n; cin >> n;
	dsu uf(n + 1);
	for(int i = 0; i < n; i++) {
		int a, b; cin >> a >> b;
		memo[{a, b}] = i + 1;
		for(int j = 0; j < 6; j++) {
			int nx = a + di[j];
			int ny = b + dj[j];
			if(memo[{nx, ny}] > 0) {
				uf.merge(memo[{nx, ny}], i + 1);
			}
		}
	}
	cout << uf.groups().size() -1<< endl;
}

int main(void){
  solve();
  return 0;
}
