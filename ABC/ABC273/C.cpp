#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;

int main(void){
  int n; cin >> n;
	vector<int> a(n), b(n);
	set<int> st;
  for(int i = 0; i < n; i++) {
		cin >> a[i];
		b[i] = a[i];
		st.insert(a[i]);
	}

	vector<int> a2;
	for(auto s: st) {
		a2.push_back(s);
	}
	vector<int> ans(n + 10, 0);
	sort(a2.begin(), a2.end());
	for(int i = 0; i < n; i++) {
		int d = a2.end()- upper_bound(a2.begin(), a2.end(), b[i]);
		ans[d]++;
	}
	for(int i = 0; i < n; i++) {
		cout << ans[i] << endl;
	}
  return 0;
}
