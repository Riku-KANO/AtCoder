#include<bits/stdc++.h>
using namespace std;

const int N=200005;
int a[N], b[N], pa[N], pb[N], lc[N], rc[N], c;

int wk(int l,int r){
	if(l > r) return 0;
	c++;
	int m=pb[a[c]],x=a[c];
	if(m<l||m>r){
        cout << -1 << endl;
        exit(0);
    }
	lc[x] = wk(l, m-1);
	rc[x] = wk(m+1, r);
	return x;
}

int main(){
	int n; cin>>n;
	for(int i=1; i<=n; i++) cin >> a[i], pa[a[i]]=i;
	for(int i=1; i<=n; i++) cin >> b[i], pb[b[i]]=i;
	if(wk(1, n) != 1){
		cout << -1 << endl;
		return 0;
	}
	for(int i=1; i<=n; i++) cout << lc[i] << " " << rc[i] <<endl;
	return 0;
}