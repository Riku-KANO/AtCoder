#include<iostream>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
using namespace std;

class Solver {
  public:
    map<string, int> memo;

    void solve2() {
      int n, m; cin >> n >> m;
      int remain = 16-(n-1);
      vector<string> vs(n);
      vector<string> vt(m);
      for(int i = 0; i < n; i++) {
        cin >> vs[i];
        remain -= (int)vs[i].size();
      }
      for(int i = 0; i < m; i++) {
        cin >> vt[i];
        this->memo[vt[i]] = 1;
      }
      if(n == 1 and vs.size() < 3) {
        cout << -1 << endl;
        return ;
      }
      if(!dfs2(vs, remain)) {
        cout << -1 << endl;
      }
    }

    bool dfs2(vector<string> v, int r) {
      do {
        std::string now = "";
        if(r < 0) return false;
        for(int i = 0; i < (int)v.size(); i++){
          if(i == 0) now += v[i];
          else now += "_" + v[i]; 
        }
        if(now.size() >= 3 and now.size() <= 16 and this->memo[now] != 1 and now[0] != '_') {
          cout << now << endl;
          return true;
        }
      } while (next_permutation(v.begin(), v.end()));

      std::vector<std::string> nx = v;
      for(int i = 0; i < (int)v.size(); i++) {
        string buf = nx[i];
        nx[i] = "_" + nx[i];
        if(dfs2(nx, r-1)) return true;
        nx[i] = buf;
      }
      return false;
    }

    void solve() {
      int n, m; cin >> n >> m;
      int sz = n-1;
      vector<string> vs(n);
      vector<string> vt(m);
      for(int i = 0; i < n; i++) {
        cin >> vs[i];
        sz += vs[i].size();
      }
      for(int i = 0; i < m; i++) {
        cin >> vt[i];
        this->memo[vt[i]] = 1;
      }
      if(n == 1 and vs[0].size() < 3) {
        cout << -1 << "\n";
        return;
      }
      sort(vs.begin(), vs.end());
      do {
        vector<string> nx = vs;
        string now = "";
        for(int i = 0; i < (int)nx.size(); i++) {
          if(i > 0) nx[i] = '_' + vs[i];
          now += nx[i];
        }

        if(this->dfs(nx, 16-sz)) return;
      } while(next_permutation(vs.begin(), vs.end()));
      cout << -1 << endl;
    }

    bool dfs(vector<string> v, int r) {
      string now = "";
      for(string s: v) now += s;
      if(this->memo[now] != 1 and now.size() >= 3 and now.size() <= 16) {
        cout << now << "\n";
        return true;
      }
      if(r == 0) return false;
      bool flag = false;
      vector<string> vnx = v;
      for(int i = 1; i < (int)v.size(); i++) {
        vnx[i] = "_" + v[i];
        string buf = v[i];
        flag = dfs(vnx, r-1);
        vnx[i] = buf;
        if(flag) return true;
      }
      return false;
    }
};

int main () {
  Solver solver;
  solver.solve();
  return 0;
}