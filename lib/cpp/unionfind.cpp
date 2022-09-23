template<class T>
class UnionFind {
  public:
    T n;
    std::vector<T> par;
    UnionFind(){}
    UnionFind(T _n): n(_n) {
      par.resize(_n, -1);
    }

    void merge(T a, T b) {
      T a_leader = this->leader(a);
      T b_leader = this->leader(b);
      par[a_leader] = b_leader;
    }

    bool same(T a, T b) {
      return this->leader(a) == this->leader(b);
    }

    T leader(T a) {
      while(par[a] != -1) {
        a = par[a];
      }
      return a;
    }
  
  std::vector<std::vector<T>> groups() {
    std::vector<std::vector<T>> ret;
    std::map<T, std::vector<int>> memo;
    for(T i = 0; i < n; i++) {
      if(par[i] == -1) {
        ret.push_back(std::vector<T> i);
      } else {
        memo[par[i]].push_back(i);
      }
    }
    for(std::pair<T, std::vector<int>> p: memo) {
      ret.push_back(p.second);
    }
    return ret;
  }
};