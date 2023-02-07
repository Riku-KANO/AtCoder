#include<bits/stdc++.h>
#include<atcoder/lazysegtree>
using namespace std;
using namespace atcoder;
using ll = long long;

struct Query {
  int t;
  int x, y;
};

struct Card {
  ll a, b;
};

struct S {
  ll a;
  int size;
  S() {
    this->a = 0;
    this->size = 1;
  }
  S(ll _a, int _size): a(_a), size(_size){}
};

struct F {
  ll a, b;
};
S op(S l, S r) { return S{l.a + r.a, l.size + r.size}; }
S e() { return S{0LL, 0LL}; }
S mapping(F l, S r) { return S{r.a * l.a + r.size * l.b, r.size}; }
F composition(F f, F g) { return F{g.a * f.a, g.b * f.a + f.b}; }
F id() { return {1, 0}; }

int main() {
  int n;
  cin >> n;
  ll a[n];
  ll b[n];
  vector<Card> cards;
  vector<ll> vals;
  for(int i = 0; i < n; i++) {
    cin >> a[i] >> b[i];
    cards.push_back({a[i], b[i]});
    vals.push_back(a[i]);
  }
  int q;
  cin >> q;
  vector<Query> queries(q);
  for(int i = 0; i < q; i++) {
    int t; cin >> t;
    if(t == 1 or t == 2) {
      int x, y; cin >> x >> y;
      x--;
      queries[i] = {t, x, y};

      if(t == 1) {
        vals.push_back(y);
      }
    } else {
      int x; cin >> x;
      queries[i] = {t, x};
    }

  }
  sort(vals.begin(), vals.end());
  vector<ll> vc = vals;
  vc.erase(unique(vc.begin(), vc.end()), vc.end());
  map<int, int> memo;
  int r = vc.size();
  for(int i = 0; i < r; i++) memo[vc[i]] = i;

  vector<S> v(r), vsum(r);
  for(int i = 0; i < n; i++) {
    cards[i].a = memo[a[i]];
    int idx = memo[a[i]];
    v[idx].a += b[i];
    vsum[idx].a += a[i] * b[i];
  }
  lazy_segtree<S, op, e, F, mapping, composition, id> seg(v), seg_sum(vsum);
  for(int i = 0; i < q; i++) {
    Query query = queries[i];
    if(query.t == 1) {
      int x = query.x;
      int y = query.y;
      int pre = cards[x].a;
      int next = memo[y];
      ll num = cards[x].b;
      auto tmp = seg.get(pre);
      tmp.a -= num;
      seg.set(pre, tmp);
      tmp = seg.get(next);
      tmp.a += num;
      seg.set(next, tmp);
      tmp = seg_sum.get(pre);
      tmp.a -= num * vc[pre];
      seg_sum.set(pre, tmp);
      tmp = seg_sum.get(next);
      tmp.a += num * vc[next];
      seg_sum.set(next, tmp);
      cards[x].a = next;
    } else if(query.t == 2) {
      int x = query.x;
      ll next_num = query.y;
      ll pre_num = cards[x].b;
      int pos = cards[x].a;
      auto tmp = seg.get(pos);
      tmp.a += next_num - pre_num;
      seg.set(pos, tmp);
      tmp = seg_sum.get(pos);
      tmp.a += vc[pos] * (next_num - pre_num);
      seg_sum.set(pos, tmp);
      cards[x].b = next_num;
    } else if(query.t == 3) {
      ll x = query.x;
      auto g = [&x](const S s) {
        return x > s.a;
      };

      int l = seg.min_left(r, g);
      S s = seg.prod(l, r);
      ll num = s.a;
      S t = seg_sum.prod(l, r);
      ll ans = t.a;
      if(l == 0) {
        cout << -1 << endl;
      } else {
        ans += (x - num) * vc[l-1];
        cout << ans << endl;
      }
    }
    
  }
  return 0;
}