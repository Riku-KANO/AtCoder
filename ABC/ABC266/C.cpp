#include<iostream>
using namespace std;

struct Point {
  int x, y;
};

bool f(Point&a, Point&b, Point&c, Point&d) {
  int dx = a.x-b.x;
  int dy = a.y-b.y;
  int vc = dx * (c.y - a.y) - dy * (c.x - a.x);
  int vd = dx * (d.y - a.y) - dy * (d.x - a.x);
  if(vc == 0 || vd == 0) return false;
  if(vc * vd > 0) return true;
  return false;
}

int main() {
  Point a, b, c, d;
  cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;
  bool ans = true;
  ans &= f(a, b, c, d);
  ans &= f(b, c, d, a);
  ans &= f(c, d, a, b);
  ans &= f(d, a, b, c);
  if(ans) cout << "Yes" << endl;
  else cout << "No" << endl;
  return 0;
}