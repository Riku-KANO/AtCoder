use proconio::input;
use std::ops;

#[derive(Debug)]
struct Point {
  x: i32,
  y: i32
}

fn outer(a: &Point, b: &Point, o: &Point) -> bool {
  let oa = Point {x: a.x - o.x, y: a.y - o.y,};
  let ob = Point {x: b.x - o.x, y: b.y - o.y,};
  if oa.x * ob.y - oa.y * ob.x <= 0 {return false;}
  else {return true;}
}

fn main() {
  let mut v: Vec<Point> = Vec::new();
  for i in 0..4 {
    input! {x: i32, y: i32}
    v.push(Point{x: x, y: y});
  }

  let mut ans: bool = true;
  ans &= outer(&v[2], &v[0], &v[1]);
  ans &= outer(&v[3], &v[1], &v[2]);
  ans &= outer(&v[0], &v[2], &v[3]);
  ans &= outer(&v[1], &v[3], &v[0]);
  if ans {
    println!("Yes");
  } else {
    println!("No");
  }

}