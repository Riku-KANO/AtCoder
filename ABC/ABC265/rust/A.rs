use proconio::input;

fn main() {
  input! {
    x: i32, y: i32, n: i32
  }
  if x * 3 > y {
    let ans = n/3*y+(n%3)*x;
    println!("{}", ans);
  } else {
    println!("{}", x*n)
  }
}