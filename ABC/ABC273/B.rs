use proconio::input;

fn rec(x: u32) -> u32 {
  return if x == 0 {1} else {rec(x-1) * x};
}

fn main() {
  input! {
    mut n: u64, k: u32,
  }
  for i in 1..=k {
    let d = 10_u64.pow(i);
    let r = n % d;
    n /= d;
    n *= d;
    if r >= 5 * d / 10 { n += d }
  }
  println!("{}", n);
}