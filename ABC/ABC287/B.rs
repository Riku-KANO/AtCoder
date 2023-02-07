use proconio::{input, fastout};
use num::Float;
use std::collections::HashSet;

fn main() -> Result<(), ()> {
  input! {
    n: usize, m: usize,
    s: [u32; n],
    t: [u32; m]
  }
  let mut ans = 0;
  let vs: Vec<u32> = s.iter().map(|&x|x%1000).collect();
  let memo = t.iter().collect::<HashSet<_>>();
  for ss in &vs {
    if memo.contains(ss) {ans += 1}
  }
  println!("{}", ans);
  Ok(())
}