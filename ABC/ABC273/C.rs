use proconio::input;
use std::collections::{HashMap, BTreeMap};

fn main() {
  input! {
    n: usize,
    arr: [usize; n],
  }
  
  let mut ans = vec![0; n];
  let mut map: BTreeMap<usize, usize> = BTreeMap::new();
  for val in arr {
    map.entry(val)
      .and_modify(|e| *e += 1)
      .or_insert(1);
  }
  let mut num = 0;
  for (k, v) in map.iter().rev() {
    ans[num] += v;
    num += 1;
  }
  for k in 0..n {
    println!("{}", ans[k]);
  }
}