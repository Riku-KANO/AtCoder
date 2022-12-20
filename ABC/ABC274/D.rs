use proconio::{input, marker::Chars, fastout};
use std::collections::HashSet;

const OFFSET: usize = 10020 as usize;

#[fastout]
fn main() -> Result<(), ()> {
  input! {
    n:usize,
    x:i32,
    y:i32,
    a:[i32;n]
  }

  let mut odds = vec![];
  let mut evens = vec![];

  for i in 1..n {
    if i % 2 == 0 {
      evens.push(a[i]);
    } else {
      odds.push(a[i]);
    }
  }

  let mut set = HashSet::new();
  set.insert(a[0]);
  for v in evens {
    let mut new_set = HashSet::new();
    for cv in set {
      new_set.insert(cv+v);
      new_set.insert(cv-v);
    }
    set = new_set;
  }

  if !set.contains(&x) {
    println!("No");
    return Ok(())
  }

  let mut set = HashSet::new();
  set.insert(0);
  for v in odds {
    let mut new_set = HashSet::new();
    for cv in set {
      new_set.insert(cv+v);
      new_set.insert(cv-v);
    }
    set = new_set;
  }

  if !set.contains(&y) {
    println!("No");
    return Ok(())
  }

  println!("Yes");
  Ok(())
}