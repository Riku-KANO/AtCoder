use proconio::{input, marker::Chars, fastout};

#[fastout]
fn main() -> Result<(), ()> {
  input! {
    n: usize,
    arr: [usize; n]
  }
  let mut par = vec![0; 2 * n + 2];
  for i in (1..=n) {
    par[2 * i] = par[arr[i-1]] + 1;
    par[2 * i + 1] = par[arr[i-1]] + 1;
  }
  for k in (1..=2 * n + 1) {
    println!("{}", par[k]);
  }
  Ok(())
}