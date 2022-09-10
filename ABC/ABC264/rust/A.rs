use proconio::input;

fn main() {
  input! {
    l: usize, r: usize
  }
  let s: &str = "atcoder";
  println!("{}", &s[l-1..r]);
}