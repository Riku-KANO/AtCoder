use proconio::input;

fn main() {
  input! {s: String}
  let mut ans: i32;
  if s == "Monday" {
    ans = 5;
  } else if s == "Tuesday" {
    ans = 4;
  } else if s == "Wednesday" {
    ans = 3;
  } else if s == "Thursday" {
    ans = 2;
  } else {
    ans = 1;
  }
  println("{}", ans);
}