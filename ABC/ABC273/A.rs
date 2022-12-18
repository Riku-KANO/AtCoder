use proconio::input;

fn rec(x: u32) -> u32 {
  return if x == 0 {1} else {rec(x-1) * x};
}

fn main() {
  let mut buf = String::new();
  std::io::stdin().read_line(&mut buf).ok();
  buf = buf.trim().to_string();
  let n = buf.parse::<u32>().unwrap();
  println!("{}", rec(n));
}