use proconio::input;

fn main() -> Result<(), ()> {
  let mut buf = String::new();
  std::io::stdin().read_line(&mut buf).ok();
  let line = buf.split_whitespace().map(|x|x.parse().unwrap()).collect::<Vec<i32>>();
  let a = line[0];
  let b = line[1];
  println!("{:.3}", b as f32 /a as f32);
  Ok(())
}