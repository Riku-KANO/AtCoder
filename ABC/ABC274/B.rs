use proconio::{input, marker::Chars};

fn main() -> Result<(), ()> {
  input! {
    h: usize, w: usize,
    grid: [String; h]
  }
  let grid: Vec<Vec<char>> = grid
        .iter()
        .map(|row| row.chars().collect::<Vec<char>>())
        .collect();
  let grid_t: Vec<String> = (0..w).map(|j|{
    let str: String = (0..h).map(|i|{
      grid[i][j]
    }).collect();
    str
  }).collect();
  let mut ans = vec![0; w];
  for i in (0..w) {
    ans[i] = grid_t[i].matches("#").count();
  }
  println!("{}", ans.iter().map(|s| s.to_string()).collect::<Vec<String>>().join(" "));
  Ok(())
}