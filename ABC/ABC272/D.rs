use proconio::input;
use std::collections::VecDeque;

fn in_grid(i: i64, j: i64, n: usize) -> bool {
    return i >= 0 && j >= 0 && i < n as i64 && j < n as i64;
}

fn main() {
    input! {
        n: usize, m: i64
    }
    let mut dq = VecDeque::new();
    let mut nx :Vec<(i64,i64)> = Vec::new();
    for i in -400..=400 {
        for j in -400..=400 {
            if i * i + j * j == m {
                nx.push((i,j));
            }
        }
    }
    let mut dist: Vec<Vec<i64>> = vec![vec![-1; n]; n];
    dist[0][0] = 0;
    dq.push_back((0,0));
    while ! dq.is_empty() {
        let (ci, cj) = dq.pop_front().unwrap();
        for (di,dj) in &nx {
            let ni = ci + di;
            let nj = cj + dj;
            if in_grid(ni, nj, n){
                if dist[ni as usize][nj as usize] == -1 {
                    dist[ni as usize][nj as usize] = dist[ci as usize][cj as usize] + 1;
                    dq.push_back((ni, nj));
                }
            }
        }
    }
    for d in dist {
        for e in d {
            print!("{} ", e);
        }
        println!("");
    }
}