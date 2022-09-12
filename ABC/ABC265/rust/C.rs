use proconio::{input, marker::Chars};

fn main() {
    input!{
        h: usize,
        w: usize,
        G: [Chars; h]
    }
    let mut visited = vec![vec![false; w]; h];
    let mut ci:usize= 0;
    let mut cj:usize= 0;
    loop {
        visited[ci][cj] = true;
        let mut ni: usize = ci;
        let mut nj: usize = cj;
        match G[ci][cj] {
            'L' => nj = cj - 1,
            'R' => nj = cj + 1,
            'U' => ni = ci - 1,
            'D' => ni = ci + 1,
            _ => (),
        }
        if (ni as i32) < 0 || (nj as i32) < 0 || ni == h || nj == w {
            println!("{} {}", ci + 1, cj + 1);
            return;
        } else if visited[ni][nj] {
            println!("{}", -1);
            return;
        } else {
            ci = ni;
            cj = nj;
        }
    }
}