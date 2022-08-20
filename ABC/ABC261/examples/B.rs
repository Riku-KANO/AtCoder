use proconio::{input, marker::Chars};

fn main() {
    input! {
        n: usize,
        s: [Chars; n],
    }
    let mut ans: bool = true;
    for i in 0..n {
        for j in i+1..n {
            if s[i][j] == 'W' && s[j][i] != 'L' {ans = false};
            if s[i][j] == 'L' && s[j][j] != 'W' {ans = false};
            if s[i][j] == 'D' && s[j][i] != 'D' {ans = false};
        }
    }
    if ans {println!("correct");}
    else {println!("incorrect");}
}