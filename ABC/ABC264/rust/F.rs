use proconio::{input, marker::Chars};
use std::cmp;

fn chtoi(ch: char) -> usize {
    let ret: i32 = ch as i32 - 48;
    return ret as usize;
}

fn main() {
    input! {
        H: usize, W: usize,
        R: [i64; H],
        C: [i64; W],
        A: [Chars; H],
    }
    const INF: i64 = 1<<60;
    let mut dp = vec![vec![vec![vec![INF as i64; W+1]; H+1]; 2]; 2];
    dp[0][0][0][0] = 0;
    dp[1][0][0][0] = R[0];
    dp[0][1][0][0] = C[0];
    dp[1][1][0][0] = R[0] + C[0];
    for i in 0..H {
        for j in 0..W {
            for r in 0..2 {
                for c in 0..2 {
                    if i > 0 {
                        for nr in 0..2 {
                            if nr == 1 && (nr + r + chtoi(A[i][j]) + chtoi(A[i-1][j])) % 2 == 0 {
                                dp[nr][c][i][j] = cmp::min(dp[nr][c][i][j], dp[r][c][i-1][j] + R[i]);
                            } else if nr == 0 && (r + chtoi(A[i][j]) + chtoi(A[i-1][j])) % 2 == 0 {
                                dp[nr][c][i][j] = cmp::min(dp[nr][c][i][j], dp[r][c][i-1][j]);
                            }
                        }
                    }
                    if j > 0 {
                        for nc in 0..2 {
                            if nc == 1 && (nc + c + chtoi(A[i][j]) + chtoi(A[i][j-1])) % 2 == 0 {
                                dp[r][nc][i][j] = cmp::min(dp[r][nc][i][j], dp[r][c][i][j-1] + C[j]);
                            } else if nc == 0 && (c + chtoi(A[i][j]) + chtoi(A[i][j-1])) % 2 == 0{
                                dp[r][nc][i][j] = cmp::min(dp[r][nc][i][j], dp[r][c][i][j-1]);
                            }
                        }
                    }
                }
            }
        }
    }
    let mut ans = INF;
    for r in 0..2 {
        for c in 0..2 {
            ans = cmp::min(ans, dp[r][c][H-1][W-1]);
        }
    }
    println!("{}", ans);
}