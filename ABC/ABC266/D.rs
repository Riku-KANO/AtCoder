use proconio::input;
use std::cmp::max;

const MIN: i64 = std::i64::MIN;

fn main() {
    input!{
        n: usize,
        txa: [(usize, usize, i64); n],
    }
    let mut x = vec![0; 100_001];
    let mut a = vec![0; 100_001];
    for (t, xi, ai) in txa {
        x[t] = xi;
        a[t] = ai;
    }

    let mut dp = vec![vec![MIN; 100_001]; 5];
    dp[0][0] = 0;
    for t in 1..=100_000 {
        for i in 0..5 {
            dp[i][t] = dp[i][t-1];
            if i > 0 {
                dp[i][t] = max(dp[i][t], dp[i-1][t-1]);
            }
            if i < 4 {
                dp[i][t] = max(dp[i][t], dp[i+1][t-1]);
            }
        }
        dp[x[t]][t] += a[t];
    }

    let mut ans = 0;
    for i in 0..5 {
        ans = max(ans, dp[i][100_000]);
    }
    println!("{}", ans);
}