use proconio::input;
use std::cmp::max;

const MIN: i64 = std::i64::MIN;

fn main() {
    input!{
        N: usize,
    }
    let mut snuke = vec![vec![0; 5]; 100_005];
    let mut Tmax = 0;
    for _ in 0..N {
        input! {
          T: usize, X: usize, A:i64
        }
        snuke[T][X] = A;
        Tmax = max(Tmax, T);
    }

    let mut dp = vec![vec![MIN; 5]; 100_005];
    dp[0][0] = 0;
    for t in 0..Tmax+1 {
        for i in 0..5 {
            if dp[t][i] == MIN {continue;}
            dp[t][i] += snuke[t][i];
            dp[t+1][i] = max(dp[t+1][i], dp[t][i]);
            if i > 0 {
                dp[t+1][i-1] = max(dp[t+1][i-1], dp[t][i]);
            }
            if i < 4 {
                dp[t+1][i+1] = max(dp[t+1][i+1], dp[t][i]);
            }
        }
    }

    let mut ans = 0;
    for i in 0..5 {
        ans = max(ans, dp[Tmax][i]);
    }
    println!("{}", ans);
  
   
}
