use::proconio::input;

fn main() {
    let MOD: usize = 998244353;
    input! {
        N: usize,
        A: [usize; N],
    }
    let mut dp = vec![vec![vec![vec![0; N + 1]; N + 1]; N + 1]; N + 1]; // 先頭からi個からj個選んだ時にkで割った余りがlであるものの個数.
    for k in 0..N + 1 {
        dp[0][0][k][0] = 1;
    }
    for i in 0..N {
        for j in 0..i+1 {
            for k in 1..N+1 {
                for l in 0..k {
                    dp[i+1][j+1][k][(l + A[i])%k] += dp[i][j][k][l]; //選ぶ
                    dp[i+1][j][k][l] += dp[i][j][k][l]; // 選ばない
                    dp[i+1][j+1][k][(l+A[i])%k] %= MOD;
                    dp[i+1][j][k][l] %= MOD;
                }
            }
        }
    }
    let mut ans: usize = 0;
    for i in 1..N+1 {
        ans += dp[N][i][i][0];
        ans %= MOD;
    }
    println!("{}", ans);
    
}