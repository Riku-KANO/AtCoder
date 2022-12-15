use proconio::input;

// 遅いアルゴリズム
#[warn(dead_code)]
fn dfs(d: usize, sum: usize, n: usize, S: usize, str: &str, coin: &Vec<Vec<usize>>) -> String {
    // eprintln!("{}, {}, {}", d, sum, str);
    let mut ret = "".to_string();
    if d == n {
        if sum == S {
            eprintln!("{}, {}", str, sum);
            return str.to_string();
        } else {
            return "".to_string();
        }
    }
    let str1: &str = &(str.to_string() + "H");
    let str2: &str = &(str.to_string() + "T");
    let ret1 = dfs(d + 1, sum + coin[d][0], n, S, str1, &coin);
    let ret2 = dfs(d + 1, sum + coin[d][1], n, S, str2, &coin);
    if ret1 != "" {
        ret = ret1;
    }
    if ret2 != "" {
        ret = ret2;
    }
    return ret.to_string();
}


fn main() {
    let mut dp: Vec<Vec<usize>> = vec![vec![0; 10605]; 105];
    dp[0][0] = 1;
    input! {
        n: usize, s: usize,
        coin:[[usize; 2]; n]
    }
    for i in 0..n {
        for sum in 0..s {
            if dp[i][sum] == 0 {
                continue;
            }
            dp[i+1][sum + coin[i][0]] = 1;
            dp[i+1][sum + coin[i][1]] = 1;
        }
    }
    if dp[n][s] == 0 {
        println!("No");
    } else {
        let mut ans: String = "".to_string();
        let mut total = s;

        for j in (0..n).rev() {
            if total >= coin[j][0] {
                if dp[j][total - coin[j][0]] == 1 {
                    total -= coin[j][0];
                    let mut tmp = "H".to_string();
                    tmp.push_str(&ans);
                    ans = tmp;
                    continue;
                }
            }
            if total >= coin[j][1] {
                if dp[j][total - coin[j][1]] == 1 {
                    total -= coin[j][1];
                    let mut tmp = "T".to_string();
                    tmp.push_str(&ans);
                    ans = tmp;
                }
            }
        }
        println!("Yes\n{}", ans);
    }
}