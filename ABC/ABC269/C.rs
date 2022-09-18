use proconio::input;

fn main() {
    input!{
        n: i64
    }
    let mut ans: Vec<i64> = vec![0];
    for i in 0..60 {
        if (n>>i)&1 == 1 {
            let size = ans.len();
            for j in 0..size {
                ans.push(ans[j] + (1<<i));
            }
        }
    }
    ans.sort();
    for a in ans {
        println!("{}", a);
    }
}