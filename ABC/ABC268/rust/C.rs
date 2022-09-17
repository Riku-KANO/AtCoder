use proconio::input;

fn main() {
    input! {
        n: usize,
        p: [usize; n]
    }
    let mut cnt = vec![0; n];
    for i in 0..n {
        cnt[(p[i]+n-i)%n] += 1;
        cnt[(p[i]+n+1-i)%n] += 1;
        cnt[(p[i]+n-1-i)%n] += 1;
    }
    let mut ans: i32 = 0;
    for c in cnt {
        ans = std::cmp::max(ans, c);
    }
    println!("{}", ans);
}