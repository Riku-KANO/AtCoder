use proconio::input;

fn main() {
    input! {
        n: usize,
        arr: [usize; n]
    }
    let mut evens: Vec<usize> = arr.iter().cloned().filter(|&x| x%2==0).collect();
    let mut odds: Vec<usize> = arr.iter().cloned().filter(|&x| x%2==1).collect();
    evens.sort();
    odds.sort();
    let mut ans: i64 = -1;
    let ne = evens.len();
    let no = odds.len();
    if ne >= 2 {
        ans = std::cmp::max(ans, (evens[ne-1] + evens[ne-2]) as i64);
    }
    if no >= 2 {
        ans = std::cmp::max(ans, (odds[no-1] + odds[no-2]) as i64);
    }
    println!("{}", ans);
}