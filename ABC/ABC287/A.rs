use proconio::input;
use proconio::marker::Chars;
fn main() {
    input! {
        n: usize, m: usize,
        s: [String; n],
        t: [String; m]
    }
    let c = s.iter().filter(|&x|x=="For").count();
    println!("{}", if c > n / 2 {"Yes"} else {"No"});
}
