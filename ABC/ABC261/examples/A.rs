use proconio::input;
use std::cmp;

fn main() {
    input!{
        l1: i32, r1: i32,
        l2: i32, r2: i32,
    }
    let l = cmp::max(l1, l2);
    let r = cmp::min(r1, r2);
    println!("{}", cmp::max(0, r-l));
}
