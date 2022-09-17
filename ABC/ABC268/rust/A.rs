use proconio::input;
use std::collections::HashSet;

fn main() {
    input! {
        nums: [i32; 5]
    }
    let s: HashSet<_> = nums.iter().collect();
    println!("{}", s.len());
}