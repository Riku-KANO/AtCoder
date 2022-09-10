use proconio::input;
use std::cmp;

fn main() {
    input! {
        r: i32, c: i32
    }
    if cmp::max((8-r).abs(), (8-c).abs()) % 2 == 0 {
        println!("white");
    } else {
        println!("black");
    }
}