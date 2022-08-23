use::proconio::input;

fn main() {
    input! {
        Y: i32
    }
    let mut y = Y;
    loop {
        if y % 4 == 2 {
            println!("{}", y);
            return;
        }
        y += 1;
    }
}