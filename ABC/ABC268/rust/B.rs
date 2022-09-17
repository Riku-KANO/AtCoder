use proconio::input;

fn main() {
    input! {
        s: String, t: String
    }
    let hoge: Vec<_> = t.match_indices(&s).collect();
    if s.len() > t.len() || hoge.len() == 0 {
        println!("No");
        return;
    }
    if hoge[0].0 == 0 {
        println!("Yes");
    } else {
        println!("No");
    }
}