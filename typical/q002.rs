use proconio::{input, fastout, derive_readable, marker::Chars};

fn rec(s: String, n: i32, l: i32, r: i32, v: &mut Vec<String>) ->Vec<String> {
    if n == 0 {
        v.push(s);
        return v.to_vec();
    } else {
        if l < n / 2 {
            v.extend(rec(s+"(", n, l + 1, r, &mut v).iter().cloned());
        }
        if r < n / 2 && r < l{
            v.extend(rec(s+")", n, l, r+1, &mut v).iter().cloned());
        }
        return v.to_vec();
    }
}

#[fastout]
fn main() {
    input! {
        N: i32
    } 
    if N % 2 == 1 {
        return;
    }
    let mut ans: Vec<String> = Vec::new();
    let s: String = "".to_string();
    rec(s , N, 0, 0, &mut ans);
    ans.sort();
    for s in ans {
        println!("{}", s);
    }
}