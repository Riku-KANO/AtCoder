use proconio::{input, fastout};

fn main() {
    input! {
        N: usize, L: i32, K: i32,
        mut A: [i32; N]
    } 
    A.push(L);
    let mut l: i32 = 0;
    let mut r: i32 = 1e9 as i32;
    let mut m: i32;
    while l + 1!= r {
        m = (l + r) / 2;
        let mut from: i32 = 0;
        let mut cnt: i32 = 0;
        for i in 0..N+1 as usize{
            if A[i]-from >= m {
                cnt += 1;
                from = A[i];
            }
        }
        if cnt > K {
            l = m;
        } else {
            r = m;
        }
    }
    println!("{}", l);
}