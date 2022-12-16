use proconio::input;
use std::collections::HashSet;

fn solve1() {
    input! {
        n: usize, m: usize,
        arr: [i64; n]
    }
    let mut v: Vec<Vec<usize>> = vec![vec![]; m + 1];
    for i in 1..=n {
        let a = arr[i-1];
        let l = std::cmp::max((0 - a + i as i64 - 1) / i as i64, 0);
        let r = std::cmp::min(((n as i64 + 1 - a) / i as i64), m as i64);
        if r < 0 || l > m as i64{
            continue;
        }
        // eprintln!("l: {}, r: {}", l, r);
        for idx in l..=r {
            v[idx as usize].push((a + (i * idx as usize) as i64) as usize);
        }
    }
    for i in 1..=m {
        v[i].sort();
        // eprintln!("{:?}", v[i]);
        let mut ans = 0;
        for e in v[i].iter() {
            if *e == ans {
                ans = ans + 1;
            } else if *e + 1 == ans {
                continue;
            } else if *e > ans {
                break;
            }
        }
        println!("{}", ans);
    }
}

fn solve2() {
    input! {
        n: usize, m: usize,
        arr: [i32; n],
    }
    let mut v: Vec<HashSet<usize>> = vec![HashSet::new(); m];
    for i in 0..n {
        let a = arr[i];
        let i = i + 1;
        let mut s = 1;
        if a < 0 {
            s = (a.abs() as usize + i - 1) / i;
        }
        for k in s..=m {
            let x = (a + (i * k) as i32) as usize;
            if x > n {
                break;
            }
            v[k-1].insert(x);
        }
    }
    for vv in &v {
        for i in 0.. {
            if vv.contains(&i) {continue}
            println!("{}", i);
            break;
        }
    }
}

fn main() {
    solve2();
}