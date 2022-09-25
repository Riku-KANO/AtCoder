use proconio::{input, fastout, marker::Chars};
use std::{io, fmt::Debug, str::FromStr};
use std::io::{stdout, Write, stdin};
use std::collections::{HashMap, HashSet};


#[fastout]
fn main() -> io::Result<()> {
    input! {
        n: usize, k: usize,
        s: Chars,
    }
    let mut char_memo: Vec<Vec<i32>> = vec![vec![-1; 26]; n];
    for i in 0..n {
        for j in 0..26 {
            let ch_id: usize = s[n-1-i] as usize - 'a' as usize;
            if ch_id == j {
                char_memo[n-1-i][ch_id] = (n-1-i) as i32;
            } else {
                if i == 0 {
                    continue
                }
                char_memo[n-i-1][j] = char_memo[n-i][j];
            }
        }
    }
    let mut ans = String::new();
    let mut idx = 0;
    // println!("{:?}", char_memo);
    while ans.len() != k {
        // println!("{}, {}, {}", idx, ans, k);
        let sz: i32 = ans.len() as i32;
        for i in 0..26 {
            if char_memo[idx][i] == -1 || n as i32 - char_memo[idx][i] < k as i32 - sz {
                continue;
            }
            ans.push((i as u8 + 'a' as u8) as char);
            idx = (char_memo[idx][i] + 1) as usize;
            break;
        }
    }
    println!("{}", ans);
    Ok(())
}

#[allow(dead_code)]
struct BIT<T> {
    n: usize,
    bit: Vec<T>
}

// 1-indexed
#[allow(dead_code)]
impl BIT<i32> {
    fn new(n: usize) -> Self {
        let bit: Vec<i32> = vec![0; n + 1];
        return Self{n: n + 1, bit: bit};
    }

    fn init(&mut self, n: usize) -> () {
        self.n = n + 1;
        self.bit = vec![0; n + 1];
    }

    fn add(&mut self, pos: usize, x: i32) -> () {
        let mut idx: usize = pos;
        let mut idx2: i32 = pos as i32;
        loop {
            if idx >= self.n {
                break;
            }
            self.bit[idx] += x;
            idx2 += idx2 & (-idx2);
            idx = idx2 as usize;
        }
    }

    // return v[1] + v[2] + ... + v[pos]
    fn sum(&mut self, pos: usize) -> i32 {
        let mut sum = 0;
        let mut idx: usize = pos;
        let mut idx2: i32 = pos as i32;
        loop {
            if idx2 <= 0 {
                break;
            }
            sum += self.bit[idx];
            idx2 -= idx2 & (-idx2);
            idx = idx2 as usize;
        }
        return sum;
    }
}

#[allow(dead_code)]
impl BIT<i64> {
    fn new(n: usize) -> Self {
        let bit: Vec<i64> = vec![0; n + 1];
        return Self{n: n + 1, bit: bit};
    }

    fn init(&mut self, n: usize) -> () {
        self.n = n + 1;
        self.bit = vec![0; n + 1];
    }

    fn add(&mut self, pos: usize, x: i64) -> () {
        let mut idx: usize = pos;
        let mut idx2: i64 = pos as i64;
        loop {
            if idx >= self.n {
                break;
            }
            self.bit[idx] += x;
            idx2 += idx2 & (-idx2);
            idx = idx2 as usize;
        }
    }

    // return v[1] + v[2] + ... + v[pos]
    fn sum(&mut self, pos: usize) -> i64 {
        
        let mut sum = 0;
        let mut idx: usize = pos;
        let mut idx2: i64 = pos as i64;
        loop {
            if idx2 <= 0 {
                break;
            }
            sum += self.bit[idx];
            idx2 -= idx2 & (-idx2);
            idx = idx2 as usize;
        }
        return sum;
    }
}

#[allow(dead_code)]
struct UnionFind<T> {
    n: usize,
    par: Vec<T>
}

#[allow(dead_code)]
impl UnionFind<i32> {
    fn new(n: usize) -> Self {
        let par: Vec<i32> = vec![-1; n];
        return Self{n, par};
    }
    
    fn init(&mut self, n: usize) {
        self.n = n;
        self.par.resize(self.n, -1);
    }

    fn merge(&mut self, a: i32, b: i32) -> () {
        let a_parent: usize = self.leader(a);
        let b_parent: usize = self.leader(b);
        self.par[a_parent] = b_parent as i32;
    }

    fn leader(&mut self, a: i32) -> usize {
        if self.par[a as usize] == -1 {
            return a as usize;
        }
        self.par[a as usize] = self.leader(self.par[a as usize]) as i32;
        return self.par[a as usize] as usize;
    }

    fn same(&mut self, a: i32, b: i32) -> bool {
        return self.leader(a) == self.leader(b);
    }

    fn groups(&mut self) -> Vec<Vec<i32>> {
        todo!();
    }
}


#[allow(dead_code)]
fn read_tokens<F: FromStr>() -> io::Result<Vec<F>>
where
    <F as FromStr>::Err: Debug,
{
    let mut line = String::new();
    io::stdin().read_line(&mut line)?;
    let mut result: Vec<F> = Vec::new();
    for token in line.trim().split(" ") {
        let value: F = token
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("{:?}", e)))?;
        result.push(value);
    }
    Ok(result)
}