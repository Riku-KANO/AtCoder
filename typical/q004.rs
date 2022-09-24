use proconio::{input, fastout, marker::Chars};
use std::{io, fmt::Debug, str::FromStr};
use std::io::{stdout, Write, stdin};
use std::collections::{HashMap, HashSet};

#[fastout]
fn main() -> io::Result<()> {
    input! {
        h: usize, w: usize,
        a: [[i32; w]; h]  
    }
    let mut sumh = vec![0; h];
    let mut sumw = vec![0; w];
    for i in 0..h {
        let mut s = 0;
        for j in 0..w {
            s += a[i][j];
        }
        sumh[i] = s;
    }
    for i in 0..w {
        let mut s = 0;
        for j in 0..h {
            s += a[j][i];
        }
        sumw[i] = s;
    }
    for i in 0..h {
        for j in 0..w {
            print!("{} ", sumh[i] + sumw[j] - a[i][j]);
        }
        println!();
    }
    Ok(())
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