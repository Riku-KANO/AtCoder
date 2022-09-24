use proconio::{input, fastout, marker::Chars};
use std::{io, fmt::Debug, str::FromStr};
use std::io::{stdout, Write, stdin};
use std::collections::{HashMap, HashSet};

fn rec(n: usize, r: usize, mut s: String, v: &mut Vec<String>){
    // println!("{} {} {} {}", n, r, s, s.len());
    if s.len() == n {
        v.push(s);
    } else {
        let sz = s.len();
        if r == n / 2{
            let sa = &mut s;
            sa.push(')');
            rec(n, r, sa.to_string(), v);
        } else if r == s.len() / 2 {
            let sa = &mut s;
            sa.push('(');
            rec(n, r + 1, sa.to_string(), v);
        } else {
            let sa = &mut s;
            (*sa).push('(');
            rec(n, r + 1, sa.to_string(), v);
            sa.pop();
            let sb = &mut s;
            (*sb).push(')');
            rec(n, r, sb.to_string(), v);
        }
    }
}

#[fastout]
fn main() -> io::Result<()> {
    input!{
        n: usize
    }
    if n % 2 == 1 {
        Ok(())
    } else {
        let mut v = Vec::new();
        rec(n, 0, "".to_string(), &mut v);
        v.sort();
        for ans in v {
            println!("{}", ans);
        }
        Ok(())
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