use proconio::{input, fastout, marker::Chars};
use std::{io, fmt::Debug, str::FromStr};
use std::io::{stdout, Write, stdin};
use std::collections::{HashMap, HashSet};

fn dfs(start: usize, pre: usize, goal: usize, graph: &Vec<Vec<usize>>, path: &mut Vec<i32>) {
    if path[start] != -1 {
        return;
    }
    if start == goal {
        path[start] = pre as i32;
        return;
    }
    path[start] = pre as i32;
    for nx in graph[start].iter() {
        dfs(*nx, start, goal, graph, path);
    }
}

#[fastout]
fn main() -> io::Result<()> {
    input! {
        n: usize, x: usize, y: usize
    }
    let mut G: Vec<Vec<usize>> = vec![vec![]; n + 1];
    for _i in 0..n-1 {
        input!{
            u: usize, v: usize
        }
        G[u].push(v);
        G[v].push(u);
    }
    let mut path = vec![-1; n + 1];
    dfs(x, x, y, &G, &mut path);
    let mut now = y;
    let mut ans = vec![];
    loop {
        if path[now] == now as i32{
            break;
        }
        ans.push(now);
        now = path[now] as usize;
    }
    ans.push(x);
    ans.reverse();
    for a in ans {
        println!("{}", a);
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