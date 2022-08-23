use::proconio::input;

fn main() {
    input! {
        N: usize, M: usize,
    }
    let mut adj = vec![vec![false; N + 1]; N + 1];
    for _ in 0..M {
        input!{u:usize, v:usize}
        adj[u][v] = true;
        adj[v][u] = true;
    }
    let mut ans = 0;
    for i in 1..N+1 {
        for j in 1..i {
            for k in 1..j {
                if adj[i][j]&adj[j][k]&adj[k][i] {ans += 1};
            }
        }
    }
    println!("{}", ans);
}