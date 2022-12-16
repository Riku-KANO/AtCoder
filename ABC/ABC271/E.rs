use proconio::input;

const INF: usize = 1e18 as usize;   

fn main() {
    input! {
        n: usize, m: usize, k: usize,
        roads: [[usize; 3]; m],
    }
    let mut dist = vec![INF; n + 1];
    dist[1] = 0;
    for i in 0..k {
        input! {
            mut e: usize
        }
        e -= 1;
        let from = roads[e][0];
        let to = roads[e][1];
        let len = roads[e][2];
        if dist[from] != INF {
            dist[to] = std::cmp::min(dist[to], dist[from]+len);
        }
    }
    let ans: i64 = if dist[n] == INF {-1} else {dist[n] as i64};
    println!("{}", ans);
}