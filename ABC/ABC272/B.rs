use proconio::input;

fn main() {
    input! {
        n: usize, m: usize
    }
    let mut adj = vec![vec![0; 101]; 101];
    for i in 0..m {
        input!{
            k: usize,
            x: [usize; k]
        }
        for i in x.clone() {
            for j in x.clone() {
                if i == j {
                    continue;
                }
                adj[i][j] = 1;
                adj[j][i] = 1;
            }
        }
    }
    let mut judge = true;
    for i in 1..=n {
        for j in i+1..=n {
            if adj[i][j] == 0 {
                judge = false;
            }
        }
    }
    if judge {
        println!("Yes");
    } else {
        println!("No");
    }
}