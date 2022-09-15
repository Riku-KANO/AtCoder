use proconio::{input, fastout, derive_readable};
use std::collections::HashMap;
use std::hash::Hash;


#[derive_readable]
#[derive(PartialEq, Eq, Hash, Debug)]
struct Pair {
    a: i32,
    b: i32,
}

// impl FromStr for Pair {
//     type Err = ParseIntError;
//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         let pair: Vec<&str> = s.trim().split(' ').collect();
//         let a_from: i32 = pair[0].parse::<i32>()?;
//         let b_from: i32 = pair[1].parse::<i32>()?;
//         println!("{} {}",a_from, b_from);
//         Ok(Pair {a: a_from, b: b_from})
//     }
// }

#[fastout]
fn main() {
    let dx: Vec<i32> = vec![-1, 0, 1, -1, 1, -1, 0, 1, 0];
    let dy: Vec<i32> = vec![1, 1, 1, 0, 0, -1, -1, -1, 0];
    input! {
        h: i64, w: i64, n: usize,
        p: [Pair; n]
    }
    let mut cnt: HashMap<Pair, i64> = HashMap::new();
    for i in 0..n {
        for j in 0..9 {
            let nx = p[i as usize].b + dx[j as usize] - 1;
            let ny: i32 = p[i as usize].a + dy[j as usize] - 1;
            if nx >= 1 && nx <= w as i32 - 2 && ny >= 1 && ny <= h as i32 - 2 {
                if cnt.contains_key(&Pair{a: nx, b: ny}) {
                    cnt.entry(Pair{a: nx, b: ny}).and_modify(|v: &mut i64|{*v+=1});
                } else {
                    cnt.insert(Pair{a: nx, b: ny}, 1);
                }
            }
        }
    }
    let mut ans:Vec<i64> = vec![0;10];
    let total: i64 = (h-2)*(w-2);
    let mut zero: i64 = total;
    for (k, v) in cnt {
        println!("k: ({},{}), v: {}", k.a, k.b, v);
        ans[v as usize] += 1;
        zero = zero - 1;
    }
    ans[0] = zero;
    for i in 0..10 {
        println!("{}", ans[i as usize]);
    }
}