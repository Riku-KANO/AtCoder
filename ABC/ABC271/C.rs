use proconio::input;
use std::collections::{VecDeque, HashSet};

fn main() {
    input! {
        n: usize,
        mut books: [usize; n]
    }
    books.sort();
    let mut set: HashSet<usize> = HashSet::new();
    for b in books {
        set.insert(b);
    }
    let n_doubled = n-set.len();
    let mut books2 = Vec::new();
    for b in set {
        books2.push(b);
    }
    books2.sort();
    let mut dq: VecDeque<usize> = VecDeque::from(books2);
    let mut ans: usize = 0;
    for i in 0..n_doubled {
        dq.push_back(1e9 as usize);
    }
    loop {
        let front = dq.pop_front().unwrap();

        if front == ans + 1 {
            ans = ans + 1;
        } else if front == ans {
            dq.push_back(1e9 as usize + 5);
        } else {
            dq.push_front(front);
            let _hoge = dq.pop_back();
            if dq.is_empty() {
                break;
            }
            let _foo = dq.pop_back();
            ans = ans + 1;   
        }
        if dq.is_empty() {
            break;
        }
     }
    println!("{}", ans);
}