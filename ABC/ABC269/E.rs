use proconio::{input};
use std::{io, fmt::Debug, str::FromStr};
use std::io::{stdout, Write, stdin};

fn main() -> io::Result<()> {
    let n = read_tokens::<i32>().map(|v|v[0])?;
    let mut ceil: i32 = n;
    let mut bottom: i32= 1;
    let mut middle: i32;
    while bottom != ceil {
        middle = (ceil + bottom) / 2;
        // let mut str = String::new();
        println!("? {} {} {} {}", bottom, middle, 1, n);
        // stdout().flush().unwrap();
        let t = read_tokens::<i32>().map(|v| v[0])?;
        // stdin().read_line(&mut str);
        // let t: i32 = str.trim().parse().unwrap();
        if t == middle - bottom + 1 {
            bottom = middle + 1;
        } else {
            ceil = middle;
        }
    }

    let mut left: i32 = 1;
    let mut right: i32 = n;
    while left != right {
        middle = (left + right) / 2;
        // let mut str = String::new();
        println!("? {} {} {} {}", 1, n, left, middle);
        // stdout().flush().unwrap();
        let t: i32 = read_tokens::<i32>().map(|v| v[0])?;
        // stdin().read_line(&mut str);
        // let t: i32 = str.trim().parse().unwrap();
        if t == middle - left + 1 {
            left = middle + 1;
        } else {
            right = middle;
        }
    }
    println!("!{} {}", bottom, left);

    Ok(())
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