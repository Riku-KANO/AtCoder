use clap::Parser;
use std::io::BufRead;
use proconio::input;

#[derive(Parser)]
pub struct Args {
    pub input_path: Option<String>,

    #[arg(default_value = "2.0")]
    time_limit: f32,

    #[arg(short, long, default_value = "1")]
    pub workers: usize,
}

#[derive(Debug, Copy, Clone)]
pub struct Props {
    time_limit: f32,
}

impl Default for Props {
    fn default() -> Self {
        Self { time_limit: 2.0 }
    }
}

pub fn read_props(args: &Args) -> Props {
    Props {
        time_limit: args.time_limit,
    }
}

pub struct Input {
    x: i32,
    y: i32,
}

pub fn read_input(file_path: Option<&str>) -> Input {
    match file_path {
        Some(file_path) => {
            let file = std::fs::File::open(file_path).unwrap();
            let reader = std::io::BufReader::new(file);

            let mut lines = reader.lines();

            Input { x: 2, y: 3 }
        }
        None => {
            input! {
                x: i32,
                y: i32,
            }

            Input { x, y }
        }
    }
}

pub struct Metric {
    elapsed: f32,
    score: u64,
}

