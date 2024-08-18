use clap::Parser;

use io::{read_input, read_props, Args, Input, Props};
use solve::solve;
use threadpool::ThreadPool;

mod io;
mod solve;

fn main() -> Result<(), ()> {
    let args: Args = Args::parse();
    let props: Props = read_props(&args);

    let result = match args.input_path {
        Some(input_path) => {
            let dir = std::path::Path::new(&input_path);
            if !dir.is_dir() {
                eprintln!("{} is not a directory", dir.display());
                return Err(());
            }

            let files: Vec<Input> = std::fs::read_dir(dir)
                .unwrap()
                .filter_map(|entry| entry.ok())
                .filter(|entry| entry.path().is_file())
                .map(|entry| read_input(entry.path().to_str()))
                .collect();

            let pool = ThreadPool::new(args.workers);
            
            let (tx, rx) = std::sync::mpsc::channel();
            files.into_iter().for_each(|input| {
                let tx = tx.clone();
                pool.execute(move || {
                    tx.send(solve(&input, &props)).expect("send failed");
                });
            });

            pool.join();

            let results: Vec<Result<io::Metric, String>> = rx.iter().collect();

            results 
        }
        None => {
            let input = read_input(None);
            let result = vec![solve(&input, &props)];

            result
        }
    };

    Ok(())
}
