import os
import datetime
import pathlib
import argparse
import re
import sys
import logging
import concurrent.futures

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--src_file", type=str)
    parser.add_argument("-e", "--exe_file", type=str)
    parser.add_argument("-o", "--out_path", type=str, default=None)
    parser.add_argument("-b", "--build", type=bool, default=False)
    parser.add_argument("-j", "--job_worker", type=int, default=1)
    parser.add_argument("-d", "--default", type=bool, default=False)
    parser.add_argument("--num_test", type=int, default=10)
    parser.add_argument("--input_path", type=str, default="./in/normal")
    args = parser.parse_args()
    return args

def build(src_file: str, exe_file: str) -> None:
    os.system(f"g++ ./src/{src_file} -o bin/{exe_file} -std=c++17")

def job(i: int, input_path: pathlib.Path, output_path: pathlib.Path) -> None:
    os.system(f"./bin/{exe_file} < {input_path}/{i:04d}.txt 1> {output_path}/{i:04d}/out.txt 2> {output_path}/{i:04d}/err.txt")

if __name__ == "__main__":
    args = get_arg()
    src_file = args.src_file
    exe_file = args.exe_file
    out_path = pathlib.Path(args.out_path)
    input_path = pathlib.Path(args.input_path)

    if not input_path.exists():
        logger.error("Input path does not exist.")
        sys.exit(1)

    logger.info(out_path)
    do_build = args.build
    job_worker = os.cpu_count() if args.job_worker == -1 else args.job_worker
    num_test = args.num_test
    is_default_job = args.default
    if do_build:
        build(src_file, exe_file)

    if out_path is None:
        t_delta = datetime.timedelta(hours=9)
        JST = datetime.timezone(t_delta, 'JST')
        now = datetime.datetime.now(JST)
        out_path = pathlib.Path("./out/{}".format(now.strftime('%Y%m%d_%H%M%S')))

    out_path.mkdir(parents=True, exist_ok=True)
    for i in range(num_test):
        path = out_path / f"{i:04d}"
        path.mkdir(parents=True, exist_ok=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=job_worker) as tpe:
        futures = [tpe.submit(job, i, input_path, out_path) for i in range(num_test)]
        concurrent.futures.wait(futures)

    total_score = 0
    for i in range(num_test):
        with open(f"{out_path}/{i:04d}/err.txt") as f:
            text = f.read()
            result = re.findall(".*Score\s*:\s*([0-9]*)\s*", text)
            score = int(result[0])
            total_score += score
            logger.info(f"[case {i:04d}]: score: {score:10d}")

    logger.info(f"total_score: {total_score}")
