import os
import re
import subprocess

EXE_MAIN = "./bin/solve"
REF = "out/sub16"
INPUT_DIR="./in/normal"

def run_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    return result.stdout, result.stderr

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for i in range(100):
    print(f"### CASE {i:03d} ######")

    output_directory = f"./out/simple/{i:04d}"
    create_directory_if_not_exists(output_directory)

    command = f"{EXE_MAIN} < {INPUT_DIR}/{i:04d}.txt 1> {output_directory}/out.txt 2> {output_directory}/err.txt"
    stdout, stderr = run_command(command)

    with open(f"{output_directory}/err.txt", "r") as f:
        s_target = f.read()

    with open(f"{REF}/{i:04d}/err.txt", "r") as f:
        s_ref = f.read()

    score_pattern = r".*Score\s*:\s*(\d+)"
    match_target = re.search(score_pattern, s_target)
    match_ref = re.search(score_pattern, s_ref)

    score_target = int(match_target.group(1)) if match_target else -1
    score_ref = int(match_ref.group(1)) if match_ref else -1

    print(f"score(main): {score_target:10d}")
    print(f"score(ref ): {score_ref:10d}")
    print(f"ratio: {score_ref / score_target if score_target != 0 else 'N/A'}\n")
