NAME="0000"
NUM_GRID=31
GROUP="sample/N$NUM_GRID"
INPUT="../input/sample/N$NUM_GRID/in$NAME.txt"
OUTPUT="../output/$GROUP/out$NAME.txt"
g++ solve.cpp -o solve
./solve < $INPUT 1> $OUTPUT 2> ./simple_job_debug.txt