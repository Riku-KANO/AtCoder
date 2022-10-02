NAME="0000"
NUM_GRID=31
GROUP="sample/N$NUM_GRID"
INPUT="../input/sample/N$NUM_GRID/in$NAME.txt"
OUTPUT="./out_trial.txt"
g++ -D_TEST=1 -o solve solve.cpp 
date > ./trial_info.txt;
date
./solve  < $INPUT 1> $OUTPUT 2>> ./trial_info.txt
date