NAME="0000"
NUM_GRID=31
GROUP="sample/N$NUM_GRID"
INPUT="../input/sample/N$NUM_GRID/in$NAME.txt"
SPECIAL_INPUT="../input/special/all.txt"
OUTPUT="./out_debug.txt"
g++ -D_DEBUG=1 -o solve solve.cpp 
./solve  < $SPECIAL_INPUT 1> $OUTPUT 2> ./debug_info.txt