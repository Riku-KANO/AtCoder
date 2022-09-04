NAME="0024"
FILE="$NAME.txt"
GROUP="sample"
INPUT="../in/$FILE"
OUTPUT="../output/$GROUP/out$FILE"
g++ solve2.cpp
./judge ./a.out < $INPUT 1> $OUTPUT
