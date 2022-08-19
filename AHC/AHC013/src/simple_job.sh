NAME="sample_K2_N15_s0"
FILE="$NAME.txt"
GROUP="sample"
INPUT="../input/$GROUP/$FILE"
OUTPUT="../output/$GROUP/$FILE"
cat $INPUT | ./solve_copy > $OUTPUT