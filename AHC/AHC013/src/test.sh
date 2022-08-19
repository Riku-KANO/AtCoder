SRC_FILE="./solve.cpp"
SCORE_SRC_FILE="./score.cpp"
g++ $SRC_FILE
g++ -o score $SCORE_SRC_FILE
EXE_FILE="./a"
SCORE_FILE="./score"
OUTPUT_PATH="../output"
SCORE_RESULT="$OUTPUT_PATH/result.txt"

if [ ! -d $OUTPUT_PATH ]; then
mkdir $OUTPUT_PATH
fi

for i in 0..19; do
    printf "#TEST: %03d\n" $i
    INPUT_FILE="../input/input$(printf "%04d" $i)"
    OUTPUT_FILE="$OUTPUT_PATH/output$(printf "%04d" $i)"
    cat $INPUT_FILE | $EXE_FILE > $OUTPUT_FILE
    cat $OUTPUT_FILE | $SCORE_FILE >> $SCORE_RESULT
done

printf "########## TOTAL SCORE ##########\n"

cat $SCORE_RESULT | awk '{sum+=$2} END {printf "score: %10d\n\n" sum}'

printf "#################################\n"