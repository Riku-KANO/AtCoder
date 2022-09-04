SRC_FILE="./solve2.cpp"
SCORE_SRC_FILE="./judge.cpp"
g++ $SRC_FILE
g++ --std=c++2a -o judge Problem.o $SCORE_SRC_FILE
EXE_FILE="./a.out"
SCORE_FILE="./judge"
OUTPUT_PATH="../output/"
SCORE_RESULT="$OUTPUT_PATH/result.txt"

if [ "$OUTPUT_PATH" = "../output/" ]; then
NOW=$(date "+%Y%m%d_%H%M%S")
OUTPUT_PATH="$OUTPUT_PATH$NOW"
fi

echo "output path is $OUTPUT_PATH"

if [ ! -d $OUTPUT_PATH ]; then
mkdir $OUTPUT_PATH
else 
echo "THE FILE ALREADY EXISTS"
exit 0
fi

### description 
touch $OUTPUT_PATH/memo.txt
echo "" >> $OUTPUT_PATH/memo.txt

##### ファイルのコピー
cp $SRC_FILE $OUTPUT_PATH/solve.cpp

########
TOTAL_SCORE=0
for i in {0..49}; do
    printf "#TEST: %03d" $i
    INPUT_FILE="../in/$(printf "%04d" $i)"
    OUTPUT_FILE="$OUTPUT_PATH/out$(printf "%04d" $i).txt"
    ./judge $EXE_FILE < ../in/$(printf "%04d" $i).txt 1> $OUTPUT_FILE 2> /dev/null
    SCORE=$(head -n 1 $OUTPUT_FILE)
    printf ", SCORE: %12d\n" $SCORE
    TOTAL_SCORE=$(($TOTAL_SCORE + $SCORE))
done

printf "########## TOTAL SCORE ##########\n\n"

printf "score: %12d\n\n" $TOTAL_SCORE
# cat $SCORE_RESULT | awk '{sum+=$2} END {printf "score: %10d\n\n" sum}'

printf "#################################\n"
