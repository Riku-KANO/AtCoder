SRC_FILE="./solve.cpp"
EXE_FILE_NAME="solve3"
EXE_FILE="./solve3"
g++ -std=c++17 $SRC_FILE -o $EXE_FILE_NAME
g++ calc_score.cpp -o calc_score
OUTPUT_PATH="../out/"
INPUT_PATH="../in/"
SCORE_RESULT="$OUTPUT_PATH/result.txt"

if [ "$OUTPUT_PATH" = "../out/" ]; then
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
echo "baseline" >> $OUTPUT_PATH/memo.txt

##### ファイルのコピー
cp $SRC_FILE $OUTPUT_PATH/solve.cpp

########
TOTAL_SCORE=0

for i in {0..19}; do
    echo "*************** TEST CASE: N = $(pritf %04d $i) *****************"
    mkdir $OUTPUT_PATH/
    touch $OUTPUT_PATH_SECOND/scores.txt
    INPUT_FILE="$INPUT_PATH$(printf "%04d" $i).txt"
    OUTPUT_FILE="$OUTPUT_PATH/out$(printf "%04d" $i).txt"
    $EXE_FILE < $INPUT_FILE 1> $OUTPUT_FILE 2> ./test_log.txt
    cat $INPUT_FILE $OUTPUT_FILE > ./calc_input.txt 
    SCORE=$(./calc_score < ./calc_input.txt)
    printf ", SCORE:\t%13d\n" $SCORE
    TOTAL_SCORE=$(($TOTAL_SCORE + $SCORE))
    
done

printf "########## TOTAL SCORE ##########\n\n"

printf "score: %13d\n\n" $TOTAL_SCORE

printf "#################################\n"

echo "TOTAL SCORE: $TOTAL_SCORE" >> $OUTPUT_PATH/memo.txt