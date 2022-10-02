SRC_FILE="./solve.cpp"
EXE_FILE_NAME="solve"
EXE_FILE="./solve"
g++ $SRC_FILE -o $EXE_FILE_NAME
OUTPUT_PATH="../output/"
INPUT_PATH="../input/test"
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
echo "baseline" >> $OUTPUT_PATH/memo.txt

##### ファイルのコピー
cp $SRC_FILE $OUTPUT_PATH/solve.cpp

########
TOTAL_SCORE=0

for i in {15..30}; do
NUM=$(($i * 2 + 1))
echo "*************** TEST CASE: N = $NUM *****************"
mkdir $OUTPUT_PATH/N$NUM
OUTPUT_PATH_SECOND=$OUTPUT_PATH/N$NUM
INPUT_PATH_SECOND=$INPUT_PATH/N$NUM
touch $OUTPUT_PATH_SECOND/scores.txt
    for j in {0..39}; do
        printf "#CASE: %03d" $j
        INPUT_FILE="$INPUT_PATH_SECOND/in$(printf "%04d" $j).txt"
        OUTPUT_FILE="$OUTPUT_PATH_SECOND/out$(printf "%04d" $j).txt"
        $EXE_FILE < $INPUT_FILE 1> $OUTPUT_FILE 2> ./test_log.txt
        SCORE=$(cat ./test_log.txt | awk '/^SCORE:/ {print $2}')
        printf ", SCORE:\t%9d\n" $SCORE
        printf "CASE %03d: $SCORE\n" >> $OUTPUT_PATH_SECOND/scores.txt
        TOTAL_SCORE=$(($TOTAL_SCORE + $SCORE))
    done
done

printf "########## TOTAL SCORE ##########\n\n"

printf "score: %12d\n\n" $TOTAL_SCORE

printf "#################################\n"

echo "TOTAL SCORE: $TOTAL_SCORE" >> $OUTPUT_PATH/memo.txt