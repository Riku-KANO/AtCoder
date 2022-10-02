OUTPUT_PATH="../input/sample/"

for i in {15..30}; do
NUM=$(($i*2+1))
DIR_NAME="N$NUM"
mkdir $OUTPUT_PATH$DIR_NAME
python3 test_case_generator.py -n $NUM -p $OUTPUT_PATH$DIR_NAME/ --size 5
done