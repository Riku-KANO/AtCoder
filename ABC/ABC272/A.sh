tail -1 | awk 'a=0;{for(i=1;i<=NF;i++){a+=$i};print a}'