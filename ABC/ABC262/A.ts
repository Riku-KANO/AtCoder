function main(input){
    var y =  Number(input);
    while(1) {
        if(y % 4 == 2) {
            console.log(y);
            break;
        }
        y++;
    }
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));