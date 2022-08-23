function main(input){
    const args = input.split("\n");
    const N: Number = parseInt(args[0], 10);
    const A: Number[] = args[1].split(" ").map( (m) => parseInt(m, 10) );
    var ans: Number = 0;
    var same: Number = 0;
    for(var i:Number = 0; i < N; ++i) {
        if(A[i] === i + 1) same ++;
        else if(A[A[i]-1] === i + 1) ans++;
    }
    ans /= 2;
    ans += same * (same - 1) / 2;
    console.log(ans);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));