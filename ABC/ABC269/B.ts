function main(input){
  let s = input.split("\n");
  let a = 10, b = 1, c = 10, d = 1;
  for(let i = 0; i < 10; i++) {
    for(let j = 0; j < 10; j++) {
      if(s[i][j] == '#') {
        a = Math.min(i + 1, a);
        b = Math.max(i + 1, b);
        c = Math.min(j + 1, c);
        d = Math.max(j + 1, d);
      }
    }
  }
  console.log(a, b);
  console.log(c, d);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));