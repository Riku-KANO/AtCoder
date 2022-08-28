function main(input){
  let s: string = input
  let n: number = s.length;
  console.log(s[Math.floor((n-1)/2)]); 
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));