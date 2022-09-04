function main(input){
  let s: string = input.trim();
  var ans: number = -1;
  if(s === 'Monday') {ans = 5;}
  else if(s === 'Tuesday') {ans = 4;}
  else if(s === 'Wednesday') {ans = 3;}
  else if(s === 'Thursday') {ans = 2;}
  else if(s === 'Friday') {ans = 1;}
  console.log(ans);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));