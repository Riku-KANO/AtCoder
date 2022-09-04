function main(input){
  input = input.split("\n");
  let nums = input[0].split(" ");
  let n: number = parseInt(nums[0], 10);
  let m: number = parseInt(nums[1], 10);
  let a: number[] = input[1].split(" ").map((val) => parseInt(val, 10));
  let cum: number[] = new Array(n+1);
  cum[0]=0;
  for(let i = 0; i < n; i++) {
    cum[i + 1] = cum[i] + a[i];
  }
  let ans: number = -1e20;
  let sum: number = 0;
  for(let i = 0; i < m; i++) {
    sum += (i+1) * a[i];
  }
  ans = sum;
  for(let i = 0; i + m < n; i++) {
    sum += a[i+m]*m-(cum[m+i]-cum[i]);
    ans = Math.max(ans, sum);
  }
  console.log(ans);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));