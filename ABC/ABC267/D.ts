function main(input){
  input = input.split("\n");
  let nums = input[0].split(" ");
  let n = parseInt(nums[0], 10);
  let m = parseInt(nums[1], 10);
  let a = input[1].split(" ").map(val => parseInt(val, 10));
  const INF = -1e20;
  let dp = new Array(2006);
  for(let i = 0; i < 2006; i++) {
    dp[i] = new Array(2006).fill(INF);
  }
  dp[0][0] = 0;
  for(let i = 0; i < n; i++) {
    for(let j = 0; j <= i; j++) {
      dp[i+1][j] = Math.max(dp[i+1][j], dp[i][j]);
      dp[i+1][j+1] = Math.max(dp[i+1][j+1], dp[i][j] + (j + 1) * a[i]);
    }
  }
  console.log(dp[n][m]);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));