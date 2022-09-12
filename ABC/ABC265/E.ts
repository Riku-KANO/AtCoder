function main(input) {
  input = input.split("\n");
  let nums = input[0].split(" ").map(val => parseInt(val, 10));
  let n = nums[0];
  let m = nums[1];
  let d = input[1].split(" ").map(val => parseInt(val, 10));
  const MOD = 998244353;
  let xy = new Array(m);
  let mp = new Map();
  for(let i = 0; i < m; i++) {
    xy[i] = input[i+2].split(" ").map(val=>parseInt(val, 10));
    mp[xy[i]] = 1;
  }
  let dp = new Array(n+1);
  for(let i = 0; i <= n; i ++) {
    dp[i] = new Array(n+1);
    for(let j = 0; j <= n; j++) {
      dp[i][j] = new Array(n + 1);
      for(let k = 0; k <= n; k++) dp[i][j][k] = 0;
    }
  }
  dp[0][0][0] = 1;
  for(let num = 1; num <= n; num++) {
    for(let i = 0; i < num; i++) {
      for(let j = 0; j + i < num; j++) {
        let k = num-i-j-1;
        if(k<0) break;
        for(let b = 0; b < 3; b++) {
          let cx = i * d[0] + j * d[2] + k * d[4];
          let cy = i * d[1] + j * d[3] + k * d[5];
          let nx = cx + d[b*2];
          let ny = cy + d[b*2+1];
          if(mp[[nx, ny]] === undefined) {
            if(b==0) {
              dp[i+1][j][k] += dp[i][j][k];
              dp[i+1][j][k] %= MOD;
            }
            if(b==1) {
              dp[i][j+1][k] += dp[i][j][k];
              dp[i][j+1][k] %= MOD;
            }
            if(b==2) {
              dp[i][j][k+1] += dp[i][j][k];
              dp[i][j][k+1] %= MOD;
            }   
          }
         
          if(mp[[nx, ny]] == 1) {
              if(b == 0) dp[i+1][j][k] = 0;
              else if(b == 1) dp[i][j+1][k] = 0;
              else dp[i][j][k+1] = 0;
              continue; 
          }
        }
      }
    }
  }
  let ans = 0;
  for(let i = 0; i <= n; i++) {
    for(let j = 0; j <= n; j++) {
      if(dp[i][j][n-i-j] === undefined) continue;
      ans += dp[i][j][n-i-j];
      ans %= MOD;
    }
  }
  console.log(ans);
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));