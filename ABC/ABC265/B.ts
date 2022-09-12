function main(input) {
  input = input.split("\n");
  let nums = input[0];
  let n = nums[0];
  let m = nums[1];
  let t = nums[2];
  let a = input[1].split(" ").map(val=>parseInt(val, 10));
  let bonus = new Array(m+1);
  bonus.fill(0);
  for(let i = 0; i < m; i++) {
    let xy = input[i+2].split(" ").map(val=>parseInt(val, 10));
    bonus[xy[0]] = xy[1];
  }
  let room = new Array(n + 1);
  for(let i = 0; i < n; i++) {
    t -= a[i];
    if(t <= 0) {
      console.log("No");
      return;
    }
    t += bonus[i+1];
  }
  console.log("Yes");
}

main(require("fs").readFileSync("/dev/stdin", "utf8"));