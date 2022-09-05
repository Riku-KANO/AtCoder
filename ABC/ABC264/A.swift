let nums = readLine()!.split(separator: " ").map { Int($0)! - 1 }
let l = nums[0]
let r = nums[1]
let s = Array("atcoder")
print(String(s[l...r]))
