// simple BIT
// implemented for i32 and i64
// 1-indexed
// sum method returns closed interval summation
#[allow(dead_code)]
struct BIT<T> {
    n: usize,
    bit: Vec<T>
}

// 1-indexed
#[allow(dead_code)]
impl BIT<i32> {
    fn new(n: usize) -> Self {
        let bit: Vec<i32> = vec![0; n + 1];
        return Self{n: n + 1, bit: bit};
    }

    fn init(&mut self, n: usize) -> () {
        self.n = n + 1;
        self.bit = vec![0; n + 1];
    }

    fn add(&mut self, pos: usize, x: i32) -> () {
        let mut idx: usize = pos;
        let mut idx2: i32 = pos as i32;
        loop {
            if idx >= self.n {
                break;
            }
            self.bit[idx] += x;
            idx2 += idx2 & (-idx2);
            idx = idx2 as usize;
        }
    }

    // return v[1] + v[2] + ... + v[pos]
    fn sum(&mut self, pos: usize) -> i32 {
        let mut sum = 0;
        let mut idx: usize = pos;
        let mut idx2: i32 = pos as i32;
        loop {
            if idx2 <= 0 {
                break;
            }
            sum += self.bit[idx];
            idx2 -= idx2 & (-idx2);
            idx = idx2 as usize;
        }
        return sum;
    }
}

#[allow(dead_code)]
impl BIT<i64> {
    fn new(n: usize) -> Self {
        let bit: Vec<i64> = vec![0; n + 1];
        return Self{n: n + 1, bit: bit};
    }

    fn init(&mut self, n: usize) -> () {
        self.n = n + 1;
        self.bit = vec![0; n + 1];
    }

    fn add(&mut self, pos: usize, x: i64) -> () {
        let mut idx: usize = pos;
        let mut idx2: i64 = pos as i64;
        loop {
            if idx >= self.n {
                break;
            }
            self.bit[idx] += x;
            idx2 += idx2 & (-idx2);
            idx = idx2 as usize;
        }
    }

    // return v[1] + v[2] + ... + v[pos]
    fn sum(&mut self, pos: usize) -> i64 {
        
        let mut sum = 0;
        let mut idx: usize = pos;
        let mut idx2: i64 = pos as i64;
        loop {
            if idx2 <= 0 {
                break;
            }
            sum += self.bit[idx];
            idx2 -= idx2 & (-idx2);
            idx = idx2 as usize;
        }
        return sum;
    }
}